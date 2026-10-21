function frf = estimate_dut_frf( ...
    z, x, fs_z, fs_x, fl, fh, use_window, max_cycles, ...
    decim_coarse, fine_win, loc_win, use_single)
% estimate_dut_frf
%   非同期クロック前提で DUT の周波数応答 H(f) を高速推定。
%   - 各周期ごとに先頭検出（FFT相関：粗→微）
%   - 各周期を線形タイムワープで N 点に整列
%   - DC除去＋同一窓（任意）＋正規化
%   - H1 推定（Sxy/Sxx）で平均化
%
% 入力:
%   z           : 受信信号（DUT出力, 列ベクトル）※複数周期分を含む
%   x           : 送信 1周期チャープ（列ベクトル）
%   fs_z        : z のサンプリング周波数 [Hz]        (既定: 2.5e9)
%   fs_x        : x のサンプリング周波数 [Hz]        (既定: 2.5e9)
%   fl          : 解析下限 [Hz]                       (既定: 100e3)
%   fh          : 解析上限 [Hz]                       (既定: 0.98e9)
%   use_window  : 窓適用（Tukey 2%）                  (既定: true)
%   max_cycles  : 平均に使う最大周期数 K             (既定: 4)
%   decim_coarse: 粗相関のデシメーション率 D         (既定: 32)
%   fine_win    : 微調整の相関窓半幅 [samples]        (既定: 200000)
%   loc_win     : 2周期目以降の局所探索半幅 [samples] (既定: 5000)
%   use_single  : 計算を単精度化（メモリ削減＆高速化）(既定: false)
%
% 出力:
%   frf : struct
%     .f, .H, .mag, .mag_dB, .phase
%     .fl, .fh, .fs_x
%     .num_cycles_used
%     .s0_list         （検出した各周期先頭インデックス）
%     .notes           （処理ノート）

  % ===== 既定値 =====
  if ~exist('fs_z','var') || isempty(fs_z), fs_z = 2.5e9; end
  if ~exist('fs_x','var') || isempty(fs_x), fs_x = 2.5e9; end
  if ~exist('fl','var')   || isempty(fl),   fl   = 100e3; end
  if ~exist('fh','var')   || isempty(fh),   fh   = 0.98e9; end
  if ~exist('use_window','var') || isempty(use_window), use_window = true; end
  if ~exist('max_cycles','var') || isempty(max_cycles), max_cycles = 4; end
  if ~exist('decim_coarse','var') || isempty(decim_coarse), decim_coarse = 32; end
  if ~exist('fine_win','var') || isempty(fine_win), fine_win = 200000; end
  if ~exist('loc_win','var') || isempty(loc_win), loc_win = 5000; end
  if ~exist('use_single','var') || isempty(use_single), use_single = false; end

  % ===== 型・前処理 =====
  z = z(:); x = x(:);
  if use_single
    z = single(z); x = single(x);
  end

  N = length(x);
  if N < 2, error('x length is too short.'); end
  if fs_x <= 0 || fs_z <= 0, error('fs_x, fs_z must be > 0'); end
  if ~(fl < fh), error('Require fl < fh.'); end

  % 1周期の理論長
  T  = N / fs_x;
  Nz = max(1, round(T * fs_z));  % 受信1周期のサンプル数（丸め）

  notes = {};
  notes{end+1} = sprintf('async-fast: T=%.3fus, N=%d, Nz≈%d, fs_x=%.3g, fs_z=%.3g', ...
                         T*1e6, N, Nz, fs_x, fs_z);

  % ===== x を z 時間軸へテンプレ化（相関用）=====
  t_x = (0:N-1).' / fs_x;
  tz  = (0:Nz-1).' / fs_z;
  x0  = x - mean(x);                 % DC除去
  x0  = x0 / max(1e-12, max(abs(x0)));
  xz  = interp1(t_x, x0, tz, 'linear', 0);
  if use_single, xz = single(xz); end

  % ===== 周期先頭検出：粗（デシメーション）→ 微（近傍のみ等倍）=====
  s0_list = detect_period_starts_async(z, xz, Nz, decim_coarse, fine_win, loc_win);
  if isempty(s0_list)
    error('Failed to detect any full period. Increase capture length or adjust params.');
  end

  % ===== 送信用 x の前処理（窓・正規化）=====
  if use_window
    wN = local_tukeywin(N, 0.02);
  else
    wN = ones(N,1);
  end
  if use_single, wN = single(wN); end

  x_use = x - mean(x);               % DC除去
  x_use = x_use .* wN;
  x_use = x_use / max(1e-12, max(abs(x_use)));
  X = fft(x_use, N);                 % 送信のFFT（固定）

  %%% PATCH: 低域の扱いを強化（分解能に合わせて丸め上げ & 低域ビンを積算に寄与させない）
  df = fs_x / N;
  fl_eff = ceil(fl / df) * df;       % 下限周波数をビン境界に丸め上げ
  kmin = ceil(fl_eff / df);          % 有効開始ビン番号（0Hz=0番とした場合の整数）
  % 安全性のためのクリップ
  if kmin < 1, kmin = 1; end
  if kmin > N-1, kmin = N-1; end
  % 低域ビンを0化（以後の Sxy, Sxx に寄与させない）
  X(1:kmin) = 0;

  % ===== 周期ごとのH1積算 =====
  K_available = floor(min( numel(s0_list), max_cycles ));
  Sxy = zeros(N,1, class(X));
  Sxx = zeros(N,1, class(X));

  used = 0;
  for k = 1:K_available
    s0 = s0_list(k);
    if s0 + Nz - 1 > length(z), break; end
    % 該当周期
    yk_z = z(s0 : s0 + Nz - 1);

    % 線形タイムワープで N 点へ（開始0と終了Tを合わせる）
    yk = warp_period_to_N(yk_z, fs_z, N, fs_x);
    % DC・窓・正規化
    yk = yk - mean(yk);
    yk = yk .* wN;
    yk = yk / max(1e-12, max(abs(yk)));

    Yk = fft(yk, N);
    %%% PATCH: 受信スペクトル側も低域を寄与させたくない場合は以下を有効化可
    % Yk(1:kmin) = 0;

    Sxy = Sxy + (Yk .* conj(X));
    Sxx = Sxx + (X  .* conj(X));
    used = used + 1;
  end

  if used < 1
    error('No full usable period found within z.');
  end

  H = Sxy ./ (Sxx + eps_like(Sxx));

  % ===== 片側スペクトル =====
  M = floor(N/2) + 1;
  f_full = (0:N-1).' * df;
  f = (0:M-1).' * df;
  H1 = H(1:M);

  %%% PATCH: 返却前に帯域外を“配列ごと削除”（NaN化しない）
  fl_eff = ceil(fl / df) * df;         % 念のため再計算（可視化・返却用）
  mask = (f >= fl_eff) & (f <= fh) & (f > 0);
  f_sel   = f(mask);
  H1_sel  = H1(mask);

  mag_sel    = abs(H1_sel);
  mag_dB_sel = 20*log10(mag_sel);
  phase_sel  = unwrap(angle(H1_sel));

  % 結果格納（帯域外は完全に排除される）
  frf.f      = double(f_sel);
  frf.H      = H1_sel;
  frf.mag    = mag_sel;
  frf.mag_dB = mag_dB_sel;
  frf.phase  = phase_sel;

  frf.fl     = fl_eff;   % 丸め上げ後の実際の下限
  frf.fh     = fh;
  frf.fs_x   = fs_x;

  frf.num_cycles_used = used;
  frf.s0_list = s0_list(1:used);
  notes{end+1} = sprintf('cycles used: %d / detected: %d', used, numel(s0_list));
  notes{end+1} = sprintf('df=%.3f Hz, fl_req=%.3f Hz, fl_eff=%.3f Hz', df, fl, fl_eff);
  frf.notes  = notes;
end

% ====== 内部関数：粗→微の2段相関で各周期先頭を検出 ======
function s0_list = detect_period_starts_async(z, xz, Nz, D, fine_win, loc_win)
  % 粗相関（D倍デシメーション）
  z_d  = z(1:D:end);
  xt_d = xz(1:D:end);
  Lc   = 2^nextpow2(length(z_d) + length(xt_d) - 1);

  Rcoarse = ifft( fft(z_d, Lc) .* conj(fft(xt_d, Lc)) );
  Rcoarse = real(Rcoarse(1 + length(xt_d) - 1 : length(z_d)));
  [~, k1] = max(abs(Rcoarse));
  s_first_coarse = (k1 - 1) * D + 1;

  % 微調整（先頭近傍のみ等倍）
  a = max(1, s_first_coarse - fine_win);
  b = min(length(z), s_first_coarse + fine_win + Nz - 1);
  if b - a + 1 < Nz
    s0_list = [];
    return;
  end

  L2 = 2^nextpow2((b-a+1) + Nz - 1);
  Rfine = ifft( fft(z(a:b), L2) .* conj(fft(xz, L2)) );
  Rfine = real(Rfine(1 + Nz - 1 : (b - a + 1)));
  [~, kf] = max(abs(Rfine));
  s0_1 = a + kf - 1;

  % 以降は Nz 間隔の予測をベースに局所探索
  s0_list = zeros( min(64, floor(length(z)/Nz) + 2), 1); % 適当上限
  s0_list(1) = s0_1;
  count = 1;

  while true
    s_pred = s0_list(count) + Nz;
    if s_pred + Nz - 1 > length(z), break; end

    a = max(1, s_pred - loc_win);
    b = min(length(z), s_pred + loc_win + Nz - 1);
    if b - a + 1 < Nz, break; end

    L3 = 2^nextpow2((b-a+1) + Nz - 1);
    Rloc = ifft( fft(z(a:b), L3) .* conj(fft(xz, L3)) );
    Rloc = real(Rloc(1 + Nz - 1 : (b - a + 1)));
    [~, kl] = max(abs(Rloc));
    s0_new = a + kl - 1;

    count = count + 1;
    s0_list(count) = s0_new;

    % 安全停止（無限ループ回避）
    if count >= numel(s0_list)
      s0_list(count+1) = 0;
      break;
    end
  end

  s0_list = s0_list(s0_list>0);
end

% ====== 内部関数：受信1周期を線形タイムワープで N 点へ ======
function yN = warp_period_to_N(yk_z, fs_z, N, fs_x)
  Nz = length(yk_z);
  tz = (0:Nz-1).' / fs_z;
  tN = (0:N-1).'  / fs_x;

  % 線形マップ: tz_end -> tN_end を一致させる
  a = ((N-1)/fs_x) / ((Nz-1)/fs_z);
  t_map = a * tz;  % 0→0, end→end に線形伸縮
  yN = interp1(t_map, yk_z, tN, 'linear', 'extrap');
end

% ====== 内部関数：eps（型に合わせる） ======
function e = eps_like(X)
  if isa(X, 'single'), e = single(1e-7); else, e = 1e-14; end
end
