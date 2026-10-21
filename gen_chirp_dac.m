function x = gen_chirp_dac(fl, fh, fs, N, use_window, amp)
% gen_chirp_dac
%   DAC 用の線形 FM（周波数チャープ）波形を生成し、振幅ベクトル x を返します。
%   既定値（省略可）:
%     fl         - 最低周波数 [Hz]    (既定: 100e3)
%     fh         - 最高周波数 [Hz]    (既定: 0.98e9)
%     fs         - サンプリング周波数 [Hz] (既定: 2.5e9)
%     N          - サンプル数（メモリ長） (既定: 2^20)
%     use_window - 両端のクリック抑制用 Tukey 窓を使うか (既定: true)
%     amp        - クリップ回避のための最大振幅 (0<amp<=1, 既定: 0.95)
%
% 出力:
%   x - 1xN の時系列振幅ベクトル（-amp〜+amp の範囲に正規化）
%
% 例:
%   x = gen_chirp_dac();                % 既定値（100 kHz → 1.15 GHz）
%   x = gen_chirp_dac(1e6, 1.0e9);      % 最低1 MHz, 最高1.0 GHz（他は既定）
%   x = gen_chirp_dac([], [], 2.5e9, 2^18, false, 0.9);
%
% 備考:
%  - fh は fs/2 を超えないように自動でクリップします。
%  - fl >= 0、fl < fh を満たす必要があります。
%  - 必要に応じて int16 などへ量子化して DAC に渡してください。
%
% 作者メモ: Kazuhiko さんの測定条件に合わせた初期値を設定

  % 既定値（空や未指定なら適用）
  if ~exist('fl','var') || isempty(fl), fl = 100e3; end
  if ~exist('fh','var') || isempty(fh), fh = 0.98e9; end
  if ~exist('fs','var') || isempty(fs), fs = 2.5e9; end
  if ~exist('N','var')  || isempty(N),  N  = 2^20;  end
  if ~exist('use_window','var') || isempty(use_window), use_window = true; end
  if ~exist('amp','var') || isempty(amp), amp = 0.95; end

  % 入力チェックと調整
  if fs <= 0, error('fs must be > 0'); end
  if N <= 1,  error('N must be > 1');   end
  nyq = fs/2;
  if fh > nyq
    warning('fh=%.3g Hz exceeds Nyquist (%.3g Hz). Clipping to Nyquist.', fh, nyq);
    fh = nyq;
  end
  if fl < 0
    warning('fl=%.3g Hz < 0. Clipping to 0 Hz.', fl);
    fl = 0;
  end
  if ~(fl < fh)
    error('Invalid band: require fl < fh (got fl=%.3g, fh=%.3g).', fl, fh);
  end
  if ~(amp > 0 && amp <= 1)
    error('amp must be in (0, 1].');
  end

  % 時間軸とチャープ位相
  dt = 1/fs;
  t  = (0:N-1) * dt;
  T  = N / fs;
  k  = (fh - fl) / T;                  % チャープレート [Hz/s]
  phi = 2*pi*( fl*t + 0.5*k.*(t.^2) ); % φ(t) = 2π(fl t + 0.5 k t^2)

  % 生のチャープ
  x = sin(phi);

  % 端のクリック低減用 窓（任意）
  if use_window
    alpha = 0.02; % 2% フェード
    w = local_tukeywin(N, alpha).';
    x = x .* w;
  end

  % クリップ回避のためのスケーリング
  x = amp * (x / max(abs(x)));

  % （必要なら量子化例）
  % x = int16(round(x * 32767));
end

