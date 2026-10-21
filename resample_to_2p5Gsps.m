function [waveform_25g, info] = resample_to_2p5Gsps(waveform, actual_params, varargin)
%RESAMPLE_TO_2P5GSPS  LeCroy取得の実数波形を2.5 Gsample/sへリサンプリングする関数．
%
% [waveform_25g, info] = resample_to_2p5Gsps(waveform, actual_params, ...)
%
% 入力:
%   waveform.X [s]        時間軸ベクトル（取得データと同じ長さ）
%   waveform.Y [V]        電圧ベクトル（実数）
%   actual_params         取得パラメータ構造体（SamplingRate, RecordLength, T0 など）
%
% 主要オプション（Name-Value）:
%   'TargetFs' (既定 2.5e9)      目標サンプリング周波数[Hz]
%   'TStartPolicy' ('keep')      出力時間軸の開始を 'keep'(入力X(1)を踏襲) または 'zero'（0開始）にする
%   'FilterOrder' (既定 80)      resample用FIRフィルタの半分タップ数*2（偶数推奨）
%   'UniformTol'  (既定 1e-6)    Xの等間隔チェック相対許容誤差
%
% 出力:
%   waveform_25g.X [s]    2.5 Gsps等間隔の時間軸
%   waveform_25g.Y [V]    リサンプリング済み電圧
%   info                  診断情報（Fs_in, Fs_out, p, q など）
%
% 依存パッケージ:
%   signal  （pkg load signal）
%
% 注意:
%   ・resampleはFIRの群遅延補償を内部で行う．
%   ・入力Xが十分に等間隔でない場合は，自動で等間隔化→resampleに切替える．
%

  % 既定値
  p = inputParser();
  p.addParameter('TargetFs', 2.5e9, @(x) isnumeric(x) && isscalar(x) && x > 0);
  p.addParameter('TStartPolicy', 'keep', @(s) ischar(s) || isstring(s));
  p.addParameter('FilterOrder', 80, @(x) isnumeric(x) && isscalar(x) && x >= 10 && mod(x,2)==0);
  p.addParameter('UniformTol', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
  p.parse(varargin{:});
  TargetFs    = p.Results.TargetFs;
  TStartPolicy= lower(char(p.Results.TStartPolicy));
  Nfir        = p.Results.FilterOrder;
  uniformTol  = p.Results.UniformTol;

  % パッケージ
  try
    pkg("load", "signal");
  catch
    error("signalパッケージがロードできません．`pkg install -forge signal; pkg load signal` を実行してください．");
  end

  % 入力の取り出し
  t_in = waveform.X(:);
  y_in = waveform.Y(:);

  if numel(t_in) ~= numel(y_in)
    error("waveform.X と waveform.Y の長さが一致しません．");
  end
  N_in = numel(y_in);

  % 入力Fsの決定（actual_params優先，無ければXから推定）
  if isfield(actual_params, 'SamplingRate') && ~isempty(actual_params.SamplingRate)
    Fs_in = actual_params.SamplingRate;
  else
    dt = median(diff(t_in));
    Fs_in = 1/dt;
  end

  % 目標Fsが既に同一ならパススルー
  if abs(Fs_in - TargetFs) <= 1e-12*TargetFs
    waveform_25g = waveform;
    info = struct('Fs_in', Fs_in, 'Fs_out', Fs_in, 'ratio', 1, ...
                  'p', 1, 'q', 1, 'method', 'passthrough', ...
                  'N_in', N_in, 'N_out', N_in, 'note', '既に目標Fs');
    return;
  end

  % Xが等間隔かチェック
  dt_vec = diff(t_in);
  dt_med = median(dt_vec);
  uniformity = max(abs(dt_vec - dt_med)) / dt_med;

  % 等間隔でなければ等間隔化
  used_preinterp = false;
  if uniformity > uniformTol
    % 入力を等間隔グリッドに補間（アンチエイリアス用途ではなく，resampleの前処理）
    % グリッドは median(dt) で生成
    t0_align = t_in(1);
    t_in_reg = (t0_align : dt_med : t_in(end)).';
    % 補間法は保形性の高い'makima'（Octaveではinterp1のメソッドに対応）
    y_in_reg = interp1(t_in, y_in, t_in_reg, 'makima', 'extrap');
    t_in = t_in_reg;
    y_in = y_in_reg;
    N_in = numel(y_in);
    Fs_in = 1/dt_med; % 再定義
    used_preinterp = true;
  end

  % レート比の有理数近似（整数p/q）
  % WavePro 16/8/4 Gsps → 2.5 Gsps は厳密に 5/32, 5/16, 5/8 になる．
  [p_up, q_down] = rat(TargetFs / Fs_in, 1e-12);
  % 極端な巨大整数を避けるための安全策（通常は小さい値になる）
  maxInt = 1e6;
  if abs(p_up) > maxInt || abs(q_down) > maxInt
    % 近似精度を少し緩める
    [p_up, q_down] = rat(TargetFs / Fs_in, 1e-6);
  end

  % resample 実行（FIR長はNfir+1程度になる．Octaveのresampleは群遅延補償済）
  % フィルタ長を指定したい場合は 'resample' 第4引数でnを与えられる（Octave: n は半分の次数）
  y_out = resample(y_in, p_up, q_down, Nfir);

  % 出力サンプル数とFs_out
  N_out = numel(y_out);
  Fs_out = (p_up / q_down) * Fs_in;
  Ts_out = 1 / Fs_out;

  % 出力時間軸の構成
  switch TStartPolicy
    case 'keep'
      % 入力の開始時刻を踏襲（群遅延補償済みのため，単純に線形グリッドを置く）
      t0_out = t_in(1);
    case 'zero'
      t0_out = 0;
    otherwise
      error("TStartPolicy は 'keep' または 'zero' を指定してください．");
  end
  t_out = t0_out + (0:N_out-1).' * Ts_out;

  % 目標2.5 Gspsに揃えるために，丸め誤差を微調整
  % 2.5e9に非常に近い場合はFs_outを名目で上書き（情報上）．
  if abs(Fs_out - TargetFs) <= 1e-9*TargetFs
    Fs_out = TargetFs;
    Ts_out = 1/Fs_out;
    % t_outの再生成（端数誤差回避）
    t_out = t0_out + (0:N_out-1).' * Ts_out;
  end

  % 出力構造体
  waveform_25g = struct();
  waveform_25g.X = t_out;
  waveform_25g.Y = y_out;

  % 診断情報
  info = struct();
  info.Fs_in        = Fs_in;
  info.Fs_out       = Fs_out;
  info.TargetFs     = TargetFs;
  info.p            = p_up;
  info.q            = q_down;
  info.ratio        = (p_up / q_down);
  info.used_preinterp = used_preinterp;
  info.N_in         = N_in;
  info.N_out        = N_out;
  info.TStartPolicy = TStartPolicy;
  info.FilterOrder  = Nfir;
  info.Uniformity   = uniformity;
end

