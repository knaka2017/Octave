function y = chirp01(N = 2^20, fs = 2.5e9, f1, f2)
  % chirp01: 境界で完全に連続するように調整された正弦波の往復チャープ信号
  %
  % 使用例：
  %   y = chirp01(2^20, 2.5e9, 100e3, 6e3);
  %
  % 入力：
  %   N  - データ長（デフォルト: 2^20）
  %   fs - サンプリング周波数 [Hz]（デフォルト: 2.5GHz）
  %   f1 - 開始周波数 [Hz]（必須, 境界で連続になるように調整）
  %   f2 - 中間周波数 [Hz]（必須, f1 > f2）
  %
  % 出力：
  %   y  - 正弦波のチャープ信号（往復: f1 -> f2 -> f1）

  % 引数チェック
  if nargin < 4
    error('f1 および f2 を指定してください。');
  end
  if f1 < f2
    error('f1 は f2 より大きい値を指定してください。');
  end

  % f1 を調整（x(N) ≈ 0 になるように）
  f1_adj = (round(f1 * N / fs + 0.5) - 0.5) * fs / N;

  % 調整後の f1 を表示
  fprintf('adjusted f1: %.6f Hz -> %.6f Hz\n', f1, f1_adj);

  % 時間ベクトル
  t = (0:N-1) / fs;  % 時間軸 [秒]

  % チャープの折り返し位置（全データ長の半分）
  N_half = floor(N / 2);

  % 周波数変化率（線形変化）
  k1 = (f2 - f1_adj) / (N_half / fs);  % f1_adj -> f2
  k2 = (f1_adj - f2) / ((N - N_half) / fs);  % f2 -> f1_adj

  % 位相計算
  phi = zeros(1, N);
  phi(1:N_half) = 2 * pi * (f1_adj * t(1:N_half) + 0.5 * k1 * t(1:N_half).^2);
  phi(N_half+1:end) = 2 * pi * (f2 * (t(N_half+1:end) - t(N_half)) + 0.5 * k2 * (t(N_half+1:end) - t(N_half)).^2) + phi(N_half);

  % チャープ信号の生成
  y = sin(phi);
end

