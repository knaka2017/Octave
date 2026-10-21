%% VNA-like Frequency Response Analyzer - 関数版
%
% 使用例:
%   [freq, H] = vna_analysis(x, fs_x, z, fs_z);
%   [freq, H] = vna_analysis(x, fs_x, z, fs_z, 'plot', true, 'n_avg', 4);

function [freq, H] = vna_analysis(x, fs_x, z, fs_z, varargin)
% VNA_ANALYSIS  DUT の周波数応答（振幅・位相特性）を計算する
%
% 入力:
%   x      - DUT への入力信号ベクトル（AWG 波形, 長さ N サンプル）
%   fs_x   - x のサンプリング周波数 [Hz]
%   z      - DUT からの出力信号ベクトル（オシロ取得, 長さ M サンプル）
%   fs_z   - z のサンプリング周波数 [Hz]
%
% オプション（名前付き引数）:
%   'n_avg'      - 平均化する周期数（デフォルト: 利用可能な最大数）
%   'threshold'  - X(f) の有効判定閾値（最大値に対する比, デフォルト: 1e-6）
%   'plot'       - グラフ表示フラグ（デフォルト: true）
%   'freq_range' - 表示する周波数範囲 [Hz]（デフォルト: [0, fs_z/2]）
%
% 出力:
%   freq - 周波数ベクトル（片側, DC〜Nyquist）[Hz]
%   H    - 周波数応答ベクトル（複素数, 励振されていない周波数は NaN）

  %% --- 入力チェック ---
  x = x(:);  % 列ベクトルに統一
  z = z(:);

  %% --- オプション引数のパース ---
  p = inputParser();
  addParameter(p, 'n_avg',      [],    @(v) isempty(v) || (isnumeric(v) && v >= 1));
  addParameter(p, 'threshold',  1e-6,  @(v) isnumeric(v) && v > 0);
  addParameter(p, 'plot',       true,  @islogical);
  addParameter(p, 'freq_range', [],    @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
  parse(p, varargin{:});
  opt = p.Results;

  %% --- Step 1: x を fs_z にリサンプルして参照スペクトル X(f) を計算 ---
  [pq_p, pq_q] = rat(fs_z / fs_x, 1e-6);
  x_ref     = resample(x, pq_p, pq_q);        % 1周期分をリサンプル
  N_ref     = length(x_ref);
  X_fft     = fft(x_ref);
  freq_full = (0:N_ref-1)' * fs_z / N_ref;    % 周波数ベクトル（全帯域）

  %% --- Step 2: z から複数周期を切り出してコヒーレント平均化 ---
  M       = length(z);
  n_avail = floor(M / N_ref);
  if n_avail < 1
    error('vna_analysis: z のサンプル数 (%d) が1周期分 (%d samples) に満たない', ...
          M, N_ref);
  end

  if isempty(opt.n_avg)
    n_avg = n_avail;
  else
    n_avg = min(floor(opt.n_avg), n_avail);
  end

  z_matrix = reshape(z(1 : n_avg*N_ref), N_ref, n_avg);
  z_avg    = mean(z_matrix, 2);    % コヒーレント平均（ノイズ低減）

  fprintf('[vna_analysis] N_ref=%d, 平均化周期数=%d / %d\n', N_ref, n_avg, n_avail);

  %% --- Step 3: 出力スペクトル Z(f) の計算 ---
  Z_fft = fft(z_avg);

  %% --- Step 4: H(f) = Z(f) / X(f) ---
  threshold_abs        = max(abs(X_fft)) * opt.threshold;
  valid                = abs(X_fft) > threshold_abs;
  H_full               = complex(NaN(N_ref, 1));
  H_full(valid)        = Z_fft(valid) ./ X_fft(valid);

  %% --- 片側スペクトルに変換（DC 〜 Nyquist）---
  half_N = floor(N_ref/2) + 1;
  freq   = freq_full(1:half_N);
  H      = H_full(1:half_N);

  %% --- Step 5: プロット ---
  if opt.plot
    vna_plot_(freq, H, fs_z, fs_x, N_ref, n_avg, opt.freq_range);
  end

end  % function vna_analysis


%% ============================================================
%  内部関数: プロット
%% ============================================================
function vna_plot_(freq, H, fs_z, fs_x, N_ref, n_avg, freq_range)

  H_mag_dB = 20*log10(abs(H));
  H_phase  = unwrap(angle(H)) * 180/pi;
  valid    = ~isnan(H);

  if isempty(freq_range)
    freq_range = [0, fs_z/2];
  end

  figure('Position', [100, 100, 900, 700], 'Name', 'VNA-like Frequency Response');

  % --- 振幅特性 ---
  subplot(2,1,1);
  plot(freq(valid)/1e6, H_mag_dB(valid), 'ro-', 'MarkerSize', 5, 'LineWidth', 1.2);
  grid on;
  xlabel('周波数 [MHz]');
  ylabel('振幅 [dB]');
  title('振幅特性 |H(f)|');
  xlim(freq_range / 1e6);
  set(gca, 'FontSize', 11);

  % --- 位相特性 ---
  subplot(2,1,2);
  plot(freq(valid)/1e6, H_phase(valid), 'bo-', 'MarkerSize', 5, 'LineWidth', 1.2);
  grid on;
  xlabel('周波数 [MHz]');
  ylabel('位相 [deg]');
  title('位相特性 ∠H(f)');
  xlim(freq_range / 1e6);
  set(gca, 'FontSize', 11);

  % sgtitle は Octave 未対応のため annotation で代替
  axes('Position', [0 0 1 1], 'Visible', 'off');
  text(0.5, 0.98, ...
       sprintf('DUT 周波数応答  (fs_x=%.0f MHz, fs_z=%.0f MHz, N=%d, 平均 %d 周期)', ...
               fs_x/1e6, fs_z/1e6, N_ref, n_avg), ...
       'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
       'FontSize', 13, 'FontWeight', 'bold', 'Units', 'normalized');

end  % function vna_plot_
