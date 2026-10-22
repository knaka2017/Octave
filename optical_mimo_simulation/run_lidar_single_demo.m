% =========================================================================
% スクリプト名: run_lidar_single_demo.m
%
% 【概要】
% 光路長差による伝搬時間遅延 (Time of Flight: TOF) を考慮した
% 2次元光空間伝搬シミュレーション関数 optical_mimo_simulation を用いた
% 1光源 (1Tx) パルスLiDARのTOF検証デモスクリプトです。
%
% 【検証ポイント】
%   - fs = 2.5 GHz (0.4 ns 分解能) によるパルス波形の時系列シミュレーション
%   - FWHM 10 ns のガウシアンパルスを、シミュレーション開始から遅らせて発射
%   - ターゲット距離 9m（往復約 18m）による光の到達遅延（約 60 ns）
%   - 受光パルスのピーク時刻から求めた TOF とターゲット距離の推定
%   - 水中伝搬（屈折率 1.33・水質ごとの消散係数）による TOF の増大と受光パルスの減衰
%   - 受光モジュール C5658 相当（検出感度 2.5e5 V/W、φ0.5 mm、0.05 MHz ～ 1 GHz）の
%     固定雑音（光入力換算 16 nW rms）による SNR と検出可否
% =========================================================================

clear;
close all;
clc;

fprintf('\n');
fprintf('  =======================================================\n');
fprintf('    1光源 パルスLiDAR 2次元空間伝搬シミュレーション (TOF遅延モデル対応)\n');
fprintf('  =======================================================\n\n');

%% ========================================================================
% 【ケース 1】ガウシアンパルスによるTOF時系列シミュレーション（fs = 2.5 GHz）
% ========================================================================
fprintf('--- [ケース 1] ガウシアンパルスTOFシミュレーション ---\n');

c0 = 2.99792458e8;              % 真空中の光速 [m/s]
target_distance = 9.0;          % ターゲット正面距離 [m]
n_medium = 1.0;                 % 媒体屈折率 (空気中: 1.0, 純水: 1.33)

% 1. 送信パルス条件
P_peak     = 1.0;               % パルスのピーク光パワー [W]
pulse_fwhm = 10.0e-9;           % パルス半値全幅 FWHM 10 ns
t_fire     = 20.0e-9;           % パルス中心の発射時刻 20 ns (シミュレーション開始から遅延)
sigma_p    = pulse_fwhm / (2 * sqrt(2 * log(2)));   % ガウス標準偏差 (約 4.25 ns)

% 2. 時間軸の定義（往復伝搬時間 + パルス裾野が収まる解析時間を自動設定）
fs = 2.5e9;                     % サンプリング周波数 2.5 GHz (0.4 ns 分解能)
dt = 1 / fs;
t_stop = t_fire + 2 * target_distance * n_medium / c0 + 4 * pulse_fwhm;  % 約 120 ns
t_vec = 0:dt:t_stop;
num_samples = length(t_vec);

% 3. 送信光強度信号 x(t) の定義（1光源: 1 x num_samples）
x_pulse = P_peak * exp(-((t_vec - t_fire).^2) / (2 * sigma_p^2));

% 4. シミュレーションパラメータの設定
params_ts = struct();
params_ts.fs = fs;
params_ts.enable_plot = true;            % 可視化を有効化
params_ts.medium_extinction = 0.0;       % 空気中
params_ts.refractive_index = n_medium;   % 媒体屈折率
params_ts.initial_state = 'zero';        % 単発パルス: t < 0 は無送信・受光回路は静止状態

% --- 送信光源 (Tx) の配置・角度 (y軸上 x=0 に配置, ターゲット面に対し垂直 0 deg) ---
params_ts.tx_array(1).pos          = [0.0,  0.50];
params_ts.tx_array(1).type         = 'LD';
params_ts.tx_array(1).pointing_deg =   0.0;
params_ts.tx_array(1).fwhm         = 25.0;

% --- 受光アレイ (Rx) の配置・受光軸 ---
params_ts.rx_array(1).pos          = [0.0, -0.25];
params_ts.rx_array(1).pointing_deg =   0.0;      % 正面 0 deg

params_ts.rx_array(2).pos          = [0.0, -1.050];
params_ts.rx_array(2).pointing_deg =   0.0;      % 正面 0 deg

% --- 受光器共通仕様 (浜松ホトニクス C5658 相当) ---
params_ts.rx_common.aperture_width = 0.5e-3;    % 受光開口幅 0.5 mm (受光面 φ0.5 mm を2次元開口幅として扱う)
params_ts.rx_common.fov_deg        = 45.0;      % 視野角 半角 45 deg
params_ts.rx_common.responsivity   = 45.0;      % APD感度 45 A/W (量子効率 70% @800 nm, M = 100 から概算)
params_ts.rx_common.tia_gain       = 2.5e5 / params_ts.rx_common.responsivity;  % 検出感度 2.5e5 V/W に合わせた変換ゲイン [V/A]
params_ts.rx_common.hpf_cutoff     = 0.05e6;    % 低域遮断 0.05 MHz (C5658, AC結合)
params_ts.rx_common.lpf_cutoff     = 1.0e9;     % 高域遮断 1 GHz (C5658)
params_ts.rx_common.compression_power = 3.0e-6; % 1 dB 圧縮点 3 uW (超過時は警告のみ)

% --- 受光器雑音 (光入力換算の固定雑音) ---
params_ts.noise_params.enable          = true;
params_ts.noise_params.model           = 'fixed';
params_ts.noise_params.input_noise_rms = 16e-9;  % C5658 検出限界 -48 dBm (16 nW rms)
params_ts.noise_params.seed            = 1;      % 再現性のため乱数シードを固定
det_threshold = 5.0;                              % パルス検出閾値 (出力雑音 rms の倍数)

% --- 反射ターゲット仕様 (距離 9.0m, 長さ 1.8m) ---
params_ts.target.distance    = target_distance; % 正面距離 9.0 m
params_ts.target.length      = 1.8;             % 長さ 1.8 m
params_ts.target.num_patches = 90;              % 90パッチ分割
params_ts.target.reflectance = 0.8;             % 反射率 0.8

% 5. 関数の呼び出し
[y_ts, H_elec_ts, H_opt_ts, info_ts] = optical_mimo_simulation(x_pulse, params_ts);

% 6. チャネル解析情報のコンソール出力
fprintf('  -> 計算完了:\n');
fprintf('     ・静的電気チャネル行列 H_elec (DC極限) [V/W]:\n');
disp(H_elec_ts);
fprintf('     ・特異値 (Singular Values)         : %s\n', ...
sprintf('%.2e ', info_ts.singular_values));
fprintf('     ・行列のランク (Rank)              : %d\n\n', info_ts.rank);

% 7. 受光パルスのピーク時刻による TOF・距離推定
c_medium = c0 / n_medium;
num_rx = size(y_ts, 1);

% 受光回路 1次LPF の群遅延 1/(2 pi f_H)（パルス帯域 << f_H での近似）
if params_ts.rx_common.lpf_cutoff > 0
tau_rx = 1 / (2 * pi * params_ts.rx_common.lpf_cutoff);
else
tau_rx = 0;
end

fprintf('  -> パルスTOF解析 (発射時刻 t_fire = %.2f ns):\n', t_fire * 1e9);
for n = 1:num_rx
dH_n   = squeeze(info_ts.patch_data.dH(n, 1, :));
path_n = squeeze(info_ts.patch_data.path_length(n, 1, :));

if sum(dH_n) <= 0
    fprintf('     ・Rx #%d: 受光なし (FOV外)\n', n);
    continue;
end

% 理論値: パッチ寄与 dH で重み付けした平均往復光路長
path_mean = sum(dH_n .* path_n) / sum(dH_n);
tof_theory = path_mean / c_medium;

fprintf('     ・Rx #%d: 受光光パワーのピーク %.3e W, ピーク SNR %.1f dB\n', ...
    n, info_ts.P_rx_peak(n), info_ts.snr_peak_db(n));

% 雑音を含む出力のピークが検出閾値に届かない場合は検出不可
if max(y_ts(n, :)) < det_threshold * info_ts.noise_rms_v
    fprintf('              検出不可 (ピークが雑音 rms の %.0f 倍未満)\n', det_threshold);
    continue;
end

% 実測値: 雑音を含む受光波形のピーク時刻（3点放物線補間でサブサンプル推定）
[~, idx] = max(y_ts(n, :));
t_peak = t_vec(idx);
if idx > 1 && idx < num_samples
    y_l = y_ts(n, idx - 1);
    y_c = y_ts(n, idx);
    y_r = y_ts(n, idx + 1);
    denom = y_l - 2 * y_c + y_r;
    if denom ~= 0
        t_peak = t_peak + 0.5 * (y_l - y_r) / denom * dt;
    end
end
tof_meas = t_peak - t_fire;

fprintf('              ピーク電圧 %.3e V (雑音 rms %.3e V)\n', max(y_ts(n, :)), info_ts.noise_rms_v);
fprintf('              TOF 理論値 %.3f ns (往復光路長 %.3f m: %.3f ns + 受光回路群遅延 %.3f ns)\n', ...
    (tof_theory + tau_rx) * 1e9, path_mean, tof_theory * 1e9, tau_rx * 1e9);
fprintf('              TOF 実測値 %.3f ns -> 推定距離 (片道換算, 受光回路遅延補正後) %.3f m\n', ...
    tof_meas * 1e9, c_medium * (tof_meas - tau_rx) / 2);


end
fprintf('\n');

%% ========================================================================
% 【ケース 2】静的DC光パワー解析（スナップショット）
% ========================================================================
fprintf('--- [ケース 2] 静的DC光パワー解析 (DCブロック無効) ---\n');

x_static = 1.0;
params_static = params_ts;
params_static.rx_common.hpf_cutoff = 0;        % DC解析のため HPF (AC結合) を無効化
params_static.noise_params.enable = false;      % DC解析のため雑音を無効化
params_static.enable_plot = false;

[y_static, H_elec_dc, ~, info_static] = optical_mimo_simulation(x_static, params_static);

fprintf('  -> 送信光パワー x_static [W]:\n');
disp(x_static');
fprintf('  -> 受光光パワー P_rx [W]:\n');
disp(info_static.P_rx');
fprintf('  -> 受光直流電圧 y_static [V]:\n');
disp(y_static');

%% ========================================================================
% 【ケース 3】水中伝搬シミュレーション（屈折率・水質ごとの消散係数）
%   ケース 1 と同じ配置・送信パルスで媒体だけを変え、Rx #1 の受光パルスを比較する。
%   消散係数 c = a + b [1/m] は波長 532 nm での代表値
%   (Adnan et al., IJDNE 16(2), 2021 / Kumar et al., IJSER 2(6), 2011)。
%   散乱光はすべて損失として扱い、多重散乱・水中後方散乱は含まない。
% ========================================================================
fprintf('--- [ケース 3] 水中伝搬シミュレーション (屈折率・消散係数) ---\n');

media_name = {'空気',  '水 (消散なし)', '清澄な海水', '沿岸海水', '濁った港湾水'};
media_label = {'Air', 'Water (Lossless)', 'Clean Ocean', 'Coastal Ocean', 'Turbid Harbor'};   % グラフ凡例用
media_tick  = {'Air', 'Water', 'Clean', 'Coastal', 'Harbor'};                                  % グラフ軸目盛用
media_n    = [1.0,     1.33,            1.33,         1.33,       1.33];
media_ext  = [0.0,     0.0,             0.15,         0.30,       2.19];   % 消散係数 [1/m]
num_media  = length(media_name);
rx_view    = 1;                 % 比較に用いる受光器番号

% 時間軸: 最も遅い媒体 (屈折率最大) の往復伝搬時間 + パルス裾野が収まるよう設定
t_stop_w = t_fire + 2 * target_distance * max(media_n) / c0 + 4 * pulse_fwhm;  % 約 140 ns
t_vec_w = 0:dt:t_stop_w;
num_samples_w = length(t_vec_w);
x_pulse_w = P_peak * exp(-((t_vec_w - t_fire).^2) / (2 * sigma_p^2));

y_media = zeros(num_media, num_samples_w);    % 雑音を含む受光モジュール出力 [V]
v_sig_media = zeros(num_media, num_samples_w);% 雑音付加前の受光モジュール出力 [V]
prx_peak = zeros(1, num_media);               % 受光光パワーのピーク [W]
snr_peak = zeros(1, num_media);               % ピーク SNR [dB]
is_comp = false(1, num_media);                % 1 dB 圧縮点の超過

noise_in_rms = params_ts.noise_params.input_noise_rms;                      % 光入力換算の雑音 rms [W]
noise_out_rms = noise_in_rms * params_ts.rx_common.responsivity * params_ts.rx_common.tia_gain;  % 出力換算 [V]
fprintf('  -> Rx #%d の受光パルス (ターゲット距離 %.1f m, 発射時刻 %.2f ns):\n', ...
rx_view, target_distance, t_fire * 1e9);
fprintf('     雑音: 光入力換算 %.3g nW rms (出力換算 %.3g mV rms), 検出閾値: 雑音 rms の %.0f 倍 (%.1f dB)\n', ...
noise_in_rms * 1e9, noise_out_rms * 1e3, det_threshold, 20 * log10(det_threshold));
fprintf('      n     c [1/m]  P_rx ピーク [W]  空気比 [dB]  SNR [dB]  TOF [ns]  推定距離 [m]  媒体\n');

% 1 dB 圧縮点の超過は表の「*」で示すため、媒体ごとの警告は一時的に抑制
warn_state = warning('query', 'optical_mimo_simulation:compression');
warning('off', 'optical_mimo_simulation:compression');
for k = 1:num_media
params_w = params_ts;
params_w.enable_plot = false;
params_w.refractive_index = media_n(k);
params_w.medium_extinction = media_ext(k);

[y_w, ~, ~, info_w] = optical_mimo_simulation(x_pulse_w, params_w);
y_media(k, :) = y_w(rx_view, :);
v_sig_media(k, :) = info_w.v_signal(rx_view, :);
prx_peak(k) = info_w.P_rx_peak(rx_view);
snr_peak(k) = info_w.snr_peak_db(rx_view);
is_comp(k) = info_w.is_compressed(rx_view);
noise_rms_w = info_w.noise_rms_v;

comp_mark = ' ';
if is_comp(k)
    comp_mark = '*';
end
row_head = sprintf('     %.2f  %5.2f    %.3e%s      %8.2f   %7.1f', media_n(k), media_ext(k), ...
    prx_peak(k), comp_mark, 10 * log10(prx_peak(k) / prx_peak(1)), snr_peak(k));

% 雑音を含む出力のピークが検出閾値に届かない場合は検出不可
[v_max, idx] = max(y_media(k, :));
if v_max < det_threshold * noise_rms_w
    fprintf('%s     ---      検出不可      %s\n', row_head, media_name{k});
    continue;
end

% 受光波形のピーク時刻（3点放物線補間でサブサンプル推定）
t_peak = t_vec_w(idx);
if idx > 1 && idx < num_samples_w
    y_l = y_media(k, idx - 1);
    y_c = y_media(k, idx);
    y_r = y_media(k, idx + 1);
    denom = y_l - 2 * y_c + y_r;
    if denom ~= 0
        t_peak = t_peak + 0.5 * (y_l - y_r) / denom * dt;
    end
end
tof_meas = t_peak - t_fire;
c_medium_k = c0 / media_n(k);

fprintf('%s   %7.3f   %7.3f       %s\n', row_head, tof_meas * 1e9, ...
    c_medium_k * (tof_meas - tau_rx) / 2, media_name{k});


end
warning(warn_state.state, 'optical_mimo_simulation:compression');
fprintf('     (空気比は受光光パワーの比 10 log10。* は 1 dB 圧縮点 %.3g W 超過で、飽和は模擬していない)\n', ...
params_ts.rx_common.compression_power);
fprintf('     (推定距離は媒体中の光速 c0/n と受光回路群遅延 %.3f ns の補正を用いた片道換算値)\n\n', tau_rx * 1e9);

% 媒体ごとの受光パルス波形（雑音を含む実際の出力電圧、媒体ごとに縦軸を自動調整）
figure('Name', 'Pulse LiDAR - Underwater Received Waveforms', ...
'NumberTitle', 'off', 'Position', [100, 60, 820, 900]);
for k = 1:num_media
subplot(num_media, 1, k);
hold on; grid on;
plot(t_vec_w * 1e9, y_media(k, :) * 1e3, 'Color', [0.2, 0.5, 0.8], 'LineWidth', 1.0, ...
    'DisplayName', 'Output with Noise');
plot(t_vec_w * 1e9, v_sig_media(k, :) * 1e3, 'k--', 'LineWidth', 1.2, ...
    'DisplayName', 'Noise-Free Signal');
plot([0, t_stop_w * 1e9], [1, 1] * det_threshold * noise_out_rms * 1e3, 'r:', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Detection Threshold (%.0f\\sigma = %.3g mV)', det_threshold, det_threshold * noise_out_rms * 1e3));
xlim([0, t_stop_w * 1e9]);
ylabel('v_{out} [mV]');
title(sprintf('%s (n = %.2f, c = %.2f /m): P_{rx} peak %.3g W, SNR %.1f dB', ...
    media_label{k}, media_n(k), media_ext(k), prx_peak(k), snr_peak(k)));
if k == 1
    legend('Location', 'northwest');
end
if k == num_media
    xlabel('Time [ns]');
end
end

% 受光パルスのピーク SNR と検出閾値（赤破線）。雑音レベルをタイトルに明示
figure('Name', 'Pulse LiDAR - Underwater Peak SNR', ...
'NumberTitle', 'off', 'Position', [940, 60, 760, 540]);
hold on;
bar(1:num_media, snr_peak, 'FaceColor', [0.2, 0.5, 0.8]);
grid on;
snr_th_db = 20 * log10(det_threshold);
y_lo = min(min(snr_peak) * 1.25, -10);
y_hi = max(max(snr_peak) * 1.25, snr_th_db * 2);
ylim([y_lo, y_hi]);
plot([0.5, num_media + 0.5], [1, 1] * snr_th_db, 'r--', 'LineWidth', 1.2);
for k = 1:num_media
if snr_peak(k) >= 0
    text(k, snr_peak(k) + 0.03 * (y_hi - y_lo), sprintf('%.1f dB', snr_peak(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
else
    text(k, snr_peak(k) - 0.03 * (y_hi - y_lo), sprintf('%.1f dB', snr_peak(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);
end
end
set(gca, 'XTick', 1:num_media);
set(gca, 'XTickLabel', media_label);
xlim([0.5, num_media + 0.5]);
xlabel('Medium');
ylabel('Peak SNR [dB]');
title({sprintf('Rx #%d Peak SNR = 20 log_{10}(v_{signal,peak} / \\sigma_{noise})', rx_view), ...
sprintf('\\sigma_{noise}: %.3g nW rms Input-Referred (%.3g mV rms at Output), Red Dashed: Threshold %.0f dB (%.0f\\sigma)', ...
noise_in_rms * 1e9, noise_out_rms * 1e3, snr_th_db, det_threshold)});

fprintf('すべてのデモケースが正常に終了しました。\n');
fprintf('====================================================\n');
