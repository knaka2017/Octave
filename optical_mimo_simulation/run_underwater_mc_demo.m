% =========================================================================
% スクリプト名: run_underwater_mc_demo.m
%
% 【概要】
% 3次元モンテカルロ法（局所推定法）による水中パルス LiDAR シミュレーション
% underwater_mc_lidar_simulation を用いて、水質ごとの散乱の影響を比較するデモ。
% （2次元・1回反射の決定論的モデル optical_mimo_simulation とは別モデル）
%
% 【検証ポイント】
%   - 散乱を含めた受光光パワー（ターゲット反射、水の後方散乱、散乱なしの弾道成分）
%   - 前方散乱による受光パワーの増加（Beer-Lambert の弾道成分との比較）
%   - 受光モジュール C5658 相当（2.5e5 V/W、0.05 MHz ～ 1 GHz、16 nW rms）での検出可否
%   - 計算時間・イベント数（計算リソースの見積もり用）
%
% 【主な仮定】
%   - 送信ビーム: 発散角 FWHM 1 deg（ビーム整形した LD を想定）
%   - 受光光学系: 開口径 25.4 mm のレンズ + φ0.5 mm の APD、視野角 半角 1 deg
%     （エテンデュ: 12.7 mm x sin(1 deg) = 0.22 mm <= 受光面半径 0.25 mm）
%   - HG 位相関数の非対称パラメータ g = 0.9（仮定値。結果は g に敏感）
%   - 窓・水面での屈折、偏光、背景光は考慮しない
% =========================================================================

clear;
close all;
clc;

fprintf('\n');
fprintf('  =======================================================\n');
fprintf('    水中パルス LiDAR 3次元モンテカルロ・シミュレーション (散乱を考慮)\n');
fprintf('  =======================================================\n\n');

c0 = 2.99792458e8;
n_water = 1.33;
target_distance = 9.0;

% 1. 送信パルス（run_lidar_single_demo と同じ）
P_peak     = 1.0;                % ピーク光パワー [W]
pulse_fwhm = 10.0e-9;            % FWHM 10 ns
t_fire     = 20.0e-9;            % 発射時刻 20 ns
sigma_p    = pulse_fwhm / (2 * sqrt(2 * log(2)));

fs = 2.5e9;
t_stop = t_fire + 2 * target_distance * n_water / c0 + 4 * pulse_fwhm;   % 約 140 ns
num_samples = round(t_stop * fs);
t_vec = (0:(num_samples - 1)) / fs;
x_pulse = P_peak * exp(-((t_vec - t_fire).^2) / (2 * sigma_p^2));

% 2. 共通パラメータ
params = struct();
params.fs = fs;
params.refractive_index = n_water;
params.enable_plot = false;

params.tx.pos          = [0.0,  0.05, 0.0];   % Tx と Rx の間隔 10 cm
params.tx.pointing_deg = 0.0;
params.tx.fwhm         = 1.0;                  % ビーム発散角 FWHM 1 deg

params.rx_array(1).pos          = [0.0, -0.05, 0.0];
params.rx_array(1).pointing_deg = 0.0;

params.rx_common.aperture_diameter = 25.4e-3;  % 受光レンズ開口径 25.4 mm
params.rx_common.fov_deg           = 1.0;      % 視野角 半角 1 deg
params.rx_common.responsivity      = 45.0;     % APD感度 45 A/W (量子効率 70% @800 nm, M = 100 から概算)
params.rx_common.tia_gain          = 2.5e5 / params.rx_common.responsivity;  % 検出感度 2.5e5 V/W
params.rx_common.hpf_cutoff        = 0.05e6;   % 低域遮断 0.05 MHz (C5658)
params.rx_common.lpf_cutoff        = 1.0e9;    % 高域遮断 1 GHz (C5658)
params.rx_common.compression_power = 3.0e-6;   % 1 dB 圧縮点 3 uW

params.target.distance    = target_distance;
params.target.width       = 1.8;
params.target.height      = 1.8;
params.target.reflectance = 0.8;

params.noise_params.enable          = true;
params.noise_params.input_noise_rms = 16e-9;  % C5658: -48 dBm (16 nW rms)
params.noise_params.seed            = 1;

params.mc.num_photons = 1e6;
params.mc.batch_size  = 2e5;
params.mc.seed        = 1;

det_threshold = 5.0;                          % 検出閾値（出力雑音 rms の倍数）

% 3. 水質（Jasman et al., IJECE 8(4), 2018 の表: Mobley 1994 / Petzold 1972 より）
water_name  = {'清澄な海水', '沿岸海水', '濁った港湾水'};
water_label = {'Clean Ocean', 'Coastal Ocean', 'Turbid Harbor'};
water_a     = [0.114, 0.179, 0.366];          % 吸収係数 [1/m]
water_b     = [0.037, 0.219, 1.824];          % 散乱係数 [1/m]
water_g     = 0.9;                            % HG 非対称パラメータ（仮定値）
num_water   = length(water_name);

% ターゲット反射の到達時刻の目安（正面往復）と検出ゲート
t_target = t_fire + 2 * target_distance * n_water / c0;
gate = abs(t_vec - t_target) <= 2 * pulse_fwhm;

res = cell(1, num_water);
y_all = zeros(num_water, num_samples);
warn_state = warning('query', 'underwater_mc_lidar_simulation:compression');
warning('off', 'underwater_mc_lidar_simulation:compression');
for k = 1:num_water
params_k = params;
params_k.water.a = water_a(k);
params_k.water.b = water_b(k);
params_k.water.g = water_g;
[y_k, info_k] = underwater_mc_lidar_simulation(x_pulse, params_k);
res{k} = info_k;
y_all(k, :) = y_k(1, :);
end
warning(warn_state.state, 'underwater_mc_lidar_simulation:compression');

% 4. 結果の表示
fprintf('  -> Rx #1 (ターゲット距離 %.1f m, 送信ピーク %.1f W, 光子数 %.2g / 水質)\n', ...
target_distance, P_peak, params.mc.num_photons);
fprintf('     雑音: 光入力換算 %.3g nW rms, 検出閾値: 雑音 rms の %.0f 倍, ゲート: %.1f ns ± %.0f ns\n\n', ...
params.noise_params.input_noise_rms * 1e9, det_threshold, t_target * 1e9, 2 * pulse_fwhm * 1e9);
fprintf('     %-14s  c [1/m]  P_rx(ゲート内) [W]  弾道成分 [W]  散乱による増加  後方散乱ピーク [W]  SNR [dB]  検出  計算時間 [s]  イベント/光子  相対誤差(全体/ターゲット)\n', '水質');
for k = 1:num_water
info_k = res{k};
P_gate = max(info_k.P_rx(1, gate));
P_bal = max(info_k.P_rx_ballistic(1, gate));
P_bs = max(info_k.P_rx_water(1, :));
v_gate = max(y_all(k, gate));
snr_gate = 20 * log10(max(info_k.v_signal(1, gate)) / info_k.noise_rms_v);
detected = v_gate >= det_threshold * info_k.noise_rms_v;
det_str = '不可';
if detected
    det_str = '可';
end
gain_str = '---';
if P_bal > 0
    gain_str = sprintf('x %.2f', P_gate / P_bal);
end
fprintf('     %-14s  %5.3f    %.3e           %.3e     %-10s     %.3e          %6.1f    %-4s  %8.1f      %8.1f      %5.1f%% / %5.1f%%\n', ...
    water_label{k}, water_a(k) + water_b(k), P_gate, P_bal, gain_str, P_bs, snr_gate, det_str, ...
    info_k.mc.runtime_s, info_k.mc.num_events / info_k.mc.num_photons, 100 * info_k.mc.relative_std_error(1), 100 * info_k.mc.relative_std_error_target(1));
end
fprintf('\n     (散乱による増加 = ゲート内の受光光パワー / 散乱なしの弾道成分。相対誤差は受光エネルギーのモンテカルロ標準誤差)\n\n');

% 5. 可視化: 水質ごとの受光光パワー（成分別、対数軸）と受光モジュール出力
figure('Name', 'Underwater MC LiDAR - Received Optical Power by Water Type', ...
'NumberTitle', 'off', 'Position', [60, 60, 1000, 820]);
colors = lines(4);
comp_names = {'Total', 'Target Return', 'Water Backscatter', 'Ballistic Target (Beer-Lambert)'};
comp_styles = {'-', '-', '-', '--'};
ax_k = zeros(1, num_water);
for k = 1:num_water
info_k = res{k};
subplot(num_water, 1, k);
comps = {info_k.P_rx(1, :), info_k.P_rx_target(1, :), info_k.P_rx_water(1, :), info_k.P_rx_ballistic(1, :)};
p_pos = [comps{:}];
p_pos = p_pos(p_pos > 0);
for j = 1:4
    pj = comps{j};
    pj(pj <= 0) = NaN;
    semilogy(t_vec * 1e9, pj, comp_styles{j}, 'Color', colors(j, :), 'LineWidth', 1.3, ...
        'DisplayName', comp_names{j});
    hold on;
end
semilogy([0, t_stop * 1e9], [1, 1] * params.noise_params.input_noise_rms, 'k:', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Noise Level (%.3g nW rms)', params.noise_params.input_noise_rms * 1e9));
grid on;
if ~isempty(p_pos)
    ylim([max(p_pos) * 1e-8, max(p_pos) * 3]);
end
xlim([0, t_stop * 1e9]);
ylabel('Rx #1 P_{rx} [W]');
set(gca, 'YMinorGrid', 'off');
title(sprintf('%s (a = %.3g, b = %.3g /m, g = %.2g)', ...
    water_label{k}, water_a(k), water_b(k), water_g));
ax_k(k) = gca;
if k == 1
    legend('Location', 'northeastoutside');
end
if k == num_water
    xlabel('Time [ns]');
end
end

% 凡例を軸の外に置いた上で、各段の時間軸の幅をそろえる
drawnow;
pos_k = cell2mat(get(ax_k, 'Position'));
for k = 1:num_water
set(ax_k(k), 'Position', [pos_k(k, 1), pos_k(k, 2), min(pos_k(:, 3)), pos_k(k, 4)]);
end

figure('Name', 'Underwater MC LiDAR - Photoreceiver Output by Water Type', ...
'NumberTitle', 'off', 'Position', [900, 60, 1000, 820]);
ax_k = zeros(1, num_water);
for k = 1:num_water
info_k = res{k};
subplot(num_water, 1, k);
plot(t_vec * 1e9, y_all(k, :) * 1e3, 'Color', [0.2, 0.5, 0.8], 'LineWidth', 1.0, 'DisplayName', 'Output with Noise');
hold on;
plot(t_vec * 1e9, info_k.v_signal(1, :) * 1e3, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Noise-Free Signal');
plot([0, t_stop * 1e9], [1, 1] * det_threshold * info_k.noise_rms_v * 1e3, 'r:', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Detection Threshold (%.0f\\sigma)', det_threshold));
grid on;
xlim([0, t_stop * 1e9]);
ylabel('v_{out} [mV]');
title(sprintf('%s: Rx #1 Output (Target Expected at %.1f ns)', water_label{k}, t_target * 1e9));
ax_k(k) = gca;
if k == 1
    legend('Location', 'northeastoutside');
end
if k == num_water
    xlabel('Time [ns]');
end
end

% 凡例を軸の外に置いた上で、各段の時間軸の幅をそろえる
drawnow;
pos_k = cell2mat(get(ax_k, 'Position'));
for k = 1:num_water
set(ax_k(k), 'Position', [pos_k(k, 1), pos_k(k, 2), min(pos_k(:, 3)), pos_k(k, 4)]);
end

fprintf('すべてのデモケースが正常に終了しました。\n');
fprintf('====================================================\n');
