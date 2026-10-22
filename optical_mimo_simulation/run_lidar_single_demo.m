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

% --- 受光器共通仕様 ---
params_ts.rx_common.aperture_width = 0.01;      % 受光開口幅 10 mm
params_ts.rx_common.fov_deg        = 45.0;      % 視野角 半角 45 deg
params_ts.rx_common.responsivity   = 0.60;      % PD感度 0.6 A/W
params_ts.rx_common.tia_gain       = 1.0e4;     % TIAゲイン 10 kOhm
params_ts.rx_common.hpf_cutoff     = 0.05e6;    % 低域遮断 0.05 MHz (C5658, AC結合)
params_ts.rx_common.lpf_cutoff     = 1.0e9;     % 高域遮断 1 GHz (C5658)

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

% 実測値: 受光波形のピーク時刻（3点放物線補間でサブサンプル推定）
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

fprintf('     ・Rx #%d: ピーク電圧 %.3e V\n', n, max(y_ts(n, :)));
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
params_static.enable_plot = false;

[y_static, H_elec_dc, ~, ~] = optical_mimo_simulation(x_static, params_static);

fprintf('  -> 送信光パワー x_static [W]:\n');
disp(x_static');
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

y_media = zeros(num_media, num_samples_w);
v_peak = zeros(1, num_media);

fprintf('  -> Rx #%d の受光パルス (ターゲット距離 %.1f m, 発射時刻 %.2f ns):\n', ...
rx_view, target_distance, t_fire * 1e9);
fprintf('      n     c [1/m]  ピーク電圧 [V]  空気比 [dB]  TOF [ns]  推定距離 [m]  媒体\n');
for k = 1:num_media
params_w = params_ts;
params_w.enable_plot = false;
params_w.refractive_index = media_n(k);
params_w.medium_extinction = media_ext(k);

y_w = optical_mimo_simulation(x_pulse_w, params_w);
y_media(k, :) = y_w(rx_view, :);

% 受光波形のピーク時刻（3点放物線補間でサブサンプル推定）
[v_peak(k), idx] = max(y_media(k, :));
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

% 受光電圧は光パワーに比例するため、光パワー比は 10 log10 で表す
fprintf('     %.2f  %5.2f    %.3e      %8.2f    %7.3f   %7.3f       %s\n', ...
    media_n(k), media_ext(k), v_peak(k), 10 * log10(v_peak(k) / v_peak(1)), ...
    tof_meas * 1e9, c_medium_k * (tof_meas - tau_rx) / 2, media_name{k});


end
fprintf('     (推定距離は媒体中の光速 c0/n と受光回路群遅延 %.3f ns の補正を用いた片道換算値)\n\n', tau_rx * 1e9);

% 媒体ごとの受光パルス比較の可視化
figure('Name', 'Pulse LiDAR - Underwater Propagation', ...
'NumberTitle', 'off', 'Position', [100, 100, 820, 560]);
colors_media = lines(num_media);

% 上段: 各媒体のピークで正規化した受光パルス波形（TOF の違い）
ax_w1 = subplot(2, 1, 1);
hold on; grid on;
for k = 1:num_media
plot(t_vec_w * 1e9, y_media(k, :) / v_peak(k), 'LineWidth', 1.5, 'Color', colors_media(k, :), ...
    'DisplayName', sprintf('%s (n = %.2f, c = %.2f /m)', media_label{k}, media_n(k), media_ext(k)));
end
xlabel('Time [ns]');
ylabel('Normalized Output');
xlim([0, t_stop_w * 1e9]);
title(sprintf('Rx #%d Photoreceiver Output v_{out}(t) (Normalized to Each Peak)', rx_view));
legend('Location', 'northeastoutside');

% 下段: 空気中を基準とした受光パルスのピーク比 [dB]（消散による減衰、光パワー比 = 10 log10）
ax_w2 = subplot(2, 1, 2);
peak_db = 10 * log10(v_peak / v_peak(1));
bar(1:num_media, peak_db, 'FaceColor', [0.2, 0.5, 0.8]);
grid on;
y_floor = min(min(peak_db) * 1.25, -10);
ylim([y_floor, -0.1 * y_floor]);
for k = 1:num_media
db_label = round(peak_db(k) * 10) / 10;
if abs(db_label) < 0.05
    db_label = 0;               % "-0.0 dB" 表示の回避
end
text(k, min(peak_db(k), 0) - 0.04 * abs(y_floor), sprintf('%.1f dB', db_label), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);
end
set(gca, 'XTick', 1:num_media);
set(gca, 'XTickLabel', media_tick);
xlabel('Medium (see legend above for n and c)');
ylabel('Peak Relative to Air [dB]');
title(sprintf('Rx #%d Peak Output Relative to Air (Target Distance %.1f m)', rx_view, target_distance));

% 上下の軸の幅をそろえる（上段の凡例は軸の外）
drawnow;
pos_w1 = get(ax_w1, 'Position');
pos_w2 = get(ax_w2, 'Position');
set(ax_w2, 'Position', [pos_w2(1), pos_w2(2), pos_w1(3), pos_w2(4)]);

fprintf('すべてのデモケースが正常に終了しました。\n');
fprintf('====================================================\n');
