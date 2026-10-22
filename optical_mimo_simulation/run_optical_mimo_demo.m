% =========================================================================
% スクリプト名: run_optical_mimo_demo.m
%
% 【概要】
% 光路長差による伝搬時間遅延 (Time of Flight: TOF) を考慮した
% 2次元光MIMO空間伝搬シミュレーション関数の高周波検証デモスクリプトです。
%
% 【検証ポイント】
%   - fs = 2 GHz, f_carrier = 200 MHz による高周波正弦波変調
%   - 連続送信の定常状態 (initial_state = 'steady') における受光モジュール
%     (C5658: 0.05 MHz ～ 1 GHz, AC結合) 出力波形
%     ※ 'zero' に変更すると、光の到達遅延（約 13.6 ns の無信号時間）と点灯直後の過渡を観測可能
%   - 光路長差（Rx #1 と Rx #2 の幾何距離差）によって生じる受光波形の目視可能な位相差
% =========================================================================

clear;
close all;
clc;

fprintf('\n');
fprintf('  =======================================================\n');
fprintf('    光MIMO 2次元空間伝搬シミュレーション (TOF遅延モデル対応)\n');
fprintf('  =======================================================\n\n');

%% ========================================================================
% 【ケース 1】高周波変調時系列シミュレーション（fs = 1 GHz, f_carrier = 50 MHz）
% ========================================================================
fprintf('--- [ケース 1] 高周波時系列シミュレーション (TOF位相差検証) ---\n');

% 1. 時間軸および送信光強度信号 x(t) の定義
fs = 2.0e9;                     % サンプリング周波数 2 GHz (0.5 ns 分解能)
dt = 1 / fs;
t_stop = 50e-9;                % 解析時間 50 ns (200 MHz の 10 周期分: 定常状態では信号周期の整数倍とする)
num_samples = round(t_stop * fs);
t_vec = (0:(num_samples - 1)) * dt;   % 周期の継ぎ目が重複しないよう終端 t_stop は含めない

% 送信搬送周波数 50 MHz（周期 20 ns）
% ※ 100 MHz や 200 MHz に引き上げると、受光波形の位相差がさらに拡大します
f_carrier = 200.0e6;

phase_tx1_rad = 0.0;            % Tx #1: 位相 0 deg
phase_tx2_rad = 10.0 * pi / 180;% Tx #2: 位相 10 deg ずらし

P_bias = 0.50;                  % 直流光バイアス [W]
A_mod  = 0.35;                  % 変調振幅 [W]

x_mimo = zeros(2, num_samples);
x_mimo(1, :) = P_bias + A_mod * sin(2 * pi * f_carrier * t_vec + phase_tx1_rad);
x_mimo(2, :) = P_bias + A_mod *0.1 * sin(2 * pi * f_carrier * t_vec + phase_tx2_rad);

% 2. シミュレーションパラメータの設定
params_ts = struct();
params_ts.fs = fs;
params_ts.enable_plot = true;            % 可視化を有効化
params_ts.medium_extinction = 0.0;       % 空気中
params_ts.refractive_index = 1.0;        % 空気中の屈折率 (水中なら 1.33)
params_ts.initial_state = 'steady';      % 連続送信の定常状態 ('zero': t < 0 は無送信)

% --- 送信アレイ (Tx) の配置・角度 (y軸上 x=0 に配置, ターゲット面に対し垂直 0 deg) ---
params_ts.tx_array(1).pos          = [0.0,  0.50];
params_ts.tx_array(1).type         = 'LD';
params_ts.tx_array(1).pointing_deg =   0.0;
params_ts.tx_array(1).fwhm         = 25.0;

params_ts.tx_array(2).pos          = [0.0,  0.25];
params_ts.tx_array(2).type         = 'LD';
params_ts.tx_array(2).pointing_deg =   0.0;
params_ts.tx_array(2).fwhm         = 25.0;

% --- 受光アレイ (Rx) の配置・受光軸 ---
params_ts.rx_array(1).pos          = [0.0, -0.25];
params_ts.rx_array(1).pointing_deg =   0.0;      % 内向き 5 deg

params_ts.rx_array(2).pos          = [0.0, -1.050];
params_ts.rx_array(2).pointing_deg =   0.0;      % 内向き 10 deg

% --- 受光器共通仕様 ---
params_ts.rx_common.aperture_width = 0.01;      % 受光開口幅 10 mm
params_ts.rx_common.fov_deg        = 45.0;      % 視野角 半角 45 deg
params_ts.rx_common.responsivity   = 0.60;      % PD感度 0.6 A/W
params_ts.rx_common.tia_gain       = 1.0e4;     % TIAゲイン 10 kOhm
params_ts.rx_common.hpf_cutoff     = 0.05e6;    % 低域遮断 0.05 MHz (C5658, AC結合)
params_ts.rx_common.lpf_cutoff     = 1.0e9;     % 高域遮断 1 GHz (C5658)

% --- 反射ターゲット仕様 (距離 2.0m, 長さ 1.8m) ---
params_ts.target.distance    = 2.0;             % 正面距離 2.0 m
params_ts.target.length      = 1.8;             % 長さ 1.8 m
params_ts.target.num_patches = 90;              % 90パッチ分割
params_ts.target.reflectance = 0.8;             % 反射率 0.8

% 3. 関数の呼び出し
[y_ts, H_elec_ts, H_opt_ts, info_ts] = optical_mimo_simulation(x_mimo, params_ts);

% 4. チャネル解析情報のコンソール出力
fprintf('  -> 計算完了:\n');
fprintf('     ・静的電気チャネル行列 H_elec (DC極限) [V/W]:\n');
disp(H_elec_ts);
fprintf('     ・チャネル条件数 (Condition Number) : %.3f\n', info_ts.condition_number);
fprintf('     ・特異値 (Singular Values)         : [%.2e, %.2e]\n', ...
info_ts.singular_values(1), info_ts.singular_values(2));
fprintf('     ・行列のランク (Rank)              : %d\n\n', info_ts.rank);

%% ========================================================================
% 【ケース 2】静的DC光パワー解析（スナップショット）
% ========================================================================
fprintf('--- [ケース 2] 静的DC光パワー解析 (DCブロック無効) ---\n');

x_static = [1.0; 0.8];
params_static = params_ts;
params_static.rx_common.hpf_cutoff = 0;        % DC解析のため HPF (AC結合) を無効化
params_static.enable_plot = false;

[y_static, H_elec_dc, ~, ~] = optical_mimo_simulation(x_static, params_static);

fprintf('  -> 送信光パワー x_static [W]:\n');
disp(x_static');
fprintf('  -> 受光直流電圧 y_static [V]:\n');
disp(y_static');

%% ========================================================================
% 【ケース 3】水中屈折率による伝搬遅延の増大検証
% ========================================================================
fprintf('--- [ケース 3] 媒体屈折率による伝搬速度・遅延時間の比較 ---\n');

c0 = 2.99792458e8;
d_roundtrip = 2.0 * 2.0; % 往復約 4m
tau_air = d_roundtrip / c0;
tau_water = d_roundtrip / (c0 / 1.333);

fprintf('  空気中 (n = 1.00): 往復 4m の光到達時間 = %.2f ns\n', tau_air * 1e9);
fprintf('  純水中 (n = 1.33): 往復 4m の光到達時間 = %.2f ns (+%.2f ns 遅延増)\n\n', ...
tau_water * 1e9, (tau_water - tau_air) * 1e9);

fprintf('すべてのデモケースが正常に終了しました。\n');
fprintf('====================================================\n');

