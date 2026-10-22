% =========================================================================
% 関数名: optical_mimo_simulation
%
% 【関数呼び出し形式】
%   [y, H_elec, H_opt, info] = optical_mimo_simulation(x, params);
%
% 【概要】
% 純粋な2次元直交座標系 (x, y) における光MIMO（インコヒーレント光・LiDAR反射型）
% 空間伝搬シミュレーション関数。
% 光路長差による伝搬時間遅延（Time of Flight: TOF）およびマルチパス広がりを
% 周波数領域複素伝達関数 H(f) により厳密にモデル化。
%
% 【2次元座標系の定義】
%   - x 軸: 前方進行方向（右向き正面が 0 度、距離 [m]）
%   - y 軸: アレイ配置・横方向軸 [m]
%   - 偏角 theta: 正面 (+x 方向) を 0 度とした xy 平面内の回転角 [deg]
%   - 光源 (Tx) / 受光器 (Rx): y 軸上 (x = 0) に配置
%   - 反射ターゲット: x = D に位置する y 軸方向の線分（微小線素パッチ分割）
%   - 空間光強度: 各座標 (x, y) における合成光パワー密度 [W/m]
%
% 【受光回路モデル】
%   - 受光モジュール (APD + 内蔵アンプ) の帯域を 1次HPF (低域遮断 hpf_cutoff) と
%     1次LPF (高域遮断 lpf_cutoff) で模擬 (例: 浜松ホトニクス C5658 = 0.05 MHz ～ 1 GHz)
%   - params.initial_state により送受信の初期状態を選択
%       'zero'  : t < 0 は無送信・受光回路は静止状態（光の到達遅延・点灯直後の過渡を観測）
%       'steady': 入力 x を1周期として繰り返し送信している定常状態
%                 （x の長さは信号周期の整数倍とすること）
% =========================================================================

function [y, H_elec, H_opt, info] = optical_mimo_simulation(x, params)

% ---------------------------------------------------------------------
% 1. 入力引数の検証およびパラメータの補完
% ---------------------------------------------------------------------
default_params = get_default_params();
if nargin < 2 || isempty(params)
params = default_params;
else
params = fill_missing_fields(params, default_params);
end

% 旧パラメータ dc_block との互換: dc_block = false が指定された場合は HPF を無効化
if isfield(params.rx_common, 'dc_block') && ~params.rx_common.dc_block
params.rx_common.hpf_cutoff = 0;
end

% ---------------------------------------------------------------------
% 2. 入力光信号 x の次元チェックと時間軸 t_vec の生成
% ---------------------------------------------------------------------
M_expected = length(params.tx_array);

if nargin < 1 || isempty(x)
% デフォルトテスト波形: 50 MHz 高周波正弦波、Tx #2 は位相 10 度ずらし
fs = params.fs;
t_vec = 0:(1/fs):100e-9;
num_samples = length(t_vec);
x = zeros(M_expected, num_samples);
f_carrier = 50.0e6;
for m = 1:M_expected
phase_offset = (m - 1) * (10.0 * pi / 180);
x(m, :) = 0.5 + 0.35 * sin(2 * pi * f_carrier * t_vec + phase_offset);
end
else
% 入力信号が転置されている場合の補正 (M x 1 を期待)
if size(x, 1) == 1 && size(x, 2) == M_expected
x = x.';
end
num_samples = size(x, 2);
if num_samples > 1
t_vec = (0:(num_samples - 1)) / params.fs;
else
t_vec = 0;  % 静的スナップショット
end
end

% ---------------------------------------------------------------------
% 3. 2次元反射ターゲット幾何（線分パッチ群）の生成
% ---------------------------------------------------------------------
target_geo = generate_target_geometry_2d(params.target);

% ---------------------------------------------------------------------
% 4. 2次元光チャネル行列 H_opt およびパッチ別遅延データの計算
% ---------------------------------------------------------------------
[H_opt, patch_data] = calculate_optical_channel_2d(params.tx_array, params.rx_array, ...
target_geo, params.rx_common, ...
params.medium_extinction);

% ---------------------------------------------------------------------
% 5. 静的電気チャネル行列 H_elec の計算 (DC / 準静的極限ゲイン)
% ---------------------------------------------------------------------
G_tia = params.rx_common.tia_gain;
R_resp = params.rx_common.responsivity;
H_elec = G_tia * R_resp * H_opt;

% ---------------------------------------------------------------------
% 6. 受光出力信号 y の算出（TOF伝搬遅延モデル対応）
% ---------------------------------------------------------------------
c_speed = 2.99792458e8 / params.refractive_index;
is_steady = strcmpi(params.initial_state, 'steady');

if num_samples > 1
% 時系列信号: パッチごとの伝搬遅延と受光回路の周波数特性を周波数領域伝達関数 H(f) で畳み込み
%   steady: 巡回畳み込み（周期送信の定常状態）、HPF も定常応答として周波数領域で適用
%   zero  : 線形畳み込み、HPF は時定数が解析区間より長いため時間領域の零初期状態 IIR で適用
dH_elec_patches = G_tia * R_resp * patch_data.dH;
tau_patches = patch_data.path_length / c_speed;
y_rx = simulate_tof_propagation(x, dH_elec_patches, tau_patches, params.fs, ...
params.rx_common, is_steady);
if ~is_steady
y_rx = apply_rx_highpass_zero_state(y_rx, params.rx_common.hpf_cutoff, params.fs);
end
else
% 静的スナップショット: 遅延は位相に影響しないため実数行列の瞬時積
% (DC解析のため、HPF 有効時は定常状態で DC 成分が遮断され出力 0)
y_rx = H_elec * x;
if params.rx_common.hpf_cutoff > 0
y_rx = zeros(size(y_rx));
end
end

y = apply_channel_noise(y_rx, params.noise_params);

% ---------------------------------------------------------------------
% 7. チャネル解析情報の収集
% ---------------------------------------------------------------------
sing_vals = svd(H_elec);
cond_num = cond(H_elec);

info.singular_values = sing_vals;
info.condition_number = cond_num;
info.rank = rank(H_elec);
info.target_geometry = target_geo;
info.H_optical = H_opt;
info.H_electrical = H_elec;
info.time_vector = t_vec;
info.patch_data = patch_data;

% ---------------------------------------------------------------------
% 8. 可視化レンダリング（enable_plot が true の場合のみ実行）
% ---------------------------------------------------------------------
if isfield(params, 'enable_plot') && params.enable_plot
render_simulation_visuals_2d(params, target_geo, H_elec, x, y, t_vec, cond_num);
end

end

% =========================================================================
% 内部関数群: 2次元デフォルトパラメータの取得 (2Tx x 2Rx)
% =========================================================================
function p = get_default_params()
p.medium_extinction = 0.0;     % 消散係数 (0.0: 空気中, >0: 水中 [1/m])
p.refractive_index = 1.0;      % 媒体屈折率 (空気中: 1.0, 純水: 1.33)
p.fs = 1.0e9;                  % サンプリング周波数 [Hz] (1 GHz)
p.initial_state = 'zero';      % 送受信の初期状態 ('zero': t<0 は無送信, 'steady': x の周期送信の定常状態)
p.enable_plot = true;          % 可視化フラグ
p.intensity_plot_scale = 'dB'; % 光強度分布図の表示スケール ('dB' または 'linear')
p.intensity_db_range = 40.0;   % dB表示時のダイナミックレンジ [dB] (最大値からの表示下限)

% 送信アレイ (Tx: x=0 の y軸上に配置, ターゲット平板に対し垂直 0 deg)
% Tx #1: (0, 0.50) [m]
p.tx_array(1).pos = [0.0,  0.50];
p.tx_array(1).type = 'LD';
p.tx_array(1).pointing_deg = 0.0;
p.tx_array(1).fwhm = 25.0;

% Tx #2: (0, 0.25) [m]
p.tx_array(2).pos = [0.0,  0.25];
p.tx_array(2).type = 'LD';
p.tx_array(2).pointing_deg = 0.0;
p.tx_array(2).fwhm = 25.0;

% 受光アレイ (Rx: x=0 の y軸上に配置)
% Rx #1: (0, -0.25) [m], 内向き +5 deg
p.rx_array(1).pos = [0.0, -0.25];
p.rx_array(1).pointing_deg =   5.0;

% Rx #2: (0, -0.50) [m], 内向き +10 deg
p.rx_array(2).pos = [0.0, -0.50];
p.rx_array(2).pointing_deg =  10.0;

% 受光器共通仕様 (2D受光開口幅)
p.rx_common.aperture_width = 0.01;         % 2D受光開口幅 [m] (10mm)
p.rx_common.fov_deg = 45.0;                % 視野角 半角 [deg]
p.rx_common.responsivity = 0.60;           % PD感度 [A/W]
p.rx_common.tia_gain = 1.0e4;              % TIAゲイン [V/A]
p.rx_common.hpf_cutoff = 0.05e6;           % 受光回路の低域遮断周波数 [Hz] (1次HPF, 0: 無効) C5658: 0.05 MHz
p.rx_common.lpf_cutoff = 0;                % 受光回路の高域遮断周波数 [Hz] (1次LPF, 0: 無効) C5658: 1 GHz

% 反射ターゲット (x = distance に置かれた y軸方向の線分壁)
p.target.distance = 2.0;                   % x 軸距離 [m]
p.target.length = 1.8;                     % y 方向の長さ [m]
p.target.num_patches = 90;                 % 分割パッチ数 (線素数)
p.target.reflectance = 0.80;               % 拡散反射率 (rho = 0.8)

% ノイズ設定
p.noise_params.enable = false;
p.noise_params.snr_db = 30.0;
end

% =========================================================================
% 内部関数: 未定義パラメータの再帰的補完
% =========================================================================
function user_s = fill_missing_fields(user_s, def_s)
if ~isstruct(user_s)
user_s = def_s;
return;
end

fnames = fieldnames(def_s);
for i = 1:length(fnames)
f = fnames{i};
if ~isfield(user_s, f)
user_s.(f) = def_s.(f);
else
if isstruct(def_s.(f)) && isstruct(user_s.(f)) && ...
numel(def_s.(f)) == 1 && numel(user_s.(f)) == 1
user_s.(f) = fill_missing_fields(user_s.(f), def_s.(f));
end
end
end
end

% =========================================================================
% 内部関数: 2次元反射ターゲット線分パッチ生成
% =========================================================================
function target = generate_target_geometry_2d(t_params)
K = t_params.num_patches;
L = t_params.length;
D = t_params.distance;

y_edges = linspace(-L / 2, L / 2, K + 1);
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
dl = L / K;

target.p = zeros(K, 2);
target.p(:, 1) = D;                        % x 座標 (距離 D)
target.p(:, 2) = y_centers(:);             % y 座標
target.n = repmat([-1.0, 0.0], K, 1);      % 法線は正面を向く (-x 方向)
target.dl = dl * ones(K, 1);               % 各パッチの線素長 [m]
target.rho = t_params.reflectance;
end

% =========================================================================
% 内部関数: 2次元光伝搬結合係数およびパッチ別光路長の計算
% =========================================================================
function [H_opt, patch_data] = calculate_optical_channel_2d(tx_array, rx_array, target, rx_common, c_ext)
M = length(tx_array);
N = length(rx_array);
K = size(target.p, 1);

H_opt = zeros(N, M);
dH_tensor = zeros(N, M, K);
path_tensor = zeros(N, M, K);

for m = 1:M
tx_pos = tx_array(m).pos(1:2);
pointing_tx = tx_array(m).pointing_deg;
v_tx = [cosd(pointing_tx), sind(pointing_tx)];
fwhm_m = tx_array(m).fwhm;
type_m = tx_array(m).type;

% Tx_m から全パッチへの幾何計算
d1_vec = target.p - repmat(tx_pos, K, 1);
r1 = sqrt(sum(d1_vec.^2, 2));
u1 = d1_vec ./ repmat(max(r1, 1e-6), 1, 2);

cos_alpha = sum(repmat(v_tx, K, 1) .* u1, 2);
cos_alpha = max(min(cos_alpha, 1.0), -1.0);
alpha_deg = acosd(cos_alpha);

cos_psi1 = sum(-u1 .* target.n, 2);
cos_psi1 = max(cos_psi1, 0.0);

R_tx = eval_tx_radiation_2d(alpha_deg, fwhm_m, type_m);
E_patch = (R_tx .* cos_psi1 ./ max(r1, 1e-6)) .* exp(-c_ext * r1);

for n = 1:N
    rx_pos = rx_array(n).pos(1:2);
    pointing_rx = rx_array(n).pointing_deg;
    v_rx = [cosd(pointing_rx), sind(pointing_rx)];

    % パッチから Rx_n への幾何計算
    d2_vec = repmat(rx_pos, K, 1) - target.p;
    r2 = sqrt(sum(d2_vec.^2, 2));
    u2 = d2_vec ./ repmat(max(r2, 1e-6), 1, 2);

    cos_theta2 = sum(u2 .* target.n, 2);
    cos_theta2 = max(cos_theta2, 0.0);

    cos_psi2 = sum(-u2 .* repmat(v_rx, K, 1), 2);
    cos_psi2 = max(min(cos_psi2, 1.0), -1.0);
    psi2_deg = acosd(cos_psi2);

    fov_mask = double(psi2_deg <= rx_common.fov_deg);

    dH = (E_patch * (target.rho / 2.0)) .* ...
         (cos_theta2 .* cos_psi2 ./ max(r2, 1e-6)) .* ...
         rx_common.aperture_width .* exp(-c_ext * r2) .* fov_mask .* target.dl;

    H_opt(n, m) = sum(dH);
    dH_tensor(n, m, :) = dH;
    path_tensor(n, m, :) = r1 + r2;   % 往復光路長 [m]
end


end

patch_data.dH = dH_tensor;
patch_data.path_length = path_tensor;
end

% =========================================================================
% 内部関数: TOF周波数領域伝達関数による時系列信号の厳密伝搬シミュレーション
%   is_periodic = true : x を1周期とする周期送信の定常状態（巡回畳み込み、HPF も周波数領域で適用）
%   is_periodic = false: t < 0 は無送信（ゼロパディングによる線形畳み込み、HPF は呼び出し側で適用）
% =========================================================================
function y_out = simulate_tof_propagation(x, dH_elec_patches, tau_patches, fs, rx_common, is_periodic)
[M, num_samples] = size(x);
N = size(dH_elec_patches, 1);

if is_periodic
n_pad = 0;
else
max_tau = max(tau_patches(:));
n_pad = ceil(max_tau * fs) + 64;
end
n_fft = num_samples + n_pad;

% 正の周波数ベクトル（DCからナイキスト周波数まで）
num_pos = floor(n_fft / 2) + 1;
f_pos = (0:(num_pos - 1)) * (fs / n_fft);

% 受光回路の周波数特性（LPF、定常状態では HPF も含む）
H_rx_half = eval_rx_frequency_response(f_pos, rx_common, is_periodic);

y_out = zeros(N, num_samples);

% 送信波形をゼロパディングしてFFT（定常状態ではパディングなしの巡回畳み込み）
X_fft = zeros(M, n_fft);
for m = 1:M
x_padded = [x(m, :), zeros(1, n_pad)];
X_fft(m, :) = fft(x_padded);
end

for n = 1:N
Y_fft_n = zeros(1, n_fft);
for m = 1:M
dH_km = squeeze(dH_elec_patches(n, m, :));     % [K x 1]
tau_km = squeeze(tau_patches(n, m, :));         % [K x 1]

    % 各パッチの位相回転 exp(-j 2*pi*f*tau) を積分した伝達関数 H(f)
    phase_mat = -2.0 * pi * (tau_km * f_pos);       % [K x num_pos]
    H_half = sum(repmat(dH_km, 1, num_pos) .* exp(1i * phase_mat), 1) .* H_rx_half;

    % エルミート対称性を持たせて全周波数ベクトルを構築
    H_full = zeros(1, n_fft);
    H_full(1:num_pos) = H_half;
    if mod(n_fft, 2) == 0
        H_full((num_pos + 1):n_fft) = conj(H_half(end-1:-1:2));
    else
        H_full((num_pos + 1):n_fft) = conj(H_half(end:-1:2));
    end

    Y_fft_n = Y_fft_n + X_fft(m, :) .* H_full;
end
y_time = real(ifft(Y_fft_n));
y_out(n, :) = y_time(1:num_samples);


end
end

% =========================================================================
% 内部関数: 2次元放射パターン評価 (LED / LD)
%   戻り値 R(alpha) [1/rad] は -90 ～ +90 deg の角度積分が 1 となるよう正規化。
%   送信光パワー P [W] に対する放射強度は P * R(alpha) [W/rad] となる。
% =========================================================================
function R = eval_tx_radiation_2d(alpha_deg, fwhm_deg, type)
if strcmpi(type, 'LED')
cos_half = cosd(min(fwhm_deg / 2.0, 89.0));
m_order = -log(2.0) / log(cos_half);
% 正規化積分: int_{-pi/2}^{pi/2} cos^m(a) da = sqrt(pi) * Gamma((m+1)/2) / Gamma(m/2+1)
norm_int = sqrt(pi) * exp(gammaln((m_order + 1) / 2) - gammaln(m_order / 2 + 1));
rad = alpha_deg * (pi / 180.0);
R = zeros(size(alpha_deg));
valid = abs(alpha_deg) <= 90.0;
R(valid) = (cos(rad(valid)).^m_order) / norm_int;
else
sigma = (fwhm_deg / 2.0) / sqrt(2 * log(2));
% 正規化積分: int_{-pi/2}^{pi/2} exp(-a^2 / (2 sigma^2)) da
sigma_rad = sigma * (pi / 180.0);
norm_int = sqrt(2 * pi) * sigma_rad * erf((pi / 2) / (sqrt(2) * sigma_rad));
R = exp(-0.5 * (alpha_deg ./ sigma).^2) / norm_int;
R(abs(alpha_deg) > 90.0) = 0.0;
end
end

% =========================================================================
% 内部関数: 受光回路の周波数特性 (1次LPF・1次HPF のアナログ伝達関数)
%   H_LPF(f) = 1 / (1 + j f / f_H),  H_HPF(f) = (j f / f_L) / (1 + j f / f_L)
% =========================================================================
function H_rx = eval_rx_frequency_response(f, rx_common, include_hpf)
H_rx = ones(size(f));

f_H = rx_common.lpf_cutoff;
if f_H > 0 && isfinite(f_H)
H_rx = H_rx ./ (1 + 1i * f / f_H);
end

f_L = rx_common.hpf_cutoff;
if include_hpf && f_L > 0
jf = 1i * f / f_L;
H_rx = H_rx .* (jf ./ (1 + jf));
end
end

% =========================================================================
% 内部関数: 受光回路 1次HPF（AC結合）の零初期状態応答
%   t < 0 で回路が静止している状態から、因果的な IIR フィルタとして時間領域で適用。
%   HPF の時定数 1/(2 pi f_L) は解析区間より十分長いため、周波数領域ではなく
%   双一次変換（プリワーピング付き）による時間領域処理とする。
% =========================================================================
function y_out = apply_rx_highpass_zero_state(y_in, f_L, fs)
y_out = y_in;
if f_L <= 0
return;
end

K = tan(pi * f_L / fs);
b = [1, -1] / (1 + K);
a = [1, (K - 1) / (1 + K)];
for n = 1:size(y_in, 1)
y_out(n, :) = filter(b, a, y_in(n, :));
end
end

% =========================================================================
% 内部関数: チャネル雑音の付加
% =========================================================================
function y_out = apply_channel_noise(y_in, noise_params)
if ~isfield(noise_params, 'enable') || ~noise_params.enable
y_out = y_in;
return;
end

y_out = y_in;
if isfield(noise_params, 'snr_db')
snr_linear = 10^(noise_params.snr_db / 10);
sig_power = mean(y_in(:).^2);
if sig_power <= 0
sig_power = 1.0;
end
noise_var = sig_power / max(snr_linear, 1e-12);
noise = sqrt(noise_var) * randn(size(y_in));
y_out = y_out + noise;
end
end

% =========================================================================
% 内部関数: 2次元可視化レンダリング（1グラフ1Figure、波形は上下1Figure）
% =========================================================================
function render_simulation_visuals_2d(params, target, H_elec, X, Y, t, cond_num)
M = length(params.tx_array);
N = length(params.rx_array);

colors_tx = lines(max(M, 3));
colors_rx = lines(max(N, 5));

D_target = target.p(1, 1);
y_min_t  = min(target.p(:, 2));
y_max_t  = max(target.p(:, 2));

% 表示範囲: ターゲット・Tx・Rx の y 座標をすべて含み、距離 D に応じた余白を付加
y_all = [y_min_t, y_max_t, ...
arrayfun(@(s) s.pos(2), params.tx_array), ...
arrayfun(@(s) s.pos(2), params.rx_array)];
y_margin = 0.15 * max(max(y_all) - min(y_all), 0.05 * D_target);
y_lim_plot = [min(y_all) - y_margin, max(y_all) + y_margin];
x_margin = 0.05 * D_target;

% ---------------------------------------------------------------------
% Figure 1: 2D幾何配置図 (2D Geometry Layout)
%   距離 D と y 方向の広がりの比が大きくても潰れないよう axis equal は使わない
%   (矢印は実際の指向方向を示すが、x と y の縮尺は一致しない)
% ---------------------------------------------------------------------
figure('Name', 'Optical MIMO 2D - Geometry Layout', ...
'NumberTitle', 'off', 'Position', [50, 480, 580, 420]);
hold on; grid on;
xlabel('x-axis (Propagation Direction) [m]');
ylabel('y-axis (Array Axis) [m]');
title('2D Geometry Layout (Tx/Rx on y-axis, Target at x=D)');

plot([D_target, D_target], [y_min_t, y_max_t], 'k-', 'LineWidth', 4.0, 'DisplayName', 'Target Line');

len_arrow_tx = 0.30 * D_target;   % Tx 矢印幹の長さ [m] (ターゲット距離に比例)
len_arrow_rx = 0.25 * D_target;   % Rx 矢印幹の長さ [m] (ターゲット距離に比例)

for m = 1:M
pos = params.tx_array(m).pos;
theta = params.tx_array(m).pointing_deg;
plot(pos(1), pos(2), '^', 'MarkerSize', 9, ...
'MarkerFaceColor', colors_tx(m, :), 'Color', 'k', ...
'DisplayName', sprintf('Tx #%d (%s, %+.0f deg)', m, params.tx_array(m).type, theta));
end

for n = 1:N
pos = params.rx_array(n).pos;
theta = params.rx_array(n).pointing_deg;
plot(pos(1), pos(2), 'o', 'MarkerSize', 8, ...
'MarkerFaceColor', colors_rx(n, :), 'Color', 'k', ...
'DisplayName', sprintf('Rx #%d (%+.0f deg)', n, theta));
end
legend('Location', 'northeastoutside');
xlim([-x_margin, D_target + x_margin]);
ylim(y_lim_plot);

% 矢印は凡例・表示範囲の確定後に描画（矢じりの形状を画面上で一定に保つため）
for m = 1:M
theta = params.tx_array(m).pointing_deg;
draw_pointing_arrow_2d(params.tx_array(m).pos, [cosd(theta), sind(theta)], len_arrow_tx, ...
    colors_tx(m, :), '-', 1.8);
end

for n = 1:N
% Rx 矢印: 細かく密なドット点線 (':')
theta = params.rx_array(n).pointing_deg;
draw_pointing_arrow_2d(params.rx_array(n).pos, [cosd(theta), sind(theta)], len_arrow_rx, ...
    colors_rx(n, :), ':', 1.5);
end

% ---------------------------------------------------------------------
% Figure 2: 空間 (x, y) における合成光強度分布ヒートマップ
% ---------------------------------------------------------------------
figure('Name', 'Optical MIMO 2D - Optical Intensity Distribution', ...
'NumberTitle', 'off', 'Position', [650, 480, 600, 420]);
hold on; grid on;

num_samples = size(X, 2);
if num_samples > 1
P_weights = mean(X, 2);
weight_mode_str = 'Time-Averaged Power';
else
P_weights = X;
weight_mode_str = 'Static Power';
end

gx = linspace(0.01, D_target, 80);
gy = linspace(y_lim_plot(1), y_lim_plot(2), 80);
[GX, GY] = meshgrid(gx, gy);
field_intensity = zeros(size(GX));

for m = 1:M
t_pos = params.tx_array(m).pos;
v_m = [cosd(params.tx_array(m).pointing_deg), sind(params.tx_array(m).pointing_deg)];
dx = GX - t_pos(1);
dy = GY - t_pos(2);
dist = sqrt(dx.^2 + dy.^2);
cos_ang = (dx * v_m(1) + dy * v_m(2)) ./ max(dist, 1e-6);
ang = acosd(max(min(cos_ang, 1.0), -1.0));
rad_pat = eval_tx_radiation_2d(ang, params.tx_array(m).fwhm, params.tx_array(m).type);

field_intensity = field_intensity + ...
    P_weights(m) * (rad_pat ./ max(dist, 0.05)) .* exp(-params.medium_extinction * dist);


end

use_db = strcmpi(params.intensity_plot_scale, 'dB');
if use_db
% dB表示: 最大値から intensity_db_range [dB] 下までを表示範囲とする
I_max = max(field_intensity(:));
if I_max <= 0
    I_max = 1e-30;
end
db_top = 10 * log10(I_max);
db_floor = db_top - params.intensity_db_range;
field_plot = 10 * log10(max(field_intensity, 10^(db_floor / 10)));
cb_label_str = 'Optical Power Density [dB(W/m)]';
scale_str = 'dB';
else
field_plot = field_intensity;
cb_label_str = 'Optical Power Density [W/m]';
scale_str = 'Linear';
end

contourf(GX, GY, field_plot, 25, 'LineStyle', 'none');
h_cb_field = colorbar;
ylabel(h_cb_field, cb_label_str);
if use_db
caxis([db_floor, db_top]);
end
try
colormap(gca, 'copper');
catch
colormap('copper');
end
plot([D_target, D_target], [y_min_t, y_max_t], 'c-', 'LineWidth', 3.0);
title(sprintf('Combined Optical Intensity Distribution (%s, %s)', weight_mode_str, scale_str));
xlabel('x-axis (Propagation Direction) [m]');
ylabel('y-axis (Array Axis) [m]');
xlim([-x_margin, D_target + x_margin]);
ylim(y_lim_plot);

% ---------------------------------------------------------------------
% Figure 3: 電気チャネル行列 H_elec ヒートマップ (DC極限)
% ---------------------------------------------------------------------
figure('Name', 'Optical MIMO 2D - Electrical Channel Matrix', ...
'NumberTitle', 'off', 'Position', [50, 40, 580, 390]);
imagesc(H_elec);
h_cb = colorbar;
ylabel(h_cb, 'Channel Gain [V/W]');
title(sprintf('Electrical Channel Matrix H_{elec} [V/W] (Cond: %.2f)', cond_num));
xlabel('Transmitter Tx (m)');
ylabel('Receiver Rx (n)');
set(gca, 'XTick', 1:M, 'YTick', 1:N);
set(gca, 'XTickLabel', arrayfun(@(i) sprintf('Tx #%d', i), 1:M, 'UniformOutput', false));
set(gca, 'YTickLabel', arrayfun(@(i) sprintf('Rx #%d', i), 1:N, 'UniformOutput', false));

h_min = min(H_elec(:));
h_max = max(H_elec(:));
h_range = max(h_max - h_min, 1e-12);
for n = 1:N
for m = 1:M
val = H_elec(n, m);
norm_val = (val - h_min) / h_range;
if norm_val > 0.55
text_col = 'k';
else
text_col = 'w';
end
text(m, n, sprintf('%.2e', val), ...
'HorizontalAlignment', 'center', 'Color', text_col, 'FontWeight', 'bold', 'FontSize', 9);
end
end

% ---------------------------------------------------------------------
% Figure 4: 送信光パワー波形（上段）および受光モジュール出力電圧波形（下段）
% ---------------------------------------------------------------------
figure('Name', 'Optical MIMO 2D - Transmitted Power and Photoreceiver Output', ...
'NumberTitle', 'off', 'Position', [650, 40, 760, 560]);

% 受光回路の帯域・初期状態の表示文字列
if strcmpi(params.initial_state, 'steady')
state_str = 'Steady State';
else
state_str = 'Zero Initial State';
end
rx_band_str = sprintf('f_L = %s, f_H = %s', ...
format_freq_str(params.rx_common.hpf_cutoff), format_freq_str(params.rx_common.lpf_cutoff));

if length(t) > 1
t_max = max(t);
if t_max < 1.0e-6
t_scale = 1.0e9;
t_unit = 'ns';
elseif t_max < 1.0e-3
t_scale = 1.0e6;
t_unit = 'us';
else
t_scale = 1.0e3;
t_unit = 'ms';
end
t_disp = t * t_scale;

% 上段: 送信光パワー波形 P_tx(t)
ax_tx = subplot(2, 1, 1);
hold on; grid on;
for m = 1:M
    plot(t_disp, X(m, :), 'LineWidth', 1.5, 'Color', colors_tx(m, :), ...
        'DisplayName', sprintf('Tx #%d Optical Power', m));
end
xlabel(sprintf('Time [%s]', t_unit));
ylabel('Optical Power [W]');
xlim([0, max(t_disp)]);
title('Transmitted Optical Power P_{tx}(t)');
legend('Location', 'northeastoutside');

% 下段: 受光モジュール出力電圧波形 v_out(t)
ax_rx = subplot(2, 1, 2);
hold on; grid on;
for n = 1:N
    plot(t_disp, Y(n, :), 'LineWidth', 1.5, 'Color', colors_rx(n, :), ...
        'DisplayName', sprintf('Rx #%d Photoreceiver Output', n));
end
xlabel(sprintf('Time [%s]', t_unit));
ylabel('Output Voltage [V]');
xlim([0, max(t_disp)]);
title({'Photoreceiver Output Voltage v_{out}(t)', sprintf('(%s, %s)', rx_band_str, state_str)});
legend('Location', 'northeastoutside');

% 凡例を軸の外に置いた上で、上下の時間軸の幅をそろえる
drawnow;
pos_tx = get(ax_tx, 'Position');
pos_rx = get(ax_rx, 'Position');
ax_width = min(pos_tx(3), pos_rx(3));
set(ax_tx, 'Position', [pos_tx(1), pos_tx(2), ax_width, pos_tx(4)]);
set(ax_rx, 'Position', [pos_rx(1), pos_rx(2), ax_width, pos_rx(4)]);


else
subplot(2, 1, 1);
bar(1:M, X, 'FaceColor', [0.85, 0.45, 0.2]);
xlabel('Transmitter Tx (m)');
ylabel('Optical Power [W]');
set(gca, 'XTick', 1:M);
set(gca, 'XTickLabel', arrayfun(@(i) sprintf('Tx #%d', i), 1:M, 'UniformOutput', false));
title('Transmitted Static Optical Power P_{tx}');

subplot(2, 1, 2);
bar(1:N, Y, 'FaceColor', [0.2, 0.5, 0.8]);
xlabel('Receiver Rx (n)');
ylabel('Output Voltage [V]');
set(gca, 'XTick', 1:N);
set(gca, 'XTickLabel', arrayfun(@(i) sprintf('Rx #%d', i), 1:N, 'UniformOutput', false));
title({'Photoreceiver Static Output Voltage v_{out}', sprintf('(%s, DC Analysis)', rx_band_str)});


end

end

% =========================================================================
% 内部関数: 指向方向矢印の描画（x, y の縮尺が異なる軸でも矢じり形状を一定に保つ）
%   矢印の幹はデータ座標上の実際の指向方向に描き、矢じりは画面ピクセル座標で
%   一定の大きさ・開き角となるよう計算してデータ座標に戻す。
% =========================================================================
function draw_pointing_arrow_2d(p0, v_dir, len, color, line_style, line_width)
p0 = p0(1:2);
tip = p0 + len * v_dir;

% 現在の軸における 1 ピクセルあたりのデータ量 [m/px]
drawnow;
ax = gca;
old_units = get(ax, 'Units');
set(ax, 'Units', 'pixels');
ax_pos = get(ax, 'Position');
set(ax, 'Units', old_units);
x_lim = get(ax, 'XLim');
y_lim = get(ax, 'YLim');
s_px = [diff(x_lim) / ax_pos(3), diff(y_lim) / ax_pos(4)];

% 矢印をピクセル座標に変換し、矢じりの2辺を求める
d_px = (tip - p0) ./ s_px;
L_px = norm(d_px);
if L_px <= 0
return;
end
u_px = d_px / L_px;
head_px = min(12.0, 0.35 * L_px);      % 矢じりの長さ [px]
head_ang = 25.0;                        % 矢じりの開き半角 [deg]
rot_p = [cosd(head_ang), -sind(head_ang); sind(head_ang), cosd(head_ang)];
rot_m = rot_p.';
w1 = tip - head_px * (rot_p * u_px.').' .* s_px;
w2 = tip - head_px * (rot_m * u_px.').' .* s_px;

plot([p0(1), tip(1)], [p0(2), tip(2)], 'Color', color, 'LineStyle', line_style, ...
'LineWidth', line_width, 'HandleVisibility', 'off');
plot([w1(1), tip(1), w2(1)], [w1(2), tip(2), w2(2)], 'Color', color, 'LineStyle', '-', ...
'LineWidth', line_width, 'HandleVisibility', 'off');
end

% =========================================================================
% 内部関数: 周波数の表示用文字列 (0 以下・非有限は 'Off')
% =========================================================================
function s = format_freq_str(f)
if ~(f > 0) || ~isfinite(f)
s = 'Off';
elseif f >= 1e9
s = sprintf('%g GHz', f / 1e9);
elseif f >= 1e6
s = sprintf('%g MHz', f / 1e6);
elseif f >= 1e3
s = sprintf('%g kHz', f / 1e3);
else
s = sprintf('%g Hz', f);
end
end
