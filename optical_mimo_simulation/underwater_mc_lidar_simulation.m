% =========================================================================
% 関数名: underwater_mc_lidar_simulation
%
% 【関数呼び出し形式】
%   [y, info] = underwater_mc_lidar_simulation(x, params);
%
% 【概要】
% 3次元のモンテカルロ法（半解析的モンテカルロ法: 局所推定法）による
% 水中パルス LiDAR の光伝搬シミュレーション関数。
% 吸収・散乱（Henyey-Greenstein 位相関数）・ターゲットでのランベルト反射を
% 光子の追跡で扱い、受光器ごとのインパルス応答 h(t) を求める。
% 送信光パワー波形 x(t) との畳み込みで受光光パワー P_rx(t) を求め、
% 受光回路（R * G、1次HPF/LPF、光入力換算の固定雑音）を通して出力電圧 y(t) を得る。
% optical_mimo_simulation.m（2次元・1回反射の決定論的モデル）とは独立した別モデル。
%
% 【3次元座標系】
%   - x 軸: 前方（ターゲット方向）、y 軸: 横方向、z 軸: 上下方向 [m]
%   - Tx / Rx: 面 x = 0 上に配置（x < 0 側は機器の前面とみなし、到達した光子は終了）
%   - ターゲット: x = D の平面上の矩形（幅 width: y 方向、高さ height: z 方向）、ランベルト反射
%   - 指向角 pointing_deg: xy 平面内の方位角（正面 +x が 0 deg）
%
% 【モンテカルロ法】
%   - 自由行程: 減衰係数 c = a + b の指数分布からサンプリング
%   - 散乱: 重みに単一散乱アルベド b / c を掛け、新しい方向は HG 位相関数からサンプリング
%   - 局所推定: 散乱点・ターゲット反射点ごとに、そこから受光開口へ直接届く割合
%     （位相関数またはランベルト反射 × 開口の立体角 × exp(-c d)）を到達時刻のビンに加算
%   - 低重みの光子はロシアンルーレットで確率的に終了
%   - 寄与の分類: target    = ターゲットで反射した後の光子による寄与
%                 water     = ターゲット到達前の水の散乱による寄与（後方散乱クラッタ）
%                 ballistic = 散乱を経ずにターゲットで1回反射した成分（target の一部）
%   - 送信は単発（t < 0 は無送信）、受光回路は零初期状態
% =========================================================================

function [y, info] = underwater_mc_lidar_simulation(x, params)

% ---------------------------------------------------------------------
% 1. パラメータの補完と入力の整形
% ---------------------------------------------------------------------
default_params = get_default_params();
if nargin < 2 || isempty(params)
params = default_params;
else
params = fill_missing_fields(params, default_params);
end

fs = params.fs;
x = x(:).';                                  % 送信光パワー波形 [W] (1 x サンプル数)
num_samples = length(x);
t_vec = (0:(num_samples - 1)) / fs;
N = length(params.rx_array);

% ---------------------------------------------------------------------
% 2. モンテカルロ法による受光器ごとのインパルス応答 h(t) の計算
% ---------------------------------------------------------------------
t_start = tic;
[h, mc_stats] = run_photon_transport(params, num_samples);
mc_stats.runtime_s = toc(t_start);

% ---------------------------------------------------------------------
% 3. 受光光パワー P_rx(t) = x(t) * h(t)（成分別）
% ---------------------------------------------------------------------
P_rx = zeros(N, num_samples);
P_rx_target = zeros(N, num_samples);
P_rx_water = zeros(N, num_samples);
P_rx_ballistic = zeros(N, num_samples);
for n = 1:N
P_rx_target(n, :) = causal_conv(x, h.target(n, :));
P_rx_water(n, :) = causal_conv(x, h.water(n, :));
P_rx_ballistic(n, :) = causal_conv(x, h.ballistic(n, :));
end
P_rx = P_rx_target + P_rx_water;

% ---------------------------------------------------------------------
% 4. 受光回路: 光電変換 (R * G) -> LPF -> HPF (零初期状態) -> 固定雑音
% ---------------------------------------------------------------------
rx_common = params.rx_common;
gain = rx_common.responsivity * rx_common.tia_gain;
v_signal = apply_rx_lowpass(gain * P_rx, rx_common.lpf_cutoff, fs);
v_signal = apply_rx_highpass_zero_state(v_signal, rx_common.hpf_cutoff, fs);
[y, noise_rms_v] = apply_receiver_noise(v_signal, params.noise_params, rx_common, fs);

% 受光光パワーが 1 dB 圧縮点を超える場合は警告（飽和特性は模擬しない）
P_rx_peak = max(P_rx, [], 2);
is_compressed = false(N, 1);
if rx_common.compression_power > 0
is_compressed = P_rx_peak > rx_common.compression_power;
for n = find(is_compressed(:).')
    warning('underwater_mc_lidar_simulation:compression', ...
        'Rx #%d: 受光光パワーのピーク %.3g W が 1 dB 圧縮点 %.3g W を超えています（飽和は模擬していません）', ...
        n, P_rx_peak(n), rx_common.compression_power);
end
end

% ---------------------------------------------------------------------
% 5. 解析情報の収集
% ---------------------------------------------------------------------
info.time_vector = t_vec;
info.h_total = h.target + h.water;           % インパルス応答 [1/ビン]（送信エネルギーに対する割合）
info.h_target = h.target;
info.h_water = h.water;
info.h_ballistic = h.ballistic;
info.P_rx = P_rx;                            % 受光光パワー [W]
info.P_rx_target = P_rx_target;
info.P_rx_water = P_rx_water;
info.P_rx_ballistic = P_rx_ballistic;
info.P_rx_peak = P_rx_peak;
info.is_compressed = is_compressed;
info.v_signal = v_signal;                    % 雑音付加前の受光モジュール出力 [V]
info.noise_rms_v = noise_rms_v;
if noise_rms_v > 0
info.snr_peak_db = 20 * log10(max(abs(v_signal), [], 2) / noise_rms_v);
else
info.snr_peak_db = Inf(N, 1);
end
info.mc = mc_stats;

% ---------------------------------------------------------------------
% 6. 可視化
% ---------------------------------------------------------------------
if params.enable_plot
render_mc_visuals(params, t_vec, x, info, y);
end

end

% =========================================================================
% 内部関数: デフォルトパラメータ（清澄な海水、1Tx x 1Rx、受光レンズ付き APD モジュール）
% =========================================================================
function p = get_default_params()
p.fs = 2.5e9;                                % サンプリング周波数 [Hz]
p.refractive_index = 1.33;                   % 水の屈折率（伝搬速度にのみ使用）
p.enable_plot = true;

% 水の固有光学特性 (清澄な海水: Jasman et al. 2018, Mobley 1994 / Petzold 1972 より)
p.water.a = 0.114;                           % 吸収係数 [1/m]
p.water.b = 0.037;                           % 散乱係数 [1/m]
p.water.g = 0.9;                             % HG 位相関数の非対称パラメータ（前方散乱が強い海水を想定した仮定値）

% 送信光源 (ガウスビーム、方位角 pointing_deg、発散角 FWHM [deg]（直交2軸それぞれ）)
p.tx.pos = [0.0, 0.50, 0.0];
p.tx.pointing_deg = 0.0;
p.tx.fwhm = 25.0;

% 受光器 (x = 0 面上、方位角 pointing_deg)
p.rx_array(1).pos = [0.0, -0.25, 0.0];
p.rx_array(1).pointing_deg = 0.0;

% 受光器共通仕様
p.rx_common.aperture_diameter = 25.4e-3;     % 受光開口径 [m]
p.rx_common.fov_deg = 1.0;                   % 視野角 半角 [deg]
p.rx_common.responsivity = 45.0;             % 感度 [A/W]
p.rx_common.tia_gain = 2.5e5 / 45.0;         % 変換ゲイン [V/A] (R * G = 2.5e5 V/W)
p.rx_common.hpf_cutoff = 0.05e6;             % 低域遮断周波数 [Hz] (0: 無効)
p.rx_common.lpf_cutoff = 1.0e9;              % 高域遮断周波数 [Hz] (0: 無効)
p.rx_common.compression_power = 0;           % 1 dB 圧縮点 [W] (0: 判定しない)

% 反射ターゲット (x = distance の平面上の矩形)
p.target.distance = 9.0;
p.target.center_yz = [0.0, 0.0];
p.target.width = 1.8;                        % y 方向 [m]
p.target.height = 1.8;                       % z 方向 [m]
p.target.reflectance = 0.8;

% モンテカルロ設定
p.mc.num_photons = 1e6;                      % 追跡する光子数
p.mc.batch_size = 2e5;                       % 1 バッチの光子数（メモリ使用量を制御）
p.mc.seed = [];                              % 乱数シード ([]: 固定しない)
p.mc.roulette_weight = 1e-3;                 % ロシアンルーレットを行う重みの閾値
p.mc.roulette_survival = 0.1;                % ロシアンルーレットの生存確率
p.mc.min_estimate_distance = 0.05;           % 局所推定の最小距離 [m]（受光器近傍の発散を抑制）
p.mc.max_events = 5000;                      % 1 光子あたりの最大イベント数（安全上限）

% 雑音 (光入力換算の固定雑音)
p.noise_params.enable = false;
p.noise_params.input_noise_rms = 16e-9;      % [W rms]
p.noise_params.seed = [];
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
elseif isstruct(def_s.(f)) && isstruct(user_s.(f)) && numel(def_s.(f)) == 1 && numel(user_s.(f)) == 1
    user_s.(f) = fill_missing_fields(user_s.(f), def_s.(f));
end
end
end

% =========================================================================
% 内部関数: 光子輸送（バッチ処理）によるインパルス応答の計算
%   h.target / h.water / h.ballistic: N x num_bins、送信エネルギーに対する割合 [1/ビン]
% =========================================================================
function [h, stats] = run_photon_transport(params, num_bins)
N = length(params.rx_array);
h.target = zeros(N, num_bins);
h.water = zeros(N, num_bins);
h.ballistic = zeros(N, num_bins);

if ~isempty(params.mc.seed)
rng(params.mc.seed);
end

num_total = round(params.mc.num_photons);
batch = max(1, round(params.mc.batch_size));
num_batches = ceil(num_total / batch);
batch_energy = zeros(num_batches, N);       % バッチごとの受光エネルギー（誤差推定用）
batch_energy_tgt = zeros(num_batches, N);   % 同上（ターゲット反射成分のみ）
num_events = 0;

for ib = 1:num_batches
nb = min(batch, num_total - (ib - 1) * batch);
[hb, ev] = transport_batch(params, nb, num_bins);
h.target = h.target + hb.target;
h.water = h.water + hb.water;
h.ballistic = h.ballistic + hb.ballistic;
batch_energy(ib, :) = (sum(hb.target, 2) + sum(hb.water, 2)).' / nb;
batch_energy_tgt(ib, :) = sum(hb.target, 2).' / nb;
num_events = num_events + ev;
end

h.target = h.target / num_total;
h.water = h.water / num_total;
h.ballistic = h.ballistic / num_total;

stats.num_photons = num_total;
stats.num_batches = num_batches;
stats.num_events = num_events;
stats.received_energy_fraction = (sum(h.target, 2) + sum(h.water, 2)).';
if num_batches > 1
% バッチ間のばらつきから求めた受光エネルギーの相対標準誤差（全体 / ターゲット反射成分）
stats.relative_std_error = std(batch_energy, 0, 1) / sqrt(num_batches) ./ max(mean(batch_energy, 1), realmin);
stats.relative_std_error_target = std(batch_energy_tgt, 0, 1) / sqrt(num_batches) ./ ...
    max(mean(batch_energy_tgt, 1), realmin);
else
stats.relative_std_error = NaN(1, N);
stats.relative_std_error_target = NaN(1, N);
end
end

% =========================================================================
% 内部関数: 1 バッチ分の光子輸送（全光子をベクトル化して同時に追跡）
% =========================================================================
function [hb, num_events] = transport_batch(params, nb, num_bins)
fs = params.fs;
v_light = 2.99792458e8 / params.refractive_index;
a = params.water.a;
b = params.water.b;
c = a + b;
albedo = b / c;
g = params.water.g;
t_max = num_bins / fs;
L_max = t_max * v_light;                     % 解析時間内に到達できる最大光路長 [m]

tgt = params.target;
D = tgt.distance;
rho = tgt.reflectance;

N = length(params.rx_array);
rx_pos = zeros(N, 3);
rx_axis = zeros(N, 3);
for n = 1:N
rx_pos(n, :) = params.rx_array(n).pos;
rx_axis(n, :) = [cosd(params.rx_array(n).pointing_deg), sind(params.rx_array(n).pointing_deg), 0];
end
A_rx = pi * (params.rx_common.aperture_diameter / 2)^2;
cos_fov = cosd(params.rx_common.fov_deg);
d_min = params.mc.min_estimate_distance;

hb.target = zeros(N, num_bins);
hb.water = zeros(N, num_bins);
hb.ballistic = zeros(N, num_bins);

% --- 光子の初期化: Tx 位置、ガウスビームの発散角で方向をサンプリング ---
sigma_rad = (params.tx.fwhm / 2) / sqrt(2 * log(2)) * (pi / 180);
th_y = sigma_rad * randn(nb, 1);
th_z = sigma_rad * randn(nb, 1);
dir_local = [ones(nb, 1), tan(th_y), tan(th_z)];
dir_local = dir_local ./ sqrt(sum(dir_local.^2, 2));
az = params.tx.pointing_deg;
R_az = [cosd(az), -sind(az), 0; sind(az), cosd(az), 0; 0, 0, 1];
dir = dir_local * R_az.';
pos = repmat(params.tx.pos, nb, 1);
w = ones(nb, 1);                             % 光子の重み（エネルギー）
L = zeros(nb, 1);                            % 累積光路長 [m]
n_scat = zeros(nb, 1);                       % 散乱回数
hit = false(nb, 1);                          % ターゲット反射を経たか

num_events = 0;
for it = 1:params.mc.max_events
if isempty(w)
    break;
end
num_events = num_events + numel(w);

% 自由行程（減衰係数 c の指数分布）
s = -log(rand(size(w))) / c;

% ターゲット平面 x = D との交差（+x 方向に進む光子）
t_plane = Inf(size(w));
fwd = dir(:, 1) > 0 & pos(:, 1) < D;
t_plane(fwd) = (D - pos(fwd, 1)) ./ dir(fwd, 1);
p_hit = pos + t_plane .* dir;
on_target = abs(p_hit(:, 2) - tgt.center_yz(1)) <= tgt.width / 2 & ...
            abs(p_hit(:, 3) - tgt.center_yz(2)) <= tgt.height / 2;
is_hit = fwd & on_target & t_plane <= s;

% ターゲット後方 (x > D) から裏面に当たる光子は吸収されたものとして終了
t_rear = Inf(size(w));
rear = dir(:, 1) < 0 & pos(:, 1) > D;
t_rear(rear) = (D - pos(rear, 1)) ./ dir(rear, 1);
p_rear = pos + min(t_rear, 1e6) .* dir;
on_rear = abs(p_rear(:, 2) - tgt.center_yz(1)) <= tgt.width / 2 & ...
          abs(p_rear(:, 3) - tgt.center_yz(2)) <= tgt.height / 2;
is_rear = rear & on_rear & t_rear <= s;

% 機器前面 x = 0 への到達（-x 方向に進む光子）は終了
t_back = Inf(size(w));
bwd = dir(:, 1) < 0;
t_back(bwd) = -pos(bwd, 1) ./ dir(bwd, 1);
is_out = ~is_hit & (is_rear | t_back <= s);

is_scat = ~is_hit & ~is_out;

% --- ターゲットでのランベルト反射 ---
ih = find(is_hit);
if ~isempty(ih)
    pos(ih, :) = p_hit(ih, :);
    L(ih) = L(ih) + t_plane(ih);
    n_tgt = repmat([-1, 0, 0], numel(ih), 1);
    for n = 1:N
        [E, t_arr, ok] = local_estimate(pos(ih, :), L(ih), rx_pos(n, :), rx_axis(n, :), ...
            A_rx, cos_fov, c, v_light, d_min, tgt, D);
        cos_out = max(sum(unit_to(pos(ih, :), rx_pos(n, :)) .* n_tgt, 2), 0);
        E = w(ih) .* rho .* (cos_out / pi) .* E;
        hb.target(n, :) = hb.target(n, :) + bin_add(t_arr(ok), E(ok), fs, num_bins);
        bal = ok & n_scat(ih) == 0;
        hb.ballistic(n, :) = hb.ballistic(n, :) + bin_add(t_arr(bal), E(bal), fs, num_bins);
    end
    w(ih) = w(ih) * rho;
    hit(ih) = true;
    dir(ih, :) = sample_lambertian(numel(ih));
end

% --- 水中での散乱 ---
is_ = find(is_scat);
if ~isempty(is_)
    pos(is_, :) = pos(is_, :) + s(is_) .* dir(is_, :);
    L(is_) = L(is_) + s(is_);
    for n = 1:N
        [E, t_arr, ok] = local_estimate(pos(is_, :), L(is_), rx_pos(n, :), rx_axis(n, :), ...
            A_rx, cos_fov, c, v_light, d_min, tgt, D);
        mu = sum(dir(is_, :) .* unit_to(pos(is_, :), rx_pos(n, :)), 2);
        E = w(is_) .* albedo .* hg_phase(mu, g) .* E;
        from_tgt = hit(is_);
        hb.target(n, :) = hb.target(n, :) + bin_add(t_arr(ok & from_tgt), E(ok & from_tgt), fs, num_bins);
        hb.water(n, :) = hb.water(n, :) + bin_add(t_arr(ok & ~from_tgt), E(ok & ~from_tgt), fs, num_bins);
    end
    w(is_) = w(is_) * albedo;
    n_scat(is_) = n_scat(is_) + 1;
    dir(is_, :) = sample_hg_direction(dir(is_, :), g);
end

% --- 終了判定: 機器前面への到達、解析時間超過、ロシアンルーレット ---
alive = ~is_out & L < L_max;
low = alive & w < params.mc.roulette_weight;
if any(low)
    survive = rand(sum(low), 1) < params.mc.roulette_survival;
    il = find(low);
    w(il(survive)) = w(il(survive)) / params.mc.roulette_survival;
    alive(il(~survive)) = false;
end
pos = pos(alive, :);
dir = dir(alive, :);
w = w(alive);
L = L(alive);
n_scat = n_scat(alive);
hit = hit(alive);
end
end

% =========================================================================
% 内部関数: 局所推定の幾何項（開口の立体角 × 直達透過率）と到達時刻
%   E   = A cos(psi) / d^2 * exp(-c d)（位相関数・反射率は呼び出し側で掛ける）
%   ok  = 視野角内、最小距離以上、かつターゲットに遮られない
% =========================================================================
function [E, t_arr, ok] = local_estimate(p, L, rx_p, rx_ax, A_rx, cos_fov, c, v_light, d_min, tgt, D)
q = p - rx_p;                                % 受光器から事象点へのベクトル
d = sqrt(sum(q.^2, 2));
cos_psi = (q * rx_ax.') ./ max(d, realmin);
ok = cos_psi >= cos_fov & d >= d_min;

% ターゲット後方 (x > D) の事象点は、受光器への直線がターゲットを横切る場合に遮蔽
behind = ok & p(:, 1) > D;
if any(behind)
    ib = find(behind);
    f = (p(ib, 1) - D) ./ (p(ib, 1) - rx_p(1));
    cross = p(ib, :) + f .* (rx_p - p(ib, :));
    blocked = abs(cross(:, 2) - tgt.center_yz(1)) <= tgt.width / 2 & ...
              abs(cross(:, 3) - tgt.center_yz(2)) <= tgt.height / 2;
    ok(ib(blocked)) = false;
end

E = A_rx * max(cos_psi, 0) ./ max(d, d_min).^2 .* exp(-c * d);
t_arr = (L + d) / v_light;
end

% =========================================================================
% 内部関数: 事象点から受光器への単位ベクトル
% =========================================================================
function u = unit_to(p, rx_p)
q = rx_p - p;
u = q ./ max(sqrt(sum(q.^2, 2)), realmin);
end

% =========================================================================
% 内部関数: Henyey-Greenstein 位相関数 [1/sr]（全立体角で積分すると 1）
% =========================================================================
function p = hg_phase(mu, g)
p = (1 - g^2) ./ (4 * pi * (1 + g^2 - 2 * g * mu).^1.5);
end

% =========================================================================
% 内部関数: HG 位相関数に従う散乱後の方向のサンプリング
% =========================================================================
function d_new = sample_hg_direction(d, g)
m = size(d, 1);
xi = rand(m, 1);
if abs(g) < 1e-6
mu = 2 * xi - 1;
else
mu = (1 + g^2 - ((1 - g^2) ./ (1 - g + 2 * g * xi)).^2) / (2 * g);
end
mu = max(min(mu, 1), -1);
sin_t = sqrt(1 - mu.^2);
phi = 2 * pi * rand(m, 1);
d_new = rotate_direction(d, mu, sin_t, cos(phi), sin(phi));
end

% =========================================================================
% 内部関数: 進行方向 d から極角 (mu = cos(theta)) と方位角 phi だけ回転した方向
% =========================================================================
function d_new = rotate_direction(d, mu, sin_t, cos_p, sin_p)
ux = d(:, 1);
uy = d(:, 2);
uz = d(:, 3);
d_new = zeros(size(d));
pole = abs(uz) > 0.99999;
r = sqrt(max(1 - uz.^2, 0));
np = ~pole;
d_new(np, 1) = sin_t(np) .* (ux(np) .* uz(np) .* cos_p(np) - uy(np) .* sin_p(np)) ./ r(np) + ux(np) .* mu(np);
d_new(np, 2) = sin_t(np) .* (uy(np) .* uz(np) .* cos_p(np) + ux(np) .* sin_p(np)) ./ r(np) + uy(np) .* mu(np);
d_new(np, 3) = -sin_t(np) .* cos_p(np) .* r(np) + uz(np) .* mu(np);
d_new(pole, 1) = sin_t(pole) .* cos_p(pole);
d_new(pole, 2) = sin_t(pole) .* sin_p(pole);
d_new(pole, 3) = sign(uz(pole)) .* mu(pole);
end

% =========================================================================
% 内部関数: ターゲット面（法線 -x）からのランベルト反射方向のサンプリング
% =========================================================================
function d_new = sample_lambertian(m)
cos_t = sqrt(rand(m, 1));
sin_t = sqrt(1 - cos_t.^2);
phi = 2 * pi * rand(m, 1);
d_new = [-cos_t, sin_t .* cos(phi), sin_t .* sin(phi)];
end

% =========================================================================
% 内部関数: 到達時刻 t のエネルギー E を時間ビンに加算
% =========================================================================
function hist = bin_add(t, E, fs, num_bins)
hist = zeros(1, num_bins);
if isempty(t)
return;
end
k = floor(t * fs) + 1;
valid = k >= 1 & k <= num_bins;
if any(valid)
hist = accumarray(k(valid), E(valid), [num_bins, 1]).';
end
end

% =========================================================================
% 内部関数: 因果的な畳み込み（出力長は入力と同じ）
% =========================================================================
function y = causal_conv(x, h)
y = conv(x, h);
y = y(1:length(x));
end

% =========================================================================
% 内部関数: 受光回路 1次LPF（周波数領域、ゼロパディングによる線形畳み込み）
% =========================================================================
function y = apply_rx_lowpass(x, f_H, fs)
y = x;
if ~(f_H > 0) || ~isfinite(f_H)
return;
end
[N, num_samples] = size(x);
n_fft = num_samples + 256;
f = (0:(n_fft - 1)) * (fs / n_fft);
f(f > fs / 2) = f(f > fs / 2) - fs;          % 負の周波数
H = 1 ./ (1 + 1i * f / f_H);
Y = real(ifft(fft([x, zeros(N, n_fft - num_samples)], [], 2) .* repmat(H, N, 1), [], 2));
y = Y(:, 1:num_samples);
end

% =========================================================================
% 内部関数: 受光回路 1次HPF（AC結合）の零初期状態応答（双一次変換による IIR）
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
% 内部関数: 光入力換算の固定雑音の付加（受光回路の帯域で整形し、出力換算 rms に合わせる）
% =========================================================================
function [y_out, noise_rms_v] = apply_receiver_noise(y_in, noise_params, rx_common, fs)
y_out = y_in;
noise_rms_v = 0;
if ~noise_params.enable
return;
end
if ~isempty(noise_params.seed)
rng(noise_params.seed);
end
[N, num_samples] = size(y_in);
noise_rms_v = noise_params.input_noise_rms * rx_common.responsivity * rx_common.tia_gain;
f = (0:(num_samples - 1)) * (fs / num_samples);
f(f > fs / 2) = f(f > fs / 2) - fs;
H = ones(1, num_samples);
if rx_common.lpf_cutoff > 0 && isfinite(rx_common.lpf_cutoff)
H = H ./ (1 + 1i * f / rx_common.lpf_cutoff);
end
if rx_common.hpf_cutoff > 0
jf = 1i * f / rx_common.hpf_cutoff;
H = H .* (jf ./ (1 + jf));
end
shaped = real(ifft(fft(randn(N, num_samples), [], 2) .* repmat(H, N, 1), [], 2));
y_out = y_in + shaped / sqrt(mean(abs(H).^2)) * noise_rms_v;
end

% =========================================================================
% 内部関数: 可視化（インパルス応答の成分、送受信波形）
% =========================================================================
function render_mc_visuals(params, t, x, info, y)
N = size(info.P_rx, 1);
t_ns = t * 1e9;
colors = lines(4);
w = params.water;

% Figure 1: Rx ごとのインパルス応答（成分別、対数軸）
figure('Name', 'Underwater MC LiDAR - Impulse Response', 'NumberTitle', 'off', ...
'Position', [60, 80, 760, 260 * N + 120]);
for n = 1:N
subplot(N, 1, n);
h_all = {info.h_total(n, :), info.h_target(n, :), info.h_water(n, :), info.h_ballistic(n, :)};
names = {'Total', 'Target Return', 'Water Backscatter', 'Ballistic Target (No Scattering)'};
styles = {'-', '-', '-', '--'};
h_pos = [h_all{:}];
h_pos = h_pos(h_pos > 0);
for k = 1:4
    hk = h_all{k};
    hk(hk <= 0) = NaN;
    semilogy(t_ns, hk, styles{k}, 'Color', colors(k, :), 'LineWidth', 1.3, 'DisplayName', names{k});
    hold on;
end
grid on;
if ~isempty(h_pos)
    ylim([max(min(h_pos), max(h_pos) * 1e-12), max(h_pos) * 3]);
end
xlim([0, max(t_ns)]);
xlabel('Time [ns]');
ylabel('h(t) [1/bin]');
title(sprintf('Rx #%d Impulse Response (a = %.3g, b = %.3g /m, g = %.3g, %.2g photons)', ...
    n, w.a, w.b, w.g, info.mc.num_photons));
legend('Location', 'northeastoutside');
end

% Figure 2: 送信光パワー、受光光パワー（成分別）、受光モジュール出力
figure('Name', 'Underwater MC LiDAR - Waveforms', 'NumberTitle', 'off', ...
'Position', [840, 80, 760, 720]);
subplot(3, 1, 1);
plot(t_ns, x, 'LineWidth', 1.5);
grid on;
xlim([0, max(t_ns)]);
xlabel('Time [ns]');
ylabel('Optical Power [W]');
title('Transmitted Optical Power P_{tx}(t)');

subplot(3, 1, 2);
p_all = {info.P_rx(1, :), info.P_rx_target(1, :), info.P_rx_water(1, :), info.P_rx_ballistic(1, :)};
names = {'Total', 'Target Return', 'Water Backscatter', 'Ballistic Target'};
styles = {'-', '-', '-', '--'};
p_pos = [p_all{:}];
p_pos = p_pos(p_pos > 0);
for k = 1:4
pk = p_all{k};
pk(pk <= 0) = NaN;
semilogy(t_ns, pk, styles{k}, 'Color', colors(k, :), 'LineWidth', 1.3, 'DisplayName', names{k});
hold on;
end
if params.noise_params.enable
semilogy([0, max(t_ns)], [1, 1] * params.noise_params.input_noise_rms, 'k:', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Noise Level (%.3g nW rms)', params.noise_params.input_noise_rms * 1e9));
end
grid on;
if ~isempty(p_pos)
ylim([max(min(p_pos), max(p_pos) * 1e-9), max(p_pos) * 3]);
end
xlim([0, max(t_ns)]);
xlabel('Time [ns]');
ylabel('Optical Power [W]');
title('Rx #1 Received Optical Power P_{rx}(t) (Photoreceiver Input, Log Scale)');
legend('Location', 'northeastoutside');

subplot(3, 1, 3);
plot(t_ns, y(1, :), 'LineWidth', 1.2, 'DisplayName', 'Output with Noise');
hold on;
plot(t_ns, info.v_signal(1, :), 'k--', 'LineWidth', 1.0, 'DisplayName', 'Noise-Free');
grid on;
xlim([0, max(t_ns)]);
xlabel('Time [ns]');
ylabel('Output Voltage [V]');
title('Rx #1 Photoreceiver Output Voltage v_{out}(t)');
legend('Location', 'northeastoutside');
end
