% SiPM Pulse & Pole-Zero (PZ) Filter Simulation
% Constraints: fs=2.5GHz, BW=500MHz

clear; clf;
pkg load signal; % butter, filter, freqz関数のために必要

% --- 1. パラメータ設定 ---
fs = 2.5e9;          % サンプリング周波数: 2.5 GHz
ts = 1/fs;           % サンプリング周期: 0.4 ns
fc_bw = 500e6;       % アナログ帯域制限: 500 MHz
T_total = 150e-9;    % シミュレーション時間: 150 ns
t = 0:ts:T_total;

% --- 2. SiPMパルス生成 (二重指数関数モデル) ---
tau_rise = 1.0e-9;   % 立ち上がり時定数 (1 ns)
tau_orig = 50.0e-9;  % 元のテールの時定数 (50 ns)
% 理想パルス
raw_pulse = (exp(-t/tau_orig) - exp(-t/tau_rise));

% 500MHzの帯域制限（1次バターワース）を適用
[b_lp, a_lp] = butter(1, fc_bw/(fs/2));
input_pulse = filter(b_lp, a_lp, raw_pulse);
input_pulse = input_pulse / max(input_pulse); % 正規化

% --- 3. PZフィルターの設計 ---
tau_new = 5.0e-9;    % 目標の時定数 (5 ns)

% デジタルフィルター係数の算出
alpha = exp(-ts/tau_orig); % 零点 (Zero): 元のテールを消す
beta  = exp(-ts/tau_new);  % 極 (Pole): 新しいテールを作る

% 伝達関数: H(z) = (1 - alpha*z^-1) / (1 - beta*z^-1)
b_pz = [1, -alpha];
a_pz = [1, -beta];

% フィルタリング実行
output_pulse = filter(b_pz, a_pz, input_pulse);

% --- 4. テキスト出力 (コンソール) ---
fprintf('\n--- PZ Filter Simulation Results ---\n');
fprintf('Sampling Frequency : %.1f GHz\n', fs/1e9);
fprintf('Sampling Interval  : %.2f ns\n', ts*1e9);
fprintf('Analog Bandwidth   : %.0f MHz\n', fc_bw/1e6);
fprintf('Original Tail (tau): %.1f ns\n', tau_orig*1e9);
fprintf('Target Tail (tau)  : %.1f ns\n', tau_new*1e9);
fprintf('------------------------------------\n');
fprintf('Filter Coefficient Alpha (Zero): %.5f\n', alpha);
fprintf('Filter Coefficient Beta  (Pole): %.5f\n', beta);
fprintf('DC Gain (Low Freq) : %.4f (%.2f dB)\n', ...
        (1-alpha)/(1-beta), 20*log10((1-alpha)/(1-beta)));
fprintf('------------------------------------\n');

% --- 5. グラフ描画 ---
subplot(2,1,1);
plot(t*1e9, input_pulse, 'b', 'LineWidth', 1.5, 'DisplayName', 'Original (50ns tail)');
hold on;
plot(t*1e9, output_pulse, 'r', 'LineWidth', 1.5, 'DisplayName', 'PZ Corrected (5ns tail)');
grid on;
title(['SiPM Pulse Compensation (fs=2.5GHz, BW=500MHz)']);
xlabel('Time [ns]');
ylabel('Normalized Amplitude');
legend();
xlim([0 100]);

subplot(2,1,2);
[h, w] = freqz(b_pz, a_pz, 1024, fs);
pos_idx = w > 0;
semilogx(w(pos_idx), 20*log10(abs(h(pos_idx))), 'k', 'LineWidth', 1.5);
grid on;
title('PZ Filter Frequency Response');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xlim([1e5 fs/2]); % 100kHzからナイキスト周波数まで表示

