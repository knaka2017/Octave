pkg load signal;

% --- Data Initialization ---
if !exist('x', 'var')
    printf("Test mode: Generating 10M sample signal...\n");
    Fs = 100e6; N = 1e7; dt = 1/Fs;
    t = (0:N-1)' * dt;
    Fc_ideal = 1.23456e6;
    x = sin(2 * pi * Fc_ideal * t) + 0.03 * randn(size(t));
end

% --- Parameters ---
N = length(x);
dt = t(2) - t(1);
Fs = 1 / dt;
margin = 5000;

printf("\nAnalysis Started (Savitzky-Golay Method): N = %d, Fs = %.2f MHz\n", N, Fs/1e6);

% ---------------------------------------------------------
% 1. Analytic Signal (Hilbert Transform)
% ---------------------------------------------------------
printf("Step 1/5: Computing Hilbert Transform...\n");
tic;
z = hilbert(x);
I = real(z);
Q = imag(z);
printf("   -> Completed in %.2f seconds.\n", toc);

% ---------------------------------------------------------
% 2. Robust Derivative via Savitzky-Golay (No 'diff' used)
% ---------------------------------------------------------
printf("Step 2/5: Calculating Derivatives via Savitzky-Golay filter...\n");

% Window size for smoothing (must be odd)
% Approx 1/4 cycle is usually good for balance between noise and signal tracking
p_order = 2; % Polynomial order
w_len = floor(Fs / (Fc_ideal * 4));
if mod(w_len, 2) == 0; w_len++; end
if w_len < 5; w_len = 5; end

% sgolay() returns differentiation filters
[~, g] = sgolay(p_order, w_len);

% Differentiation filter (1st derivative is the 2nd column of g)
% Note: sgolay differentiation filter gives (1/dt^k) * derivative
h_diff = g(:, 2)' / dt;

% Apply filter (conv or filter)
dI = conv(I, h_diff, 'same');
dQ = conv(Q, h_diff, 'same');

% ---------------------------------------------------------
% 3. Instantaneous Frequency (Cross-Product Method)
% ---------------------------------------------------------
printf("Step 3/5: Computing Inst. Frequency (No-Phase, No-Diff)...\n");
% f = (1/2pi) * (I*dQ - Q*dI) / (I^2 + Q^2)
inst_freq_raw = (1 / (2 * pi)) * (I .* dQ - Q .* dI) ./ (I.^2 + Q.^2);

% ---------------------------------------------------------
% 4. Quantitative Metrics (using Trimmed data)
% ---------------------------------------------------------
printf("Step 4/5: Calculating Quantitative Metrics...\n");

total_margin = margin + w_len;
v_idx = (total_margin + 1) : (N - total_margin);

t_v = t(v_idx);
inst_freq = inst_freq_raw(v_idx);
envelope = abs(z(v_idx));

% 1. Frequency
avg_freq = mean(inst_freq);
% 2. Frequency Stability
freq_std = std(inst_freq);

% 3. Phase Stability (Still using unwrap here for final report only)
inst_phase_v = unwrap(angle(z(v_idx)));
p = polyfit(t_v, inst_phase_v, 1);
phase_residual = inst_phase_v - polyval(p, t_v);
phase_error_rms = std(phase_residual);

% 5. Anomaly Detection
env_mean = mean(envelope);
env_std = std(envelope);
anomaly_indices = find(abs(envelope - env_mean) > 6 * env_std);

% ---------------------------------------------------------
% 5. Output and Visualization
% ---------------------------------------------------------
fprintf('\n================== ANALYSIS REPORT ==================\n');
fprintf('1. Estimated Mean Frequency:    %.6f MHz\n', avg_freq/1e6);
fprintf('2. Frequency Stability (Std):   %.3f Hz (SG-filtered)\n', freq_std);
fprintf('3. Phase Stability (RMS Error): %.6f rad\n', phase_error_rms);
fprintf('4. Inst. Freq Range:            [%.3f, %.3f] MHz\n', min(inst_freq)/1e6, max(inst_freq)/1e6);
fprintf('5. Anomalies Detected:          %d points\n', length(anomaly_indices));
fprintf('=====================================================\n');

printf("Step 5/5: Generating Plots...\n");
plot_step = max(1, floor(length(t_v) / 50000));
t_plot = t_v(1:plot_step:end);

figure('Name', 'Robust Waveform Analysis (Savitzky-Golay)');
subplot(3,1,1);
plot(t_plot, inst_freq(1:plot_step:end)/1e6); grid on;
title('4. Instantaneous Frequency (Savitzky-Golay Derivative)');
ylabel('Frequency [MHz]');

subplot(3,1,2);
plot(t_plot, phase_residual(1:plot_step:end)); grid on;
title('3. Phase Stability (Residual Phase Error)');
ylabel('Phase Error [rad]');

subplot(3,1,3);
plot(t_plot, envelope(1:plot_step:end)); grid on;
title('5. Anomaly Detection (Envelope)');
ylabel('Amplitude');

printf("All processes completed.\n");

