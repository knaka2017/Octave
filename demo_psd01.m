% ============ デバッグ版：PSD計算の各ステップを表示 ============

% パラメータ設定
fs = 100000;    % サンプリング周波数 100kHz
N = 2^18;      % サンプル数
f0 = 1000;      % 信号周波数 1kHz
A = 1.0;        % 振幅 1V
Z0 = 50;        % インピーダンス 50Ω

% 時間軸の生成
t = (0:N-1)'/fs;
x = A * sin(2*pi*f0*t);

% 窓関数
window = hanning(N);
x_windowed = x .* window;

% FFT
X = fft(x_windowed);

% 基本量の計算
fprintf('========== デバッグ情報 ==========\n');
fprintf('N = %d, fs = %d Hz, df = %.6f Hz\n', N, fs, fs/N);
fprintf('A = %.1f V, Z0 = %.0f Ω\n', A, Z0);

% 総電力（理論値）
P_total_W = A^2/(2*Z0);  % = 0.01 W
fprintf('\n1. 総電力（理論値）: %.6f W\n', P_total_W);

% 窓関数のエネルギー
window_energy = sum(window.^2)/N;
fprintf('\n2. 窓関数のエネルギー: %.6f\n', window_energy);
fprintf('   補正係数 (1/window_energy): %.6f\n', 1/window_energy);

% ピーク周波数のインデックス
[~, idx] = max(abs(X));
fprintf('\n3. FFTピーク: インデックス %d, 振幅 |X| = %.6f\n', idx, abs(X(idx)));

% 電力スペクトル（両側、補正前）
P_two_sided_uncorrected = abs(X(idx)).^2 / (fs * N * Z0);
fprintf('\n4. 両側PSD（補正前）: %.10f W/Hz\n', P_two_sided_uncorrected);

% 窓関数補正後
P_two_sided_corrected = P_two_sided_uncorrected / window_energy;
fprintf('5. 両側PSD（補正後）: %.10f W/Hz\n', P_two_sided_corrected);

% 片側スペクトルへの変換
if mod(N,2)==0
    P_one_sided = P_two_sided_corrected;
    if idx > 1 && idx < N/2+1
        P_one_sided = P_two_sided_corrected * 2;
        fprintf('6. 片側PSD変換（2倍適用）: %.10f W/Hz\n', P_one_sided);
    else
        fprintf('6. 片側PSD変換（DC/ナイキスト、2倍なし）: %.10f W/Hz\n', P_one_sided);
    end
end

% 理論PSD値
psd_theoretical = P_total_W * N / fs / window_energy;
fprintf('\n7. 理論PSD値: %.10f W/Hz\n', psd_theoretical);
fprintf('   計算値/理論値: %.6f\n', P_one_sided/psd_theoretical);

% ============ 修正版関数（複数のバリエーション） ============

function [pxx, f] = psd_v1(x, fs, unit)
    % バージョン1: 元の実装
    Z0 = 50;
    N = length(x);
    window = hanning(N);
    xw = x(:) .* window;
    X = fft(xw);

    % 窓関数補正なし
    pxx_two = (abs(X).^2) / (fs * N * Z0);

    if mod(N,2)==0
        pxx = pxx_two(1:N/2+1);
        pxx(2:end-1) = 2 * pxx(2:end-1);
        f = (0:N/2)' * fs / N;
    else
        pxx = pxx_two(1:(N+1)/2);
        pxx(2:end) = 2 * pxx(2:end);
        f = (0:(N-1)/2)' * fs / N;
    end

    if unit==0, pxx = 10*log10(pxx)+30; end
end

function [pxx, f] = psd_v2(x, fs, unit)
    % バージョン2: 窓関数補正あり
    Z0 = 50;
    N = length(x);
    window = hanning(N);
    xw = x(:) .* window;
    X = fft(xw);

    window_energy = sum(window.^2)/N;

    % 窓関数補正あり
    pxx_two = (abs(X).^2) / (fs * N * Z0 * window_energy);

    if mod(N,2)==0
        pxx = pxx_two(1:N/2+1);
        pxx(2:end-1) = 2 * pxx(2:end-1);
        f = (0:N/2)' * fs / N;
    else
        pxx = pxx_two(1:(N+1)/2);
        pxx(2:end) = 2 * pxx(2:end);
        f = (0:(N-1)/2)' * fs / N;
    end

    if unit==0, pxx = 10*log10(pxx)+30; end
end

function [pxx, f] = psd_v3(x, fs, unit)
    % バージョン3: 片側変換を2回適用しない
    Z0 = 50;
    N = length(x);
    window = hanning(N);
    xw = x(:) .* window;
    X = fft(xw);

    window_energy = sum(window.^2)/N;

    % 両側PSD
    pxx_two = (abs(X).^2) / (fs * N * Z0 * window_energy);

    if mod(N,2)==0
        pxx = pxx_two(1:N/2+1);
        % 2倍しない（すでに両側PSDとして正しい）
        f = (0:N/2)' * fs / N;
    else
        pxx = pxx_two(1:(N+1)/2);
        f = (0:(N-1)/2)' * fs / N;
    end

    if unit==0, pxx = 10*log10(pxx)+30; end
end

% ============ 各バージョンのテスト ============
fprintf('\n========== 各バージョンの比較 ==========\n');

[pxx1, f1] = psd_v1(x, fs, 1);
[~, idx1] = max(pxx1);
fprintf('v1 (補正なし): %.10f W/Hz\n', pxx1(idx1));

[pxx2, f2] = psd_v2(x, fs, 1);
[~, idx2] = max(pxx2);
fprintf('v2 (補正あり): %.10f W/Hz\n', pxx2(idx2));

[pxx3, f3] = psd_v3(x, fs, 1);
[~, idx3] = max(pxx3);
fprintf('v3 (補正あり、2倍なし): %.10f W/Hz\n', pxx3(idx3));

fprintf('理論値: %.10f W/Hz\n', psd_theoretical);

% ============ 正しい実装 ============
function [pxx, f] = power_spectral_density_correct(x, fs, unit)
% 正しいPSD計算

    if nargin < 2
        error('信号とサンプリング周波数が必要です');
    end
    if nargin < 3 || isempty(unit)
        unit = 0;
    end

    Z0 = 50;
    N = length(x);

    % 窓関数
    window = hanning(N);
    x_windowed = x(:) .* window;

    % FFT
    X = fft(x_windowed);

    % 窓関数のエネルギー補正
    window_energy = sum(window.^2)/N;

    % 両側PSD [W/Hz]（電圧→電力変換、窓関数補正済み）
    psd_two_sided = (abs(X).^2) / (fs * N * Z0 * window_energy);

    % 片側スペクトル（電力保存のため2倍）
    if mod(N, 2) == 0
        % Nが偶数の場合
        pxx = psd_two_sided(1:N/2+1);
        pxx(2:end-1) = 2 * pxx(2:end-1);  % 直流とナイキスト以外を2倍
        f = (0:N/2)' * fs / N;
    else
        % Nが奇数の場合
        pxx = psd_two_sided(1:(N+1)/2);
        pxx(2:end) = 2 * pxx(2:end);  % 直流以外を2倍
        f = (0:(N-1)/2)' * fs / N;
    end

    % 単位変換
    switch unit
        case 0  % dBm
            pxx = 10 * log10(pxx) + 30;
        case 1  % W
            % そのまま
        case 2  % dBW
            pxx = 10 * log10(pxx);
    end
end

% ============ 最終テスト ============
[pxx_final, f_final] = power_spectral_density_correct(x, fs, 1);
[~, idx_final] = max(pxx_final);

fprintf('\n========== 最終結果 ==========\n');
fprintf('計算値: %.10f W/Hz @ %.2f Hz\n', pxx_final(idx_final), f_final(idx_final));
fprintf('理論値: %.10f W/Hz\n', psd_theoretical);
fprintf('誤差: %.6f %%\n', abs(pxx_final(idx_final)-psd_theoretical)/psd_theoretical*100);
