function [pxx, f] = power_spectral_density(x, fs, unit)
% [pxx, f] = power_spectral_density_fixed(x, fs, unit)
%  50 ohm load
%  unit:
%    0: dBm
%    1: W
%    2: dBW

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

    % 窓関数のエネルギー補正係数
    window_energy = sum(window.^2)/N;  % 平均2乗値
    energy_correction = 1/window_energy;

    % 両側PSD [W/Hz]（電圧→電力変換、窓関数補正済み）
    psd_two_sided = (abs(X).^2) * energy_correction / (fs * N * Z0);

    % 片側スペクトルへの変換
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
        otherwise
            error('単位は0(dBm), 1(W), 2(dBW)のみ');
    end
end

