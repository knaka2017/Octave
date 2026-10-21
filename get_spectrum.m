function [X, f] = get_spectrum(x, Fs, use_dBm, Z)
% GET_SPECTRUM  離散時間信号の周波数スペクトルを計算する。
%
% [X, f] = get_spectrum(x, Fs, use_dBm, Z)
%
% 入力:
%   x       : 離散時間信号ベクトル (N x 1, 単位: V)
%   Fs      : サンプリング周波数 (Hz)
%   use_dBm : スペクトルを dBm 単位で出力するかどうかのフラグ (true/false)
%   Z       : インピーダンス (オーム)。use_dBmがtrueの場合にのみ使用。
%
% 出力:
%   X       : 周波数スペクトルベクトル (単位: V_peak/sqrt(Hz) または dBm)
%   f       : 対応する周波数ベクトル (Hz)

    N = length(x); % 信号のサンプル数

    % 1. FFTを実行
    Y = fft(x);

    % 2. 片側スペクトルと周波数軸の作成
    % Nが偶数の場合
    if mod(N, 2) == 0
        L = N/2 + 1;
        Y_one_sided = Y(1:L);
        f = (0:(L-1)) * Fs/N;
        % DC(Y(1))とナイキスト周波数(Y(L))以外を2倍
        Y_one_sided(2:(L-1)) = 2 * Y_one_sided(2:(L-1));
    % Nが奇数の場合
    else
        L = (N+1)/2;
        Y_one_sided = Y(1:L);
        f = (0:(L-1)) * Fs/N;
        % DC(Y(1))以外を2倍
        Y_one_sided(2:L) = 2 * Y_one_sided(2:L);
    end

    % 3. スケーリングの適用
    % X_Vpeak は信号のV_peak振幅に対応 (FFTの定義に基づく)
    X_Vpeak = abs(Y_one_sided) / N;

    % 4. 出力単位の決定
    if use_dBm
        % V_peakをV_rmsに変換し、真の電力を計算する

        if nargin < 4 || isempty(Z)
            error('dBm出力にはインピーダンス Z が必要です。');
        end

        % V_rms = V_peak / sqrt(2)
        V_rms = X_Vpeak / sqrt(2);

        % 電力P[W]を計算: P = V_rms^2 / Z
        P_W = (V_rms .^ 2) / Z;

        % WをmWに変換
        P_mW = P_W * 1000;

        % dBmに変換: P_dBm = 10 * log10(P_mW)
        X = 10 * log10(P_mW);

    else
        % dBmフラグがfalseの場合は、ピーク振幅(V_peak)を出力
        X = X_Vpeak;
    end
end
