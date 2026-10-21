function [X, f] = compute_spectrum(x, Fs, use_dBm, Z)
% GET_SPECTRUM  離散時間信号の周波数スペクトルを計算する。
%
% [X, f] = compute_spectrum(x, Fs, use_dBm, Z)
%
% 入力:
%   x       : 離散時間信号ベクトル (N x 1, 単位: V)
%   Fs      : サンプリング周波数 (Hz)
%   use_dBm : スペクトルを dBm 単位で出力するかどうかのフラグ (true/false)
%   Z       : インピーダンス (オーム)。use_dBmがtrueの場合にのみ使用。
%
% 出力:
%   X       : 周波数スペクトルベクトル (単位: V_rms/sqrt(Hz) または dBm)
%   f       : 対応する周波数ベクトル (Hz)

    N = length(x); % 信号のサンプル数
    T = 1/Fs;      % サンプリング周期

    % 1. FFTを実行
    Y = fft(x);

    % 2. 片側スペクトルと周波数軸の作成
    % Nが偶数の場合
    if mod(N, 2) == 0
        L = N/2 + 1;
        % FFT結果の半分（DC成分からナイキスト周波数まで）
        Y_one_sided = Y(1:L);
        % 周波数ベクトル
        f = (0:(L-1)) * Fs/N;
        % Y(1)とY(L)（DCとナイキスト）以外を2倍して片側スペクトルを構成
        Y_one_sided(2:(L-1)) = 2 * Y_one_sided(2:(L-1));
    % Nが奇数の場合
    else
        L = (N+1)/2;
        % FFT結果の半分（DC成分からナイキスト周波数以下まで）
        Y_one_sided = Y(1:L);
        % 周波数ベクトル
        f = (0:(L-1)) * Fs/N;
        % Y(1)（DC）以外を2倍して片側スペクトルを構成
        Y_one_sided(2:L) = 2 * Y_one_sided(2:L);
    end

    % 3. スケーリングの適用
    % V_rms 単位への変換
    X_Vrms = abs(Y_one_sided) / N;

    % 4. 出力単位の決定
    if use_dBm
        % V_rms を dBm に変換: P_dBm = 10 * log10(P_mW)
        % P = V_rms^2 / Z
        % P_mW = (V_rms^2 / Z) * 1000

        if nargin < 4 || isempty(Z)
            error('dBm出力にはインピーダンス Z が必要です。');
        end

        % V_rmsを2乗して電力P[W]を計算
        P_W = (X_Vrms .^ 2) / Z;

        % WをmWに変換
        P_mW = P_W * 1000;

        % dBmに変換
        X = 10 * log10(P_mW);

    else
        % デフォルトは振幅スペクトル (V_rms)
        X = X_Vrms;
    end
end
