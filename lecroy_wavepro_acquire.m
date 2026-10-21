function [waveform, actual_params] = lecroy_wavepro_acquire(ip_address, channel, sample_rate_request, record_length_request)
% LECROY_WAVEPRO_ACQUIRE LeCroy WavePro 960オシロスコープから波形データを取得
%
% [waveform, actual_params] = lecroy_wavepro_acquire(ip_address, channel, sample_rate_request, record_length_request)
%
% LeCroy WavePro 960 (DDA260) オシロスコープから波形データを取得します。
% TCP/IP (VICP プロトコル、ポート1861) を使用し、IEEE 488.2 SCPI コマンドで制御します。
% 元のオシロスコープ設定は自動的に保存・復元されます。
%
% 入力引数:
%   ip_address           - オシロスコープのIPアドレス (文字列)
%                          例: '192.168.1.100'
%   channel              - 取得チャンネル (文字列)
%                          例: 'C1', 'C2', 'C3', 'C4', 'M1', 'M2'
%   sample_rate_request  - 要求サンプリングレート (数値または文字列)
%                          例: 16e9, '16G', '800M', 1e6
%   record_length_request - 要求レコード長、サンプル数 (数値)
%                          例: 10000, 1000000
%                          注意: オシロスコープは離散的な値のみ許可します。
%                          指定値は最も近い許可値に自動調整されます。
%                          実際の値は actual_params.RecordLength で確認できます。
%
% 出力引数:
%   waveform            - 波形データ構造体
%       .X              - 時間軸ベクトル [秒]
%       .Y              - 電圧ベクトル [ボルト]
%   actual_params       - 実際の取得パラメータ構造体
%       .SamplingRate   - 実サンプリングレート [Hz]
%       .RecordLength   - 実レコード長 [サンプル数]
%       .TimePerDiv     - 時間軸ディビジョン [秒/div]
%       .V_Gain         - 垂直軸ゲイン
%       .V_Offset       - 垂直軸オフセット [V]
%       .T0             - 時間軸オフセット [秒]
%
% 使用例:
%   % 基本的な使用方法
%   [wave, params] = lecroy_wavepro_acquire('192.168.1.100', 'C1', '16G', 10000);
%   plot(wave.X, wave.Y);
%   xlabel('Time [s]'); ylabel('Voltage [V]');
%
%   % 数値指定
%   [wave, params] = lecroy_wavepro_acquire('192.168.1.100', 'C2', 1e9, 50000);
%
% 必要なパッケージ:
%   Octave instrument-control パッケージ (pkg load instrument-control)
%
% 注意事項:
%   - 最大サンプリングレートは 16 GS/s に制限されます
%   - エラー発生時は空の構造体が返されます
%   - オシロスコープの元の設定は自動的に復元されます
%
% 参考文献:
%   LeCroy WavePro 900 Series Operator's Manual
%   LeCroy Remote Control Manual

    % --- I. 初期設定と変数定義 ---
    VICP_PORT = 1861;
    MAX_SARA = 16e9; % 最大サンプリングレート 16 GS/s
    WAVEDESC_SIZE = 346; % LeCroy WAVEDESCの標準サイズ
    t_scope = []; % TCPクライアントオブジェクト（空配列で初期化）
    original_tdiv_str = '';
    original_trmd_str = '';

    % 入力レートの処理とクリッピング
    if ischar(sample_rate_request)
        % 例: '16G', '800M' などの文字列を数値に変換
        sample_rate_request = parse_engineering_notation(sample_rate_request);
    endif

    S_user = min(sample_rate_request, MAX_SARA); % 16 GS/sに上限をクリップ

    % --- II. 状態保存と接続の try-catch ブロック ---
    try
        % 1. 接続確立
        t_scope = tcpclient(ip_address, VICP_PORT, 'Timeout', 10);

        % 終端文字の設定 (LeCroy SCPIは通常ラインフィード \n)
        configureTerminator(t_scope, 'LF');

        % 2. 既存の状態の保存

        % タイムベース設定 (TIM:MAIN:SCAL?) の保存
        write(t_scope, 'TIM:MAIN:SCAL?', 'ascii');
        original_tdiv_str = strtrim(char(read(t_scope, 100)));

%        % タイムベース設定 (TDIV) の保存
%        write(t_scope, 'TDIV?', 'ascii');
%        original_tdiv_str = strtrim(char(read(t_scope, 100)));

        % トリガモード (TRMD) の保存
        write(t_scope, 'TRMD?', 'ascii');
        original_trmd_str = strtrim(char(read(t_scope, 100)));

        % 3. パラメータ交渉と設定

        % 要求されたサンプリングレートとレコード長からTDIVを計算
        T_span_target = record_length_request / S_user;
        TDIV_target = T_span_target / 10; % 10ディビジョン表示

        % 4. オシロスコープの設定変更

        % 通信形式の設定: Definite Block, 8-bit Byte, Binary
        write(t_scope, 'COMM_FORMAT DEF9,BYTE,BIN', 'ascii');

        % バイトオーダーの設定（リトルエンディアン）
        write(t_scope, 'COMM_ORDER LO', 'ascii');

        % 水平軸の設定
        write(t_scope, sprintf('TDIV %.6E', TDIV_target), 'ascii');

        % レコード長の設定（MSIZコマンドを使用）
        write(t_scope, sprintf('MSIZ %d', record_length_request), 'ascii');

        % 5. 取得実行

        % シングルモードで取得をアームし、完了を待機
        write(t_scope, 'TRMD SINGLE', 'ascii');
        write(t_scope, 'ARM', 'ascii');
        write(t_scope, '*WAI', 'ascii'); % 前のコマンド完了待ち

        % 6. 波形データの要求（チャンネル指定を含める）

        % バイナリ波形データ転送要求
        write(t_scope, sprintf('%s:WAVEFORM? DAT1', channel), 'ascii');

        % 7. IEEE 488.2 ヘッダの解析とデータ読み出し

        % a) 固定ヘッダ + バイトカウント文字列を読み出し
        header_fixed = char(read(t_scope, 11, 'char'));
        byte_count_str = char(read(t_scope, 9, 'char'));
        byte_count = str2double(byte_count_str); % 総ペイロードバイト数

        if isnan(byte_count) || (byte_count <= WAVEDESC_SIZE)
            error('LeCroy_Acquire:ReadError', ...
                  'バイナリヘッダのバイトカウント解析に失敗しました。');
        endif

        % b) 総ペイロードバイト数を一度に読み出す
        raw_data_stream = read(t_scope, byte_count, 'uint8');

        % c) ストリームをWAVEDESCと波形データに分割
        WAVEDESC = raw_data_stream(1:WAVEDESC_SIZE);
        raw_waveform_data = raw_data_stream(WAVEDESC_SIZE + 1:end);

        % --- III. バイナリデコードと校正 ---
        % WAVEDESCから校正係数を抽出
        % リトルエンディアンを想定してバイト配列を型キャスト

        % 1. データ幅の確認 (Offset 32, 4 bytes, int32)
        data_type_code = typecast(WAVEDESC(33:36), 'int32');

        % 2. 垂直軸ゲイン (Offset 116, 8 bytes, float64)
        V_Gain = typecast(WAVEDESC(117:124), 'double');

        % 3. 垂直軸オフセット (Offset 148, 8 bytes, float64)
        V_Offset = typecast(WAVEDESC(149:156), 'double');

        % 4. サンプル間隔 dt (Offset 160, 8 bytes, float64)
        dt = typecast(WAVEDESC(161:168), 'double');

        % 5. 水平軸オフセット T0 (Offset 176, 8 bytes, float64)
        T0 = typecast(WAVEDESC(177:184), 'double');

        % 6. 実際のレコード長 N_actual (Offset 244, 4 bytes, int32)
        N_actual = typecast(WAVEDESC(245:248), 'int32');

        % データの検証（整数同士で比較）
        if N_actual != int32(length(raw_waveform_data))
            error('LeCroy_Acquire:DataMismatch', ...
                  '記述子のレコード長 (%d) と受信データ長 (%d) が一致しません。', ...
                  N_actual, length(raw_waveform_data));
        endif

        % 7. 波形データ (D) の校正

        % raw_waveform_data は uint8 配列 (8ビットデータ)
        % 符号付きデータとして扱うために int8 に型キャスト
        D = typecast(raw_waveform_data, 'int8');

        % 電圧配列 (Y) の計算: Y = (D * V_Gain) + V_Offset
        Y_vector = (double(D) .* V_Gain) + V_Offset;

        % 時間配列 (X) の計算: X = (i * dt) + T0
        indices = 0:(N_actual - 1);
        X_vector = (double(indices) .* dt) + T0;

        % 8. 出力構造体の構築
        waveform.X = X_vector(:);
        waveform.Y = Y_vector(:);

        actual_params.SamplingRate = 1/dt;
        actual_params.RecordLength = N_actual;
        actual_params.TimePerDiv = dt * N_actual / 10;
        actual_params.V_Gain = V_Gain;
        actual_params.V_Offset = V_Offset;
        actual_params.T0 = T0;

    catch err
        warning('LeCroy_Acquire:FatalError', ...
                'クリティカルなエラーが発生しました: %s', err.message);
        waveform = struct(); % 空の構造体を返す
        actual_params = struct();

    end_try_catch

    % 9. クリーンアップと状態復元
    if isobject(t_scope)
        try
            % 元のタイムベース設定を復元
            if ~isempty(original_tdiv_str)
                write(t_scope, sprintf('TIM:MAIN:SCAL %s', original_tdiv_str), 'ascii');
            endif

%            % 元のタイムベース設定を復元
%            if ~isempty(original_tdiv_str)
%                write(t_scope, original_tdiv_str, 'ascii');
%            endif

            % 元のトリガモードを復元
            if ~isempty(original_trmd_str)
                write(t_scope, original_trmd_str, 'ascii');
            endif
        catch cleanup_err
            warning('LeCroy_Acquire:CleanupError', ...
                    '設定復元中にエラーが発生しました: %s', cleanup_err.message);
        end_try_catch

        % TCP接続の解放
        clear t_scope;
    endif
endfunction
