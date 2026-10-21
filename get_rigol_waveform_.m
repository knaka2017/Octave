function [time_data, voltage_data] = get_rigol_waveform(ip_address, port_number, channel, num_samples)
    % Load required packages
    pkg load instrument-control

    % --- 1. 接続 ---
    printf('\n--- Connection ---\n');
    instr = tcpclient(ip_address, port_number);
    instr.Timeout = 60; % 転送時間を考慮して長めに設定
    instr.ByteOrder = "big-endian";
    flush(instr);

    % --- 2. 状態取得と最高レート設定 ---
    write(instr, ":TIMebase:MAIN:SCALe?\n");
    OriginalTimeScale = str2num(strtrim(readline(instr)));

    write(instr, ":STOP\n"); pause(0.5);
    write(instr, ":ACQuire:MDEPth 500M\n"); pause(0.5);

    % 3. パラメータ取得
    write(instr, "WAVeform:PREamble?\n");
    preamble_raw = readline(instr);
    preamble_vals = str2num(strrep(strtrim(preamble_raw), ',', ' '));

    XINC = preamble_vals(5);
    XORG = preamble_vals(6);
    YINC = preamble_vals(8); YORG = preamble_vals(9); YREF = preamble_vals(10);

    write(instr, ":ACQuire:MDEPth?\n");
    MDepth = str2num(strtrim(readline(instr)));

    % --- 4. インデックス計算 (0sを中心に配置) ---
    TriggerPoint_In_Memory = round((0 - XORG) / XINC) + 1;
    StartIndex_Final = round(TriggerPoint_In_Memory - (num_samples / 2));
    StopIndex_Final  = StartIndex_Final + num_samples - 1;

    % 境界ガード
    if StartIndex_Final < 1, StartIndex_Final = 1; StopIndex_Final = min(MDepth, num_samples); end
    if StopIndex_Final > MDepth, StopIndex_Final = MDepth; StartIndex_Final = max(1, StopIndex_Final - num_samples + 1); end

    Actual_Samples_Read = StopIndex_Final - StartIndex_Final + 1;

    % --- 5. 分割読み込みによるデータダウンロード ---
    write(instr, "WAVeform:MODE RAW\n");
    write(instr, [":WAVeform:STARt ", num2str(StartIndex_Final), "\n"]);
    write(instr, [":WAVeform:STOP ", num2str(StopIndex_Final), "\n"]);

    flush(instr);
    write(instr, "WAVeform:DATA?\n");

    % ヘッダ解析 (#N DDDD...)
    header1 = read(instr, 2, 'char');
    N = str2num(char(header1(2)));
    header2 = read(instr, N, 'char');
    total_expected_bytes = str2num(char(header2));

    printf('Download Start: Total %d bytes expected (%d samples).\n', total_expected_bytes, Actual_Samples_Read);

    % ★分割読み込みロジック (10Mチャンク設定)
    raw_data = uint16([]);
    chunk_size_samples = 10000000; % 10M点に設定
    samples_read = 0;

    while samples_read < Actual_Samples_Read
        remaining = Actual_Samples_Read - samples_read;
        current_batch = min(chunk_size_samples, remaining);

        % チャンク読み込み
        batch_data = read(instr, current_batch, 'uint16');

        if isempty(batch_data)
            error('Download error: Received empty data chunk. Check network or timeout.');
        end

        % 結合
        raw_data = [raw_data; batch_data(:)];
        samples_read = length(raw_data);

        printf('Progress: %d / %d samples read (%.1f%%).\n', ...
               samples_read, Actual_Samples_Read, (samples_read/Actual_Samples_Read)*100);
    end

    read(instr, 1, 'char'); % 終端LF
    printf('Download complete.\n');

    % --- 6. クリーンアップと復元 ---
    write(instr, [":TIMebase:MAIN:SCALe ", num2str(OriginalTimeScale), "\n"]);
    write(instr, ":ACQuire:MDEPth AUTO\n");
    write(instr, ":RUN\n");
    clear instr;

    % --- 7. スケーリングと時間軸生成 ---
    voltage_data = (double(raw_data) - YREF) .* YINC + YORG;
    Time_Start = XORG + (StartIndex_Final - 1) * XINC;
    time_data = (0:(length(voltage_data)-1))' .* XINC + Time_Start;

    % --- 8. プロット ---
    Screen_Left = -5 * OriginalTimeScale;
    Screen_Right = 5 * OriginalTimeScale;
    indices_in_screen = (time_data >= Screen_Left) & (time_data <= Screen_Right);

    figure(1);
    if any(indices_in_screen)
        plot(time_data(indices_in_screen), voltage_data(indices_in_screen), 'b');
    else
        plot(time_data, voltage_data, 'r');
    end

    xlim([Screen_Left, Screen_Right]);
    line([0, 0], [min(voltage_data), max(voltage_data)], 'Color', 'r', 'LineStyle', '--');
    grid on;
    xlabel('Time [s]'); ylabel('Voltage [V]');
    title('Waveform (10M Chunk Download)');

    printf('Process finished.\n');
end
