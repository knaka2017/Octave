function [time_data, voltage_data, raw_data] = get_rigol_waveform2(ip_address, port_number, channel, num_samples)
    % Load required packages
    pkg load instrument-control

    % --- 1. 接続 ---
    instr = tcpclient(ip_address, port_number);
    instr.Timeout = 60;
    instr.ByteOrder = "big-endian";
    flush(instr);

    % --- 2. 実行前の状態を記憶 (復元用) ---
    printf('\n--- Backing up instrument state ---\n');
    % 時間軸のバックアップ
    write(instr, ":TIMebase:MAIN:SCALe?\n");
    OriginalTimeScale = str2num(strtrim(readline(instr)));

    % チャンネルON/OFF状態のバックアップ (CH1-CH4)
    ch_states = zeros(1, 4);
    for n = 1:4
        write(instr, [":CHANnel", num2str(n), ":DISPlay?\n"]);
        ch_states(n) = str2num(strtrim(readline(instr)));
    end

    % --- 3. 計測用セットアップ (500M / 4GSa/s モード) ---
    write(instr, ":STOP\n"); pause(0.5);
    write(instr, ":ACQuire:MDEPth 500M\n"); pause(0.5);

    % パラメータ取得
    write(instr, "WAVeform:PREamble?\n");
    preamble_raw = readline(instr);
    preamble_vals = str2num(strrep(strtrim(preamble_raw), ',', ' '));
    XINC = preamble_vals(5);
    XORG = preamble_vals(6);
    YINC = preamble_vals(8); YORG = preamble_vals(9); YREF = preamble_vals(10);

    write(instr, ":ACQuire:MDEPth?\n");
    MDepth = str2num(strtrim(readline(instr)));

    % --- 4. ターゲットインデックス計算 ---
    TriggerPoint_In_Memory = round((0 - XORG) / XINC) + 1;
    StartIndex_Final = round(TriggerPoint_In_Memory - (num_samples / 2));
    StopIndex_Final  = StartIndex_Final + num_samples - 1;

    if StartIndex_Final < 1, StartIndex_Final = 1; StopIndex_Final = min(MDepth, num_samples); end
    if StopIndex_Final > MDepth, StopIndex_Final = MDepth; StartIndex_Final = max(1, StopIndex_Final - num_samples + 1); end
    Actual_Samples_Read = StopIndex_Final - StartIndex_Final + 1;

    % --- 5. データダウンロード (10Mチャンク) ---
    write(instr, "WAVeform:MODE RAW\n");
    write(instr, [":WAVeform:STARt ", num2str(StartIndex_Final), "\n"]);
    write(instr, [":WAVeform:STOP ", num2str(StopIndex_Final), "\n"]);

    flush(instr);
    write(instr, "WAVeform:DATA?\n");
    header1 = read(instr, 2, 'char');
    N = str2num(char(header1(2)));
    header2 = read(instr, N, 'char');

    raw_data = uint16([]);
    chunk_size = 10000000;
    samples_read = 0;
    while samples_read < Actual_Samples_Read
        current_batch = min(chunk_size, Actual_Samples_Read - samples_read);
        batch_data = read(instr, current_batch, 'uint16');
        raw_data = [raw_data; batch_data(:)];
        samples_read = length(raw_data);
        printf('Progress: %.1f%%\n', (samples_read/Actual_Samples_Read)*100);
    end
    read(instr, 1, 'char');

    % --- 6. 状態の完全復旧 ---
    printf('\n--- Restoring instrument state ---\n');

    % 時間軸の復元
    write(instr, [":TIMebase:MAIN:SCALe ", num2str(OriginalTimeScale), "\n"]);

    % チャンネルON/OFFの復元
    for n = 1:4
        write(instr, [":CHANnel", num2str(n), ":DISPlay ", num2str(ch_states(n)), "\n"]);
    end

    write(instr, ":ACQuire:MDEPth AUTO\n");
    write(instr, ":RUN\n");
    clear instr;

    % --- 7. スケーリング・プロット ---
    voltage_data = (double(raw_data) - YREF) .* YINC + YORG;
    Time_Start = XORG + (StartIndex_Final - 1) * XINC;
    time_data = (0:(length(voltage_data)-1))' .* XINC + Time_Start;

    Screen_Left = -5 * OriginalTimeScale;
    Screen_Right = 5 * OriginalTimeScale;
    indices = (time_data >= Screen_Left) & (time_data <= Screen_Right);

    figure(1);
    if any(indices), plot(time_data(indices), voltage_data(indices), 'b');
    else plot(time_data, voltage_data, 'r'); end
    xlim([Screen_Left, Screen_Right]); grid on;
    xlabel('Time [s]'); ylabel('Voltage [V]');
    line([0, 0], [min(voltage_data), max(voltage_data)], 'Color', 'r', 'LineStyle', '--');

    printf('Done.\n');
end
