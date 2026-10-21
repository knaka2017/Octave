function [time_data, voltage_data, raw_data] = get_rigol_waveform(ip_address, port_number, channel, num_samples)
    pkg load instrument-control

    % --- 1. 接続設定 ---
    printf('\n--- Connecting to %s ---\n', ip_address);
    instr = tcpclient(ip_address, port_number);
    instr.Timeout = 60;
    instr.ByteOrder = "big-endian";
    flush(instr);

    % --- 2. 状態のバックアップ ---
    write(instr, ":TIMebase:MAIN:SCALe?\n");
    OriginalTimeScale = str2num(strtrim(readline(instr)));

    ch_states = zeros(1, 4);
    for n = 1:4
        write(instr, [":CHANnel", num2str(n), ":DISPlay?\n"]);
        ch_states(n) = str2num(strtrim(readline(instr)));
    end

    % --- 3. 計測セットアップ (500M / 4GSa/s 強制) ---
    write(instr, ":STOP\n"); pause(0.2);
    write(instr, ":ACQuire:MDEPth 500M\n"); pause(0.2);
    write(instr, ":RUN\n");

    % 500Mポイントを埋めるための待機 (最低0.125s必要)
    printf('Filling 500Mpts memory (Wait 0.5s)...\n');
    pause(0.5);
    write(instr, ":STOP\n"); pause(0.2);

    % --- 4. パラメータの取得と解析 ---
    write(instr, "WAVeform:PREamble?\n");
    p = str2num(strrep(strtrim(readline(instr)), ',', ' '));
    % p(5):XINC, p(6):XORG, p(8):YINC, p(9):YORG, p(10):YREF
    XINC = p(5); XORG = p(6); YINC = p(8); YORG = p(9); YREF = p(10);

    printf('\n--- Debug Parameters ---\n');
    printf('  XINC: %e, XORG: %e\n', XINC, XORG);
    printf('  YINC: %e, YORG: %e, YREF: %e\n', YINC, YORG, YREF);

    write(instr, ":ACQuire:MDEPth?\n");
    MDepth = str2num(strtrim(readline(instr)));

    % トリガー点(0s)を基準にした抽出インデックス計算
    TriggerIdx = round((0 - XORG) / XINC) + 1;
    StartIdx = round(TriggerIdx - (num_samples / 2));
    StopIdx  = StartIdx + num_samples - 1;

    % 境界ガード
    if StartIdx < 1, StartIdx = 1; StopIdx = min(MDepth, num_samples); end
    if StopIdx > MDepth, StopIdx = MDepth; StartIdx = max(1, StopIdx - num_samples + 1); end
    ActualPoints = StopIdx - StartIdx + 1;

    % --- 5. 分割ダウンロード (10Mチャンク) ---
    write(instr, [":WAVeform:SOURce CHANnel", num2str(channel), "\n"]);
    write(instr, ":WAVeform:FORMat WORD\n");
    pause(0.1);

    write(instr, "WAVeform:MODE RAW\n");
    write(instr, [":WAVeform:STARt ", num2str(StartIdx), "\n"]);
    write(instr, [":WAVeform:STOP ", num2str(StopIdx), "\n"]);

    flush(instr);
    write(instr, "WAVeform:DATA?\n");
    header1 = read(instr, 2, 'char');
    N = str2num(char(header1(2)));
    header2 = read(instr, N, 'char');

    printf('Downloading %d samples...\n', ActualPoints);
    raw_data = uint16([]);
    chunk_size = 10000000; % 10Mチャンク
    while length(raw_data) < ActualPoints
        batch = min(chunk_size, ActualPoints - length(raw_data));
        batch_data = read(instr, batch, 'uint16');
        raw_data = [raw_data; batch_data(:)];
        printf('  Progress: %.1f%%\n', (length(raw_data)/ActualPoints)*100);
    end
    read(instr, 1, 'char'); % 終端LF

    % --- 6. 状態の復元 ---
    write(instr, [":TIMebase:MAIN:SCALe ", num2str(OriginalTimeScale), "\n"]);
    for n = 1:4
        write(instr, [":CHANnel", num2str(n), ":DISPlay ", num2str(ch_states(n)), "\n"]);
    end
    write(instr, ":ACQuire:MDEPth AUTO\n");
    write(instr, ":RUN\n");
    clear instr;

    % --- 7. スケーリングとオフセット補正 ---
    data_double = double(raw_data);

    % 基本式: (raw - YREF) * YINC + YORG
    voltage_data = (data_double - YREF) .* YINC + YORG;

    % ★オフセット自動補正ロジック
    % もし電圧が100Vを超えるような異常値（-500V等）の場合、
    % YORGをデジタル値として扱い (raw - YREF - YORG) * YINC を試みる
    if abs(mean(voltage_data)) > 100
        printf('Applying offset correction (Large offset detected: %.1fV)\n', mean(voltage_data));
        voltage_data = (data_double - YREF - YORG) .* YINC;
    end

    % 時間軸の計算
    T_Start = XORG + (StartIdx - 1) * XINC;
    time_data = (0:(length(voltage_data)-1))' .* XINC + T_Start;

    % --- 8. 結果表示とプロット ---
    printf('\n--- Final Statistics ---\n');
    printf('  Vmax: %.4f V, Vmin: %.4f V, Vpp: %.4f V\n', max(voltage_data), min(voltage_data), max(voltage_data)-min(voltage_data));

    Screen_L = -5 * OriginalTimeScale;
    Screen_R = 5 * OriginalTimeScale;

    figure(1);
    plot(time_data, voltage_data, 'b');
    xlim([Screen_L, Screen_R]);
    grid on;
    line([0, 0], [min(voltage_data), max(voltage_data)], 'Color', 'r', 'LineStyle', '--');
    xlabel('Time [s]'); ylabel('Voltage [V]');
    title(sprintf('CH%d Waveform (Download: %d pts)', channel, length(raw_data)));

    printf('Done.\n');
end
