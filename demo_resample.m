function demo_resample()
% =========================================================
% デモ・動作確認用スクリプト
% このファイルを直接実行すると以下のデモが走ります
% =========================================================

  pkg load signal;

  fprintf('=== サンプリングレート変換デモ ===\n\n');

  % 共通パラメータ
  Fsi   = 2e9;   % 入力サンプリング周波数 [Hz]
  Nin   = 1638;    % 入力信号長
  f_sig = 10e6;    % テスト信号周波数 [Hz]（A4音）

  % 入力信号生成（正弦波）
  t_in = (0:Nin-1) / Fsi;
  x    = sin(2 * pi * f_sig * t_in);

  % ---- ケース1: アップサンプリング、Nout=時間スパン保持 ----
  Fso1  = 2.5e9;
  Nout1 = 2048;
  y1    = resample_signal(x, Fsi, Fso1, Nout1);
  fprintf('[ケース1] Fsi=%d Hz → Fso=%d Hz, Nout=%d (アップサンプリング・時間スパン保持)\n', Fsi, Fso1, Nout1);
  fprintf('  入力長 Nin=%d, 出力長=%d\n\n', Nin, length(y1));

  % ---- ケース2: ダウンサンプリング、Nout=時間スパン保持 ----
  Fso2  = 1e9;
  Nout2 = 819;
  y2    = resample_signal(x, Fsi, Fso2, Nout2);
  fprintf('[ケース2] Fsi=%d Hz → Fso=%d Hz, Nout=%d (ダウンサンプリング・時間スパン保持)\n', Fsi, Fso2, Nout2);
  fprintf('  入力長 Nin=%d, 出力長=%d\n\n', Nin, length(y2));

  % ---- ケース3: Nout を任意に指定（ゼロパディング発生）----
  Fso3  = 2.5e9;
  Nout3 = 4096;
  y3    = resample_signal(x, Fsi, Fso3, Nout3);
  fprintf('[ケース3] Fsi=%d Hz → Fso=%d Hz, Nout=%d (ゼロパディングあり)\n', Fsi, Fso3, Nout3);
  fprintf('  入力長 Nin=%d, 出力長=%d\n\n', Nin, length(y3));

  % ---- ケース4: Nout を任意に指定（切り捨て発生）----
  Fso4  = 2.5e9;
  Nout4 = 1024;
  y4    = resample_signal(x, Fsi, Fso4, Nout4);
  fprintf('[ケース4] Fsi=%d Hz → Fso=%d Hz, Nout=%d (切り捨てあり)\n', Fsi, Fso4, Nout4);
  fprintf('  入力長 Nin=%d, 出力長=%d\n\n', Nin, length(y4));

  % ---- プロット ----
  figure('Name', 'サンプリングレート変換デモ', 'NumberTitle', 'off');

  t1_out = (0:length(y1)-1) / Fso1;
  t2_out = (0:length(y2)-1) / Fso2;
  t3_out = (0:length(y3)-1) / Fso3;

  subplot(4,1,1);
  plot(t_in  * 1000, x,  'b-o', 'MarkerSize', 3);
  title(sprintf('入力信号 (Fsi=%d Hz, Nin=%d)', Fsi, Nin));
  xlabel('時間 [ms]'); ylabel('振幅'); grid on;

  subplot(4,1,2);
  plot(t1_out * 1000, y1, 'r-o', 'MarkerSize', 3);
  title(sprintf('ケース1: アップサンプリング (Fso=%d Hz, Nout=%d)', Fso1, length(y1)));
  xlabel('時間 [ms]'); ylabel('振幅'); grid on;

  subplot(4,1,3);
  plot(t2_out * 1000, y2, 'm-o', 'MarkerSize', 3);
  title(sprintf('ケース2: ダウンサンプリング (Fso=%d Hz, Nout=%d)', Fso2, length(y2)));
  xlabel('時間 [ms]'); ylabel('振幅'); grid on;

  subplot(4,1,4);
  plot(t3_out * 1000, y3, 'k-o', 'MarkerSize', 3);
  title(sprintf('ケース3: ダウンサンプリング+ゼロパディング (Fso=%d Hz, Nout=%d)', Fso3, length(y3)));
  xlabel('時間 [ms]'); ylabel('振幅'); grid on;

  fprintf('デモ完了。グラフを確認してください。\n');

end
