% VLC OFDM ver.14
% 実験用(wExp="Y", wTxdata="Y")：
% １．送信信号ファイルの名前設定　fntxdata
% ２．実行（送信信号ファイル保存のみ．グラフ意味なし）
% ３．"tx2comb_tx 送信信号ファイル名"を実行してDAC用ファイル生成
% ４．DAC用FPGAに書き込み
% ５．実験して，受信信号ファイルを取得
% ６．受信信号ファイルの名前設定 fnadc
% ７．StartSymOfBerCalc:1, Stop...:BLEN-1に設定
% ８．実行．BER,EVM
% ９－１．Figure 3のFEQをチェックし，絶対値・偏角とも
% 　　なるべく平坦かつ長期間のところを探し，上記７を変更して再実行
% ９－２．Lcは移動平均のサンプル数で大きいほど（<=OvRate）雑音が小さくなるが，
% 　　　　　大きすぎると，隣のサンプルの影響も大きく受けるので注意（特に信号歪みが大きい場合は顕著）
% 　　　　　Lcを変更してみるとBER,EVMが大きく改善する場合有り
% ex1:Figure 5の相関波形に鋭いピークが表れていない場合，そもそも波形が
%    取得できていない可能性あり
%     corrPos でピーク位置の選択可能
% ex2:Figure 9の表示に時間がかかるため，最初は非表示(wPLOTeachSub="n")
%     そこそこ良いBER,EVMになったときに表示設定(wPLOTeachSub="y")
% ex3:chkSubCで，使用するサブキャリアを指定可能．IM/DDだと通常1:((N/2)-1)
%     例：2:10だと，サブキャリア#2～#10の9本のみ使用（残りは振幅ゼロ）
%     図4, 7, 8が目安の基準
%
%% 調整が必要なパラメータ：
% corrPos: 時間波形の先頭位置検出用．相関値のピーク位置(〇)の番号．ピーク間の間隔が取得し計算するブロック
% chkSubC: BER,EVM計算の対象とするサブキャリアの番号のベクトル．SNRの悪いサブキャリアを除くことでBER,EVM向上
% StartSymOfBerCalc, StopSymOfBerCalc: 計算に利用するブロック番号（開始，終了）．狭いほうがBER,EVMは小さくなるが，計算点数が少なくなってしまう
% Lc: 1送信サンプル分のフィルタ(1～OvRate)
% OvRate: 1送信サンプルの幅（DACサンプル間隔を1として．例:8だと同値を８回繰り返し送信．矩形波の幅が広がるイメージ）
%% OFDM主要パラメータ
% N, P, Ms
%
%% ADC動作不良？
% コマンドラインでの指定で，一度，1ch, 2ch, 4ch，(A,B,C,D)ポートの
% 設定を他の設定で実行してみると，ADC(FPGA)の内部情報がリセットされて
% 正常に戻ることもある



clear;
ofdmVersion = 14;

tic;
pkg load signal;
pkg load communications;


wExp ="y";		% 実験(Y) or シミュレーション(other)
wConsteNonly = "n";
wReadDrx = "n";
wTxdata = "y";  % 出力ファイル作成
fntxdata = "txdata251217-01_OFDM_64_QPSK_625MSps.mat";		% txdata save file name
fntxdata = "test01.mat";
switch wExp
	case {"Yes" "yes" "YES" "y" "Y"}
%%    fnadc = uigetfile("adc*.mat", "ADC file");
##fid = fopen('new_data.txt');
##fnadc = fgets(fid);
##fnadc(length(fnadc)-1)=0;
##fclose(fid);
  otherwise
    fnadc = 0;
end;
 fnadc = 0;
if (fnadc == 0)
  fnadc	= "emuCH02.mat";	% read adc data file name
end;
wPLOT="y";		% コンステ，波形表示 (表示まで数秒必要)
wPLOTeachSub="n";   % サブキャリアごとのコンステ表示（数十秒かかるので注意）
kkg = 1;     % kkg個のサブキャリアごとに平均（サブキャリア毎コンステグラフ用）
corrPos = 1;  % 相関ピーク(50%以上の振幅)の選択（相関波形の左から順に）

N = 64;						% キャリア数
P = 4;						% CP長
chkSubC = 2:(floor(N/2)-2);    % 使用するサブキャリア番号のベクトル（IM/DDでは最大(N/2)-1）

%% load("-mat", "zbzp.mat");

OvRate = 16;     % 2.5GHz(DAC)に対しての数値
BLEN = floor(408*(4/OvRate)*(80/(P+N)));				% 送信ブロック数（誤り率を正確に->大，波形を見たい->小）
%BLEN = floor(((4096/3)/OvRate)/(P+N));

StartSymOfBerCalc = 10;		% BER,EVM計算開始シンボル
StopSymOfBerCalc = 90;    % BER,EVM計算終了シンボル
StartSymOfBerCalc = 1;  StopSymOfBerCalc = BLEN;
%%LFEQ = 80;
LFEQ = StopSymOfBerCalc - StartSymOfBerCalc;	% FEQ 平均ブロック数

L = 1;						% FIRフィルタのタップ数
Lc= OvRate;						% length of chip matched filter (fsr[Hz])

Nf = 1*N;					% FFT次数
Ms = 4;					% 変調方式 (4:QPSK, 他: QAM)
Ms = ones(N/2-1,1)*Ms;		% 	サブキャリア毎に変調方式指定(IM/DDなので長さ(N/2)-1)
%Ms(chkSubC) = 2.^zb(1:length(chkSubC));
%chkSubC256 = 2:3;
%chkSubC064 = 4:7;
%chkSubC016 = 8:12;
%chkSubC004 = 13:(floor(N/2)-1);
##Ms(chkSubC256) = 256;
##Ms(chkSubC064) = 64;
##Ms(chkSubC016)= 16;
##Ms(chkSubC004)= 4;
%Ms = 4.^randint(N/2-1,1, [1 4], 1234);

fs = 2.5e9;					% DACサンプリング周波数  オーバーサンプリング比:fs/fdac
fsr = 2.5e9;			  % ADC(受信)サンプリング周波数
fdac = fs/OvRate;				% 信号サンプル速度 (<= fs) （オーバーサンプル無し）
dt = 1/fs;

Ns = 0;						% ブロックスタート開始時刻-1 (0で最初から)
%ddsmp = floor(fs/fdac)-1;	% ダウンサンプルタイミングディレイ
Nq = 14;					% DAC量子化レベル 2^Nq

h=1;
hsb=[5.6248e-05, 0, -1.1250e-04, 0, 5.6248e-05];
hsa=[1.0000, -3.9786, 5.9361, -3.9364, 0.9789];
hsb=1;
hsa=1;
w=1;
ws=1;
%% sinc pulse
%t = 0:dt:((500-1)*dt);
%ws = sinc(2*(fdac/2)*(t-(2/fdac)));
%% raised cosine pulse
%t = 0:dt:((32-1)*dt); t = t - max(t)/2; TT = 1/fdac; Roff = 0.4; % Roff = 0.4;
%ws = (sin(pi*t/TT)./(pi*t/TT)).*((cos(pi*Roff*t/TT))./(1-4*((Roff*t/TT).^2)));
SEEDs = 12344;

%% FEQ係数，送信側 Dtx，受信側 Drx
Dtx = zeros(1, (N/2)-1);
Dtx(chkSubC) = 1;	% 対角要素1 （FEQ無しと等価）
%Dtx(chkSubC256) = 1.1   + 0.1   * [1:length(chkSubC256)];
%Dtx(chkSubC064) = 0.5   + 0.1   * [1:length(chkSubC064)];
%Dtx(chkSubC016) = 0.25  + 0.1  * [1:length(chkSubC016)];
%Dtx(chkSubC004) = 0.125 + 0.1 * [1:length(chkSubC004)];
%Dtx(chkSubC) = 10.^(zp/20);
%Dtx(400:500) = 4;
%%Dtx = ones(1,(N/2)-1);
%%load("-mat", "Dtx01.mat");
Dtx = diag([0, Dtx, 0, conj(fliplr(Dtx))]);
              % 中央の"1"は伝送帯域幅の上限の明示用．通常は"0"

%Dtx = diag(linspace(1,2,N)');

%%%% 通信ブロック

%% generate data symbol
%sb = randint(N, BLEN, [0, Ms-1], SEEDs);
%sb = randint(N/2-1, BLEN, [0, Ms-1], SEEDs);	% DMT
sb = randintM(N/2-1, BLEN, Ms, SEEDs);

%% pilot signal
%sb(:,1) = ones(N/2-1,1)*18; 	% pilot signal
%sb(:,1) = randsrc(N/2-1,1,[0 7 56 63], SEEDs); 	% pilot signal
%sb(:,1) = randsrc(N/2-1,1,[35 36 27 28], SEEDs); 	% pilot signal

s = qammodM(sb, Ms);
%s = ones(size(s)) + j*ones(size(s));
%s = abs(s);
s0 = s;

disp(size(zeros(1,BLEN)));
disp(size(s));
disp(size(ones(1,BLEN)));
disp(size(conj(flipud(s))));
s = [zeros(1,BLEN); s; ones(1,BLEN); conj(flipud(s))]; %DMT

s = Dtx*s;								% pre-FEQ

x = ifft(s, Nf);                        % IFFT
% x = repmat([0:(Nf-1)]', 1,BLEN)/Nf;		% ramp wave 線形性チェック用

u = [x((Nf-P+1):Nf, :); x];				% ins CP
uu = reshape(u, (P+Nf)*BLEN, 1);		% P/S
uu0 = uu;

%uuK = 20;
%uu = atan((real(uu)-10)/uuK)*uuK+10 + i*(atan((imag(uu)-10)/uuK)*uuK+10);	% nonlinear

Upow = mean(abs(uu).^2);				%

%% TEQ
%uu = filter(w, 1, uu);					% Pre-TEQ

%uu = filter(1, h, uu);					% zero-forcing IIR filter (transmitter)


%% DAC前処理 (量子化）2^Nq
[uui, uusi, uusi_imp, LSB] = DACout03(uu, Nq, fs, fdac, 1);
uui = real(uui);	% DMT: fs = fdac
uusi= real(uusi);	% DMT
uusi_imp = real(uusi_imp);

% ADCサンプル(fsr[Hz])の場合 (受信側リファレンス，確認用)
[uuir, uusir, uusir_imp, LSBr] = DACout03(uu, Nq, fsr, fdac, 1);
uuir = real(uuir);	% DMT: fs = fdac
uusir= real(uusir);	% DMT
uusir_imp = real(uusir_imp);


% 正規化
UUpow = mean(abs(uui-mean(uui)).^2);
UUpow = max(abs(uui));
uu = uui/UUpow;
uus= uusi/UUpow;

uur = uuir/UUpow;
uusr= uusir/UUpow;

uus = filter(ws, 1, uus);

%% channel

%rr = filter(h, 1, uu);                  % channel
rrs = filter(hsb, hsa, uus);

%rrs = awgn(rrs, 300, 'measured', 'dB');			% add Noise


%% ADC

% pseudo drift start position ## for sim.
rrs = filter([zeros(1,401),1], 1, [rrs,rrs]);
rrs(1:401) = randn(1,401);

% upsample fs -> fsr ## for sim.
rrs0s = rrs;
rrs = resample(rrs, fsr, fs);

switch wExp
	case {"Yes" "yes" "YES" "y" "Y"}

	%% save tx_data

%%	temp_uus = [uus uus uus uus uus uus uus uus uus uus uus uus uus]; %% x12 repeat
  temp_uus = repmat(uus, 1, ceil((2^20)/length(uus)));
	temp_uus = temp_uus - min(temp_uus);
	temp_uus = temp_uus / max(abs(temp_uus));
	temp_uus = temp_uus * (2^Nq-1);
	tx_data = floor(temp_uus);
	tx_data = tx_data(1:2^20);
  switch wTxdata
	  case {"Yes" "yes" "YES" "y" "Y"}
      save("-mat", fntxdata, "tx_data");
    otherwise
%%      Hwarndlg = warndlg("Do not generate TX file");
  end;

	%% load adc data

	if exist(fnadc) == 2
		load("-mat", fnadc);
		if (exist("yyyA") == 1)
			disp("*** 2.5GHz A ch.");
			rrs = yyyA;
		elseif (exist("yyyB") == 1)
			disp("*** 2.5GHz B ch.");
			rrs = yyyB;
		else
			disp("*** 5GHz 1ch");
			rrs = yyy;
		end
##cd(fnadc);
##y0 = dlmread('adc0.txt');
##y1 = dlmread('adc1.txt');
##y2 = dlmread('adc2.txt');
##y3 = dlmread('adc3.txt');
##rrs = reshape([y2'; y3'], 1,4096);
##rrs0 = rrs;
##rrs = (rrs - 127)/128;
##cd ..;

		rrs = -rrs;		% invert (for inverting AMP.)

	else
		disp(sprintf("file not found:%s", fnadc));
	end;
	otherwise
  rrs = awgn(rrs, 300, 'measured', 'dB');			% add Noise
end;

rrs0 = rrs;		% 受信信号 fsr [sps]
rrs = filter(ones(Lc,1)/Lc,1,rrs);	% chip-matched filter
rrs00=rrs;

% detect start position
[R, lag] = xcorr(rrs, uusr, 'none');	% 送受信信号の相関
[pks, idx] = findpeaks(abs(R), "MinPeakHeight", 0.3*max(abs(R)), "MinPeakDistance", length(uu)*OvRate/2);
                                              % 最大相関値とその位置
ddsmp = abs(lag(idx(corrPos)+1));				% 信号ブロック開始位置
switch wConsteNonly
	case {"Yes" "yes" "YES" "y" "Y"}
    figure(5);  %% 相関波形（鋭いピークが3本以上出ていないとNG）
  clf;
  plot(lag, R, "-", lag(idx), pks, "o");
  title(sprintf("Corr.(tx and rx signals) (with %d peaks)", length(idx)));
  xlabel("lag [samples]"); ylabel("Correlation");
otherwise
end;

rr = downsample(rrs(ddsmp:end), floor(fsr/fdac)); % ADC
if (length(rr) < length(uu))
  disp("Error: \"corrPos\" requires to be changed.");
  disp(sprintf("current: corrPos = %d", corrPos));
##  close(Hwarndlg);
end;
%rr = decimate(rrs, floor(fs/fdac));
rr = rr(1:length(uu));

%Rpow = mean(abs(rr-mean(rr)).^2);

rx = rr;
%nr = (randn(size(rr)) + j*randn(size(rr))) *sqrt(Rpow*sv/2);
%rx = rr + nr;
%         rx = awgn(rx, 40, 'measured', 'dB');			% add Noise


%% receiver side

rx0 = rx;
rx = [rx((Ns+1):length(rx)); zeros(Ns,1)];
rx1 = rx;

% rx = filter(1, h, rx);                % zero-forcing IIR filter (receiver)

rx = filter(w, 1, rx);                  % Post TEQ

r = reshape(rx, P+Nf, BLEN);             % S/P

xh = r((P+1):(P+Nf), :);                 % del CP

sh = fft(xh, Nf);                        % FFT

%sh = diag(exp(-j*(2*pi/Nf)*(0:(Nf-1))*Ns)) * sh;

sh = sh(2:(N/2),:);

sh0= sh;

%% calc. FEQ coeff.
switch wReadDrx
	case {"Yes" "yes" "YES" "y" "Y"}
  load("-mat", "Drx.mat");
otherwise
%Drx = diag(s0(:,1) ./ sh0(:,1));	% use only first pilot symble
xDrx00 = (s0 ./ sh0);
xDrx0 = xDrx00(:, StartSymOfBerCalc:StopSymOfBerCalc);
iy = 1:(StopSymOfBerCalc-StartSymOfBerCalc+1);
xDrx = zeros(size(xDrx0));
for ix = 1:((N/2)-1)
	pa = polyfit(iy, abs(xDrx0(ix,:)),4);
	ps = polyfit(iy, unwrap(angle(xDrx0(ix,:))),5);
	xDrx(ix,:) = polyval(pa, iy).*exp(j*polyval(ps, iy));
end
Drx = diag((mean(real(xDrx(:,1:LFEQ))')+j*mean(imag(xDrx(:,1:LFEQ))')));
%Drx = diag(mean(xDrx)');
%Drx = diag(ones(length(Drx),1))*6.7657;
save("-mat", "Drx.mat", "Drx", "xDrx");
end;
sh = Drx*sh;                    % FEQ
sh(:,StartSymOfBerCalc:StopSymOfBerCalc) = xDrx.*(inv(Drx)*sh(:,StartSymOfBerCalc:StopSymOfBerCalc));

sd = qamdemodM(sh, Ms);			% QAM復調

pause(0.1);


%%% calc.

%Nd = (N/2-2)*(BLEN-1)*log2(Ms);             % BER計算用総ビット数
Nd = (StopSymOfBerCalc-StartSymOfBerCalc+1)*sum(log2(Ms(chkSubC)));
chk01 = biterr(sb(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc), sd(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc))/Nd;  % データ(bit)チェック
PAPR_0 = max(abs(uu0).^2)/mean(abs(uu0).^2);	% w/o filter
PAPR_1 = max(max(abs(uu-mean(mean(uu))).^2)) /mean(mean(abs(uu-mean(mean(uu))).^2));	% w/ filter
%PAPR_r = max(abs(rr).^2) /mean(abs(rr).^2);     % received
%PAPR_rd = max(abs(rr-mean(rr)).^2) /mean(abs(rr-mean(rr)).^2);     % received
EVMrms = sqrt( mean( mean( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) - sh(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 ) ) / mean(mean( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 )) ); % 平均電力で正規化
EVMrmssb = sqrt( mean( ( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) - sh(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 )' ) ./ mean( ( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 )') ); % 平均電力で正規化
EVMrmspeak = sqrt( mean( mean( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) - sh(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 ) ) / max(max( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 )) );  % 最大QAMシンボルのパワーで正規化
EVMpeak = sqrt( max( max( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) - sh(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 ) ) / max(max( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 )) );
Prr = mean(rr);


disp(sprintf("rand seed s=%d\n", SEEDs));
disp(sprintf("N=%d, Nq=%d", N, Nq));
%disp(Ms');

disp(sprintf('Data Rate = %f [Mbps]', (sum(log2(Ms(chkSubC)))/((P+N)*dt*OvRate))*1e-6));
disp(sprintf('OFDM sample speed = %f [Msps]', fdac*1e-6));
disp(sprintf('BER = %g, NumOfBits. = %ld, ', chk01, Nd));
%disp(sprintf('SSNR = %g[dB](%g), \n0: w/o filter, 1: w/ filter', 10*log10(SSNR), SSNR));
disp(sprintf(' PAPR_0 = %g, PAPR_1 = %g,',   PAPR_0, PAPR_1));
%disp(sprintf(' PeakPow0 = %g, PeakPow1 = %g,', max(abs(uu0).^2), max(abs(uu).^2)));
%disp(sprintf(' Avg.Pow0 = %g, Avg.Pow1 = %g', mean(abs(uu0).^2), mean(abs(uu).^2)));
%disp(sprintf('SSNR at ws = %g[dB](%g)', 10*log10(ws'*B*ws / (ws'*A*ws)), ws'*B*ws / (ws'*A*ws)));
%disp(sprintf('PAPR_r = %g (%g [dB]), PAPR_r(w/o DC) = %g (%g[dB])', PAPR_r, 10*log10(PAPR_r), PAPR_rd, 10*log10(PAPR_rd)));
disp(sprintf('EVMrms = %g[%%] (SNR = %g[dB]),  EVMrmspeak = %g[%%],  EVMpeak = %g[%%]', EVMrms*100, -20*log10(EVMrms), EVMrmspeak*100, EVMpeak*100));

disp(sprintf('Calc. only: %f[sec]', toc));


%%% プロット

switch wPLOT
	case {"Yes" "yes" "YES" "y" "Y"}

figure(1);  %% コンスタレーションマップ：元データ信号
%%	clf;
  ps0 = s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc);  ps0 = reshape(ps0, prod(size(ps0)), 1);
	fig1x = real(ps0);  fig1y = imag(ps0);
  H=plot(fig1x, fig1y, '.b');
  title('Original signal');
	xlabel('In-Phase'), ylabel('Quadrature'), ...
	xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ...
	set(H, 'MarkerSize', 8);
	set(gca, 'XTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), set(gca, 'YTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), ...
	axis square;
figure(2);  %% コンスタレーションマップ：復調結果
%%	clf;
  psh = sh(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc);  psh = reshape(psh, prod(size(psh)), 1);
	fig2x = real(psh);  fig2y = imag(psh);
  H=plot(fig2x, fig2y, '.r');
  title('demod. signal');
	xlabel('In-Phase'), ylabel('Quadrature'), ...
	xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ...
	set(H, 'MarkerSize', 8);
	set(gca, 'XTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), set(gca, 'YTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), ...
	axis square;
switch wConsteNonly
	case {"Yes" "yes" "YES" "y" "Y"}
  figure(3);  %% 送信信号・受信信号（一部），FEQの絶対値・偏角の平均（全サブキャリア）の時間変化
  clf;
%%  subplot(2,2,1), plot(1:max(size(uu0(((N+P)*(BLEN-4)):((N+P)*(BLEN-1))))), real(uu0(((N+P)*(BLEN-4)):((N+P)*(BLEN-1))))), xlabel('k'), ylabel('u(k)');
%%  subplot(2,2,2), plot(1:max(size(rr(((N+P)*(BLEN-4)):((N+P)*(BLEN-1))))), real(rr(((N+P)*(BLEN-4)):((N+P)*(BLEN-1))))), xlabel('k'), ylabel('r(k)');
  ixDrx = StartSymOfBerCalc:StopSymOfBerCalc;
  subplot(2,2,3), plot(ixDrx,abs(mean(xDrx))); ylabel('abs of FEQ'); xlabel('index of blocks');
  subplot(2,2,4), plot(ixDrx,unwrap(angle(mean(xDrx)))*180/pi); xlabel('index of blocks'); ylabel('angle of FEQ [deg.]');
figure(4);
  clf;
  fig4x = chkSubC;  fig4y = mean(abs(xDrx(chkSubC,:).'));
  plot(fig4x, fig4y/min(fig4y));
  title("normalized avg. of abs(FEQ)");
  xlabel("# of Sub-Carrier");
  ylabel("norm. avg. of abs(FEQ coeff.)");
%%figure(5);  %% 相関波形（鋭いピークが3本以上出ていないとNG）
%%  clf;
%%  plot(R);
%%  title('Corr. tx and rx signals. (Require 3 or more peaks)');
figure(6);
  clf
  Lrrs0 = length(rrs0);
  RRS0 = fft(rrs0.*hanning(Lrrs0)');
  tempRRS0 = 10*log10((1000/50)*(abs(RRS0).^2)/Lrrs0);
  fig6x = (0:(Lrrs0-1))*(fs/Lrrs0)*1e-6;  fig6y = tempRRS0-max(tempRRS0);
  plot(fig6x, fig6y);
  xlabel("Freq. [MHz]");  ylabel("Normalized Power [dB]");
  xlim([0 500]);  ylim([-60 0]);
%%  trrs0 = 0:dt:((Lrrs0-1)*dt);
%%  [frrs0, Yrrs0] = fftspe2(trrs0, rrs0, 19, 1,1);
%%  plot(frrs0*1e-6, Yrrs0);
%%  xlabel("Freq. [MHz]");  ylabel("Power [dBM]");
%%  xlim([0 500]);;
  title("Power Spectrum of Received Signal (w/ hanning)");
figure(7);
  clf
  fig7x = chkSubC;
  fig7y = sqrt( ( mean( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) - sh(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) )'.^2 ) ) / mean(mean( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 )) )*100;
  stem(fig7x, fig7y);
  xlabel("# of sub carriers");  ylabel("EVM [%]");
  title("EVMs of each sub carriers");
figure(8);
  clf
  fig8x = chkSubC;
  fig8y = -20*log10(fig7y/100);
  plot(fig8x, fig8y, "o");
  xlabel("# of sub carriers");  ylabel("SNR [dB]");
  title("Estimated SNR");
otherwise
end;

switch wPLOTeachSub
    case  {"Yes" "yes" "YES" "y" "Y"}
      figure(9);  %% サブキャリアごとのコンスタレーションマップ．総ビット数が少ないと欠けが生じる
      clf;
      kkxL = 5;
      kkyL = ceil(max(chkSubC)/(kkxL*kkg));
      for subi=min(chkSubC):kkg:max(chkSubC)
        subplot(kkyL, kkxL, ceil((subi-min(chkSubC)+1)/kkg));
        if ((subi+kkg-1) < max(chkSubC))
          subiEnd = subi+kkg-1;
        else
          subiEnd = max(chkSubC);
        end
        temp = sh(subi:subiEnd,StartSymOfBerCalc:StopSymOfBerCalc);
        plot(reshape(temp, 1,prod(size(temp))), '.');
      	xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ...
        axis square;
        axis nolabel;
        axis off;
        pause(0.01);
      end
      subplot(kkyL, kkxL, ceil(kkxL/2));
      title(sprintf("Consteration maps for each %d sub carriers ",kkg));
    otherwise;
  end;

	otherwise
##	hold on;
##	subplot(1,2,1), plot(real(s0(:,(BLEN-10):(BLEN-1))), imag(s0(:,(BLEN-10):(BLEN-1))), '.b'), xlabel('In-Phase'), ylabel('Quadrature'), ...
##    xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ...
##    set(gca, 'XTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), set(gca, 'YTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), ...
##    axis square;
##subplot(1,2,2), plot(real(sh(:,(BLEN-10):(BLEN-1))), imag(sh(:,(BLEN-10):(BLEN-1))), '.r'), xlabel('In-Phase'), ylabel('Quadrature'), ...
##    xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ...
##    set(gca, 'XTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), set(gca, 'YTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), ...
##    axis square;

end;


pause(0.1);

dataset = struct;
dataset.ofdmVersion = ofdmVersion;
dataset.timestamp = datestr(now);
dataset.seed = SEEDs;
dataset.fntxdata = fntxdata;
dataset.fnadc = fnadc;
dataset.wExp = wExp;
dataset.corrPos = corrPos;
dataset.N = N;
dataset.P = P;
dataset.chkSubC = chkSubC;
dataset.Nq = Nq;
dataset.OvRate = OvRate;
dataset.BLEN = BLEN;
dataset.StartSymOfBerCalc = StartSymOfBerCalc;
dataset.StopSymOfBerCalc = StopSymOfBerCalc;
dataset.LFEQ = LFEQ;
dataset.L = L;
dataset.Lc = Lc;
dataset.Nf = Nf;
dataset.Ms = Ms;
dataset.fs = fs;
dataset.fsr= fsr;
dataset.fdac = fs/OvRate;
dataset.dt = dt;
dataset.Ns = Ns;
dataset.Dtx = Dtx;
dataset.Drx = Drx;
dataset.xDrx = xDrx;
dataset.DataRate = (sum(log2(Ms(chkSubC)))/((P+N)*dt*OvRate));
dataset.Nd = Nd;
dataset.chk01 = chk01;
dataset.EVMrms = EVMrms*100;
dataset.SNR = -20*log10(EVMrms);
dataset.PAPR_0 = PAPR_0;
dataset.PAPR_1 = PAPR_1;
dataset.consteOrig = ps0;
dataset.consteRcv  = psh;
switch wConsteNonly
	case {"Yes" "yes" "YES" "y" "Y"}
dataset.feqXSubC = fig4x;
dataset.feqY = fig4y/min(fig4y);
dataset.pdfX = fig6x;
dataset.pdfY = fig6y;
dataset.evmXSubC = fig7x;
dataset.evmY = fig7y;
dataset.snrXSubC = fig8x;
dataset.snry = fig8y;
endswitch
dataset.s0 = s0;
dataset.sh = sh;
dataset.uu0 = uu0;
dataset.rr = rr;
dataset.xDrx = xDrx;
dataset.rrs0 = rrs0;
dataset.lag = lag;
dataset.R = R;
dataset.idx = idx;
dataset.pks = pks;


[tDIR, tNAME, tEXT] = fileparts(fnadc);
%%save("-mat", sprintf("ofdmVer%02d-%04d%02d%02d-%02d%02d%02.0f.mat", ofdmVersion,datevec(now)), "dataset");
%save("-mat", sprintf("ofdmVer%02d-%s.mat", ofdmVersion,tNAME), "dataset");

switch wPLOTeachSub
  case  {"Yes" "yes" "YES" "y" "Y"}
    figL = 9;
  otherwise
    figL = 8;
end;

%%for k=1:figL
%%  figure(k);
%%  ss=sprintf("figure%d-%s.jpg",k, tNAME);
%%  print("-djpg", ss);
%%end;

switch wExp
	case {"Yes" "yes" "YES" "y" "Y"}
  switch wTxdata
	  case {"Yes" "yes" "YES" "y" "Y"}
      disp("tx2comb_tx.exe...");
      system(sprintf("tx2comb_tx.exe %s", fntxdata));
    otherwise
  end;
  otherwise
end;

switch wTxdata
  case {"Yes" "yes" "YES" "y" "Y"}
  otherwise
  switch wExp
	  case {"Yes" "yes" "YES" "y" "Y"}
%%      close(Hwarndlg);
    otherwise
  end;
end;

disp(sprintf('include plot: %f[sec]', toc));

figure(6); hold on; load -mat "Ysp.mat"; plot(fsp*1e-6,Ysp-max(Ysp)); hold off; xlim([0 100]); ylim([-40 0]);
