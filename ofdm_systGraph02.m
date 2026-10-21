%% display graph of OFDM_syst
%% 1. create folder
%% 2. copy some files of OFDM_syst results to the folder
%% 3. exec. this code
%% 4. select the folder (*** go into the folder in dialog. don't select file)
%% Note:
%% a) While this code is running, the folder changes. 
%%    If the prog. stops due to an error, return to the previous folder.
%%    ex.: "cd .."
clear;

wPLOTeachSub="n";   % サブキャリアごとのコンステ表示（数十秒かかるので注意）
kkg = 20;     % kkg個のサブキャリアごとに平均（サブキャリア毎コンステグラフ用）
figL = 10;

rdir = uigetdir(pwd, "Go into the folder of OFDM_syst results");

cd(rdir);
filelist = dir("ofdm*.mat");
fileL = length(filelist);
cd ..;

close all;

for k=1:figL
  figure(k);
  clf;
end

rrs0 = 0;
fig7y = 0;
berr = 0;
Nd = 0;

for k=1:fileL
  cd(rdir);
  load("-mat", filelist(k).name);
  cd ..;
  N = dataset.N;
  P = dataset.P;
  BLEN = dataset.BLEN;
  chkSubC = dataset.chkSubC;
  StartSymOfBerCalc = dataset.StartSymOfBerCalc;
  StopSymOfBerCalc = dataset.StopSymOfBerCalc;
  Ms = dataset.Ms;
  s0 = dataset.s0;
  sh = dataset.sh;
  uu0 = dataset.uu0;
  rr = dataset.rr;
  xDrx = dataset.xDrx;
  rrs0 = (rrs0*(k-1) + dataset.rrs0)/k;
  R = dataset.R;
  idx = dataset.idx;
  lag = dataset.lag;
  pks = dataset.pks;
  berr = berr + dataset.chk01*dataset.Nd;
  Nd = Nd + dataset.Nd;
  fs = dataset.fs;
  fsr= dataset.fsr;
  OvRate = dataset.OvRate;
  figure(1);  %% コンスタレーションマップ：元データ信号
    hold on;
    ps0 = s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc);  ps0 = reshape(ps0, prod(size(ps0)), 1);
    fig1x = real(ps0);  fig1y = imag(ps0);
    H=plot(fig1x, fig1y, '.b');
    title('Original signal');
    xlabel('In-Phase'), ylabel('Quadrature'), ...
    xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.5), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.5), ...
    set(H, 'MarkerSize', 8);
    set(gca, 'XTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), set(gca, 'YTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), ...
    axis square;
    hold off;
  figure(2);  %% コンスタレーションマップ：復調結果
    hold on;
    psh = sh(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc);  psh = reshape(psh, prod(size(psh)), 1);
    fig2x = real(psh);  fig2y = imag(psh);
    H=plot(fig2x, fig2y, '.r');
    title('demod. signal');
    xlabel('In-Phase'), ylabel('Quadrature'), ...
    xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.5), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.5), ...
    set(H, 'MarkerSize', 8);
    set(gca, 'XTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), set(gca, 'YTick', -(sqrt(max(Ms))-1):2:(sqrt(max(Ms))-1)), ...
    axis square;
    hold off;
  figure(3);  %% 送信信号・受信信号（一部），FEQの絶対値・偏角の平均（全サブキャリア）の時間変化
    subplot(2,2,1), plot(1:max(size(uu0(((N+P)*(BLEN-4)):((N+P)*(BLEN-1))))), real(uu0(((N+P)*(BLEN-4)):((N+P)*(BLEN-1))))), xlabel('k'), ylabel('u(k)');
    subplot(2,2,2), plot(1:max(size(rr(((N+P)*(BLEN-4)):((N+P)*(BLEN-1))))), real(rr(((N+P)*(BLEN-4)):((N+P)*(BLEN-1))))), xlabel('k'), ylabel('r(k)');
    ixDrx = StartSymOfBerCalc:StopSymOfBerCalc;
    subplot(2,2,3), plot(ixDrx,abs(mean(xDrx))); ylabel('abs of FEQ'); xlabel('index of blocks');
    subplot(2,2,4), plot(ixDrx,unwrap(angle(mean(xDrx)))*180/pi); xlabel('index of blocks'); ylabel('angle of FEQ [deg.]');
  figure(4);
    fig4x = chkSubC;  fig4y = mean(abs(xDrx(chkSubC,:).'));
    plot(fig4x, fig4y/min(fig4y));
    title("normalized avg. of abs(FEQ)");
    xlabel("# of Sub-Carrier");
    ylabel("norm. avg. of abs(FEQ coeff.)");
  figure(5);  %% 相関波形（鋭いピークが3本以上出ていないとNG）
    plot(lag, abs(R), "-", lag(idx), pks, "o");
    title(sprintf("Corr.(tx and rx signals) (with %d peaks)", length(idx)));
    xlabel("lag [samples]"); ylabel("Correlation");
  figure(6);
    Lrrs0 = length(rrs0);
    RRS0 = fft(rrs0.*hanning(Lrrs0)');
    tempRRS0 = 10*log10((1000/50)*(abs(RRS0).^2)/Lrrs0);
    fig6x = (0:(Lrrs0-1))*(fs/Lrrs0)*1e-6;  fig6y = tempRRS0-max(tempRRS0);
    plot(fig6x, fig6y);
    xlabel("Freq. [MHz]");  ylabel("Normalized Power [dB]");
    xlim([0 500]);  ylim([-60 0]);
    title("Power Spectrum of Received Signal (w/ hanning)");
  figure(7);
    fig7x = chkSubC;
    temp = sqrt( ( mean( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) - sh(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) )'.^2 ) ) / mean(mean( abs( s0(chkSubC,StartSymOfBerCalc:StopSymOfBerCalc) ).^2 )) )*100;
    fig7y = (fig7y*(k-1) + temp)/k;
    stem(fig7x, fig7y);
    xlabel("# of sub carriers");  ylabel("EVM [%]");
    title("EVMs of each sub carriers");
  figure(8);
    fig8x = chkSubC;
    fig8y = -20*log10(fig7y/100);
    plot(fig8x, fig8y, "o");
    xlabel("# of sub carriers");  ylabel("SNR [dB]");
    title("Estimated SNR");
switch wPLOTeachSub
    case  {"Yes" "yes" "YES" "y" "Y"}
  figure(9);  %% サブキャリアごとのコンスタレーションマップ．総ビット数が少ないと欠けが生じる
    hold on;
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
      plot(reshape(temp, 1,prod(size(temp))), '.b');
      xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.2), ...
      axis square;
      axis nolabel;
      axis off;
      pause(0.01);
    end
    subplot(kkyL, kkxL, ceil(kkxL/2));
    title(sprintf("Consteration maps for each %d sub carriers ",kkg));
    hold off;
  otherwise;
end;
  figure(10);
    hold on;
    UMs = unique(Ms(chkSubC));
    LUMs = length(UMs);
    kkxL = ceil(sqrt(LUMs));  kkyL = kkxL;
    for l = 1:LUMs
      subplot(kkxL, kkyL, l);
      temp = [];
      for subi=min(chkSubC):max(chkSubC)
        if (Ms(subi) == UMs(l))
          temp = [temp, sh(subi, StartSymOfBerCalc:StopSymOfBerCalc)];
        end
      end
      plot(temp, '.b');
      xlabel('In-Phase'), ylabel('Quadrature'), ...
      xlim([-(sqrt(UMs(l))-1) (sqrt(UMs(l))-1)]*1.2);
      ylim([-(sqrt(UMs(l))-1) (sqrt(UMs(l))-1)]*1.2);
      axis square;
%      axis nolabel;
%      axis off;
      set(gca, 'FontSize', 8);
      set(gca, 'XTick', -(sqrt(UMs(l))-1):2:(sqrt(UMs(l))-1)), set(gca, 'YTick', -(sqrt(UMs(l))-1):2:(sqrt(UMs(l))-1)), ...
      if (UMs(l) == 4)
        title("QPSK");
      else
        title(sprintf("%d-QAM", UMs(l)));
      end;
      pause(0.01);
    end
%    subplot(kkyL, kkxL, ceil(kkxL/2));
%    xlim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.5), ylim([-(sqrt(max(Ms))-1) (sqrt(max(Ms))-1)]*1.5), ...
%    axis square;
%    title(sprintf("Consteration maps for each QAM level ",kkg));
    hold off;
end

disp(sprintf("# of Sub carriers: %d", N));
disp(sprintf("CP length: %d Samples", P));
disp(sprintf("fs@Tx: %f Msps", dataset.fs*1e-6));
disp(sprintf("fs@Rx: %f Msps", dataset.fsr*1e-6));
disp(sprintf("Over Samples: %d", dataset.OvRate));
disp(sprintf('Data Rate = %f Mbps', (sum(log2(Ms(chkSubC)))/((P+N)*dataset.dt*dataset.OvRate))*1e-6));
disp(sprintf('OFDM sample speed = %f Msps', dataset.fdac*1e-6));
disp(sprintf("BER of total: %g", berr/Nd));
disp(sprintf("EVMrms of total: %f %%", mean(fig7y)));
UMs = unique(Ms(chkSubC));
LUMs = length(UMs);
for l = 1:LUMs
  temp = [];
  for subi=min(chkSubC):max(chkSubC)
    if (Ms(subi) == UMs(l))
      temp = [temp, fig7y(subi-min(chkSubC)+1)];
    end
  end
  if (UMs(l) == 4)
    disp(sprintf("QPSK EVM: %f %%", mean(temp)));
  else
    disp(sprintf("%d-QAM EVM: %f %%", UMs(l), mean(temp)));
  end;
end;


cd(rdir);
for k=1:figL
  if (k != 9)
    figure(k);
    ss=sprintf("figure%d.jpg",k);
    print("-djpg", ss);
  else
    switch wPLOTeachSub
      case  {"Yes" "yes" "YES" "y" "Y"}
        figure(k);
        ss=sprintf("figure%d.jpg",k);
        print("-djpg", ss);
      otherwise;
    end;
  end;
end;
cd ..;
