function flag = emuCH01(fntxdata)
%% Emulate optical channel using certain optical wireless comm. syst.
%% save "emuCH01.mat" with the same format of "fnadc" file
%% The file can be used as "fnadc" in "ofdm_syst" code.
%%   fntxdata: name of tx_data file before processing tx2comb_tx.exe
%%   flag: error code

  flag = 0;

  LyyyB = 262480;

  Vpp = 0.2;

  f = [0 0.88 3.35 4.19 5 10 15 20 30 40 50 75 100 122 175 1250]/(2500/2);
  m = [0 0.006 1 1 0.552 0.245 0.085 0.045 0.015 0.0076 0.0049 0.0022 0.0014 0.0013 0 0];
  m = [0 0.006 1 1 0.8 0.7 0.55 0.43 0.35 0.26 0.18 0.1 0.07 0.04 0 0];
  B = fir2(2000, f, m);
%  B = fir1(10, 0.001);

  load("-mat", fntxdata);
  x = (tx_data/(2^14-1))*0.5-0.25;
  x = x - mean(x);
  x0= x;
%  x = awgn(x, 20, "measured");  %% Tx noise
  x = awgn(x, 20 + 32.59488);  %% Tx noise
  x = filter(B,1, x); x = (x/max(abs(x(1:(length(x)-1000)))))*(Vpp*0.5);
  ds = 0.5/((2^8)-1);
  x = floor((x+ds*0.5)/ds);  x = (x/128)*0.25;

  startPos = randint(1,1,LyyyB-1000)+1000;
startPos = 1000;
  endPos = startPos+LyyyB-1;
  yyyB = -x(startPos:endPos);

%%  t = 0:(length(x)-1);  t = t/2.5e9;
%%  plot(t, x0, t, x);

  save("-mat", "emuCH01.mat", "yyyB");
end
