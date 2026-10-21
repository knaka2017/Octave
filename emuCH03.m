function flag = emuCH03(fntxdata, SNR, hFlag)
%% Emulate channel output using certain optical wireless comm. syst.
%%
%% flag = emuCH02(fntxdata, SNR, hFlag)
%%     fntxdata: file name of tx_data
%%       (vector name "tx_data"(different from file name),
%%        1x2^20 vector,
%%        integer:0-16383(2^14-1))
%%     SNR: SNR [dB] for awgn function
%%     hFlag: if 1, emulate channel
%%
%%     flag: error code
%%       1:error, 0:correct
%%       if flag is correct, save file of emurated received signal
%%       (emulate noise and freq. response of the channel)
%%       (1 x 262480, 2.5Gsample/s ADC output)
%%
%% saved "emuCH02.mat" with the same format of "fnadc" file
%% The file can be used as "fnadc" in "ofdm_syst" code.
%%   fntxdata: name of tx_data file before processing tx2comb_tx.exe
%%   flag: error code
pkg load signal;
pkg load communications;

  flag = 0;

  LyyyB = 262480;

  Vpp = 0.2;

  f = [0 0.88 3.35 4.19 5 10 15 20 30 40 50 75 100 122 175 1250]/(2500/2);
  f = [0 0.88 3.35 4.19 5 10 15 20 30 40 50 75 100 122 175 312.5]/(2500/2/4);
  f = [0 2.5 3.4 5 6.5 10 12.3 15 20 28.8 40 50 70 100 150 300 400 500 1250]/(2500/2);
  m = [0 0.006 1 1 0.552 0.245 0.085 0.045 0.015 0.0076 0.0049 0.0022 0.0014 0.0013 0 0];
  m = [0 0.006 1 1 0.8 0.7 0.55 0.43 0.35 0.26 0.18 0.1 0.07 0.04 0 0];
  m = [0 1 0.677 0.323 0.2 0.083 0.055 0.03 0.016 0.00716 0.00324 0.00215 0.00109 0.675e-3 0.424e-3 0.266e-3 0.2e-3 0.186e-3 0];
  B = fir2(2000, f, m);
%  B = fir1(10, 0.001);
%load "w_pam4ov01_07.mat";
%B = w;
  if (hFlag ~= 1)
    B = 1;
  end;

  if (exist(fntxdata) ~= 2)
    disp(sprintf("Not exist: %s ", fntxdata))
    disp('Please see the help (>>help emuCH01)');
    flag = 1;
    return;
  end;
  load("-mat", fntxdata);
  if (exist('tx_data') ~= 1)
    disp('CHECK ERROR: name of TX vector is not "tx_data"');
    flag = 1;
    return;
  end;
  if (size(tx_data,2) ~= 2^20 || size(tx_data,1) ~= 1 )
    disp('CHECK ERROR: size of TX vector is not 1 x 2^20');
    flag = 1;
    return;
  end;
  if (max(tx_data) > 16383)
    disp('CHECK ERROR: max(tx_data) > 16383(=2^14-1)');
    flag = 1;
    return;
  end;
  if (min(tx_data) < 0)
    disp('CHECK ERROR: min(tx_data) < 0');
    flag = 1;
    return;
  end;
  if (sum(abs(tx_data-floor(tx_data)))>0)
    disp('CHECK ERROR: tx_data is not integer');
    flag = 1;
    return;
  end;

  x = (tx_data/(2^14-1))*0.5-0.25;
  x = x - mean(x);
  x0= x;
  x = filter(B,1, x); x = (x/max(abs(x(1:(length(x)-1000)))))*(Vpp*0.5);
  x = awgn(x, SNR, "measured");  %% Tx noise
%  x = awgn(x, 20 + 32.59488);  %% Tx noise
  ds = 0.5/((2^8)-1);
  x = floor((x+ds*0.5)/ds);  x = (x/128)*0.25;

  startPos = randint(1,1,LyyyB-1000)+1000;
startPos = 1000;
  endPos = startPos+LyyyB-1;
  yyyB = -x(startPos:endPos);

%%  t = 0:(length(x)-1);  t = t/2.5e9;
%%  plot(t, x0, t, x);

  save("-mat", "emuCH02.mat", "yyyB");

  disp('CHECK OK');
  disp(sprintf('output "emuCH02.mat" (emurated received signal, "yyyB"(1x%d), 2.5GSa/s)', LyyyB));
end
