function cfg = default_cfg()

  cfg = struct();

  % random seeds
  cfg.txRandomSeed = 11;

  % mode
  cfg.mode = 'simulation';   % 'simulation' or 'experiment'

  % modulation
  cfg.modType = 'qam';       % 'qam' or 'pam'
  cfg.M = 16;                % modulation order

  % non-reversible DC-balanced test bit generation
  cfg.useDcBalancedBits = true;
  cfg.runBlockLen = 5;

  % Non-reversible DC-balanced PAM symbol generation
  cfg.useDcBalancedSymbols = true;
  cfg.dcBalanceWindowLen = 16;
  cfg.dcBalanceNumCandidates = 8;
  cfg.dcBalanceMaxRun = 2;

  % sampling / symbol rate
  cfg.fst = 2e9;             % Tx sampling frequency [Hz]
  cfg.fsr = 2e9;             % Rx sampling frequency [Hz]
  cfg.Rs  = 10e6;            % symbol rate [symbols/s]

  % exact Tx frame length in samples
  cfg.N = floor(0.8*(2^20));              % fixed exact length

  % ------------------------------------------------------------
  % Pulse shape
  %   'rect'         : rectangular pulse
  %   'sine'         : sinusoidal pulse within one symbol
  %   'rc'           : raised cosine
  %   'rrc'          : root raised cosine
  %   'supergaussian': super-Gaussian
  %   'rz'           : return-to-zero rectangular pulse
  %   'manchester'   : Manchester-shaped zero-mean pulse
  % ------------------------------------------------------------
  cfg.pulseShape = 'rect';

  % for sine pulse
  cfg.sineCyclesPerSymbol = 1;

  % for raised cosine / root raised cosine
  cfg.rcRollOff = 0.25;
  cfg.rcSpanSymbols = 6;

  % for RZ pulse
  % duty ratio in one symbol, 0 < duty <= 1
  cfg.rzDuty = 0.5;

  % for super-Gaussian
  cfg.sgSpanSymbols = 6;
  cfg.sgOrder = 4;
  cfg.sgBT = 0.5;

  % repeated transmission count
  cfg.repeatCount = 3;

  % simulation channel parameters
%  cfg.channel_h = [1.0; 0.025; -0.010; 0.005];
  cfg.channel_h = [1.0];
  cfg.delaySym = 2.5;        % propagation delay [symbols]
  cfg.SNRdB = 25;            % SNR for awgn()

  % synchronization
  cfg.peakSelectMode = 'index';   % 'largest' or 'index'
  cfg.peakSelectIndex = 2;
  cfg.minPeakDistanceFactor = 0.80;

  % manual sample offset after synchronization [samples at fsr]
  % positive: extract later
  % negative: extract earlier
  cfg.syncSampleOffset = 0;

  % experiment I/O
  cfg.txFile = 'tx_waveform.mat';
  cfg.rxFile = 'rx_waveform.mat';
  cfg.saveTxFile = true;
  cfg.tx_dataFile = 'tx_datafile.mat';

  % equalizer
  cfg.enableEqualizer = false;
  cfg.eqLen = 9;
  cfg.eqReg = 1e-6;
  cfg.eqDelay = 4;

  % plotting
  cfg.enablePlots = true;
  cfg.plotNumSymbols = 6;
  cfg.eyeNumTraces = 200;

  % Skip initial symbols for evaluation / plotting
  cfg.skipInitialSymbols = 20;

end

