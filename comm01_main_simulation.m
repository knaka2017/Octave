pkg load signal
pkg load communications

cfg = default_cfg();
cfg.mode = 'simulation';

% ------------------------------------------------------------
% Basic parameters
% ------------------------------------------------------------
cfg.modType = 'pam';     % 'qam' or 'pam'
cfg.M = 2;
cfg.fst = 2e9;
cfg.fsr = 2e9;
cfg.Rs  = 1e6;
cfg.N = floor(0.8*(2^20));
cfg.repeatCount = 1;

cfg.channel_h = [1.0; 0.30; -0.12; 0.05];
cfg.delaySym = 2.5;
cfg.SNRdB = 25;

cfg.peakSelectMode = 'largest';
% cfg.peakSelectMode = 'index';
% cfg.peakSelectIndex = 2;

% ------------------------------------------------------------
% Pulse shape selection by switch
% ------------------------------------------------------------
pulseMode = 7;
% 1: rectangular pulse
% 2: sinusoidal pulse
% 3: raised cosine
% 4: super-Gaussian
% 5: root raised cosine
% 6: RZ pulse
% 7: Manchester-shaped pulse

switch pulseMode
  case 1
    cfg.pulseShape = 'rect';

  case 2
    cfg.pulseShape = 'sine';
    cfg.sineCyclesPerSymbol = 2;   % integer number of cycles per symbol

  case 3
    cfg.pulseShape = 'rc';
    cfg.rcRollOff = 0.9;
    cfg.rcSpanSymbols = 6;

  case 4
    cfg.pulseShape = 'supergaussian';
    cfg.sgSpanSymbols = 6;
    cfg.sgOrder = 4;
    cfg.sgBT = 0.5;

  case 5
    cfg.pulseShape = 'rrc';
    cfg.rcRollOff = 0.1;
    cfg.rcSpanSymbols = 6;

  case 6
    cfg.pulseShape = 'rz';
    cfg.rzDuty = 0.5;

  case 7
    cfg.pulseShape = 'manchester';

  otherwise
    error('Unsupported pulseMode.');
end

result = sc_run_link(cfg);

save('last_cfg.mat','cfg');

