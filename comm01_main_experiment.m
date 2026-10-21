pkg load signal
pkg load communications

load('last_cfg.mat');

cfg.mode = 'experiment';

cfg.txFile = 'tx_waveform.mat';
cfg.rxFile = 'rx_waveform.mat';

result = sc_run_link(cfg);

