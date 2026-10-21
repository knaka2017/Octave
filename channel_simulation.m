function ch = channel_simulation(cfg, txRepeated, fs)

  delaySamples = round(cfg.delaySym * fs / cfg.Rs);

  delayedInput = [zeros(delaySamples, 1); txRepeated(:)];

  y = filter(cfg.channel_h(:), 1, delayedInput);
  y = awgn(y, cfg.SNRdB, 'measured');

  ch = struct();
  ch.delaySamples = delaySamples;
  ch.rxWaveform = y;

end

