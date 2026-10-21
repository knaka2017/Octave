function tx = tx_generate(cfg, OvTx)

  % ------------------------------------------------------------
  % Fix random seed for transmitted bit/symbol pattern
  % ------------------------------------------------------------
  if isfield(cfg, 'txRandomSeed')
    rng(cfg.txRandomSeed);
  end

  bitsPerSym = round(log2(cfg.M));

  if abs(log2(cfg.M) - bitsPerSym) > 1e-12
    error('cfg.M must be a power of 2.');
  end

  % ------------------------------------------------------------
  % Number of symbols that fit in exactly N samples
  % ------------------------------------------------------------
  nSym = floor(cfg.N / OvTx);

  if nSym < 1
    error('No symbol fits into N samples. Check cfg.N and OvTx.');
  end

  dataLen = nSym * OvTx;
  zeroPadLen = cfg.N - dataLen;

  % ------------------------------------------------------------
  % Symbol generation
  % Keep symbol values on ideal constellation scale.
  % ------------------------------------------------------------
  useDcBalancedSymbols = false;
  if isfield(cfg, 'useDcBalancedSymbols')
    useDcBalancedSymbols = cfg.useDcBalancedSymbols;
  end

  switch lower(cfg.modType)
    case 'pam'
      if useDcBalancedSymbols
        windowLen = 16;
        numCandidates = 8;
        maxRun = 2;

        if isfield(cfg, 'dcBalanceWindowLen')
          windowLen = cfg.dcBalanceWindowLen;
        end
        if isfield(cfg, 'dcBalanceNumCandidates')
          numCandidates = cfg.dcBalanceNumCandidates;
        end
        if isfield(cfg, 'dcBalanceMaxRun')
          maxRun = cfg.dcBalanceMaxRun;
        end

        [symIdx, symbols] = generate_dc_balanced_symbols( ...
          cfg.M, nSym, ...
          'WindowLen', windowLen, ...
          'NumCandidates', numCandidates, ...
          'MaxRun', maxRun);
      else
        symIdx = randi([0, cfg.M - 1], nSym, 1);
        symbols = pammod(symIdx, cfg.M);
      end

    case 'qam'
      symIdx = randi([0, cfg.M - 1], nSym, 1);
      symbols = qammod(symIdx, cfg.M);

    otherwise
      error('Unsupported modType: %s', cfg.modType);
  end

  symbols = symbols(:);
  symIdx = symIdx(:);

  % ------------------------------------------------------------
  % Keep transmitted bits derived from the actual transmitted symbols
  % ------------------------------------------------------------
  bitMat = de2bi(symIdx, bitsPerSym, 'left-msb');
  bits = reshape(bitMat.', [], 1);

  % ------------------------------------------------------------
  % Pulse shaping
  % Use selected pulse shape and oversampling factor OvTx.
  %
  % symUpsampled:
  %   one symbol every OvTx samples
  %
  % dataWaveformRaw:
  %   pulse-shaped waveform before power normalization
  % ------------------------------------------------------------
  [pulse, pulseGroupDelay] = generate_pulse_shape(cfg, OvTx);

  symUpsampled = zeros(nSym * OvTx, 1);
  symUpsampled(1:OvTx:end) = symbols;

  dataWaveformRaw = conv(symUpsampled, pulse, 'same');

  % ------------------------------------------------------------
  % Waveform normalization
  % Symbol values are NOT normalized.
  % Only waveform is normalized to average power 1.
  % ------------------------------------------------------------
  waveScale = sqrt(mean(abs(dataWaveformRaw).^2));
  if waveScale < eps
    error('Waveform normalization scale is too small.');
  end

  dataWaveform = dataWaveformRaw / waveScale;

  % ------------------------------------------------------------
  % Exact N-sample frame with zero padding at tail
  % ------------------------------------------------------------
  waveform = [dataWaveform; zeros(zeroPadLen, 1)];

  if length(waveform) ~= cfg.N
    error('Generated Tx waveform length does not match cfg.N.');
  end

  % ------------------------------------------------------------
  % Output
  % ------------------------------------------------------------
  tx = struct();

  % bit / symbol information
  tx.bits = bits;
  tx.bitsPerSym = bitsPerSym;
  tx.symIdx = symIdx;
  tx.symbols = symbols;                 % ideal constellation scale
  tx.nSym = nSym;

  % waveform information
  tx.waveform = waveform;               % normalized full frame waveform
  tx.dataWaveform = dataWaveform;       % normalized data part only
  tx.dataWaveformRaw = dataWaveformRaw; % before waveform normalization
  tx.symUpsampled = symUpsampled;
  tx.waveScale = waveScale;

  % pulse information
  tx.pulse = pulse;
  tx.pulseGroupDelay = pulseGroupDelay;
  tx.pulseShape = cfg.pulseShape;

  % frame information
  tx.OvTx = OvTx;
  tx.dataLen = dataLen;
  tx.zeroPadLen = zeroPadLen;

end

