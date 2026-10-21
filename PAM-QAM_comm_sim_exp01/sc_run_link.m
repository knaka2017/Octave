function result = sc_run_link(cfg)

  if nargin < 1
    cfg = default_cfg();
  end

  OvTx = round(cfg.fst / cfg.Rs);
  OvRx = round(cfg.fsr / cfg.Rs);

  if abs(cfg.fst / cfg.Rs - OvTx) > 1e-10
    warning('fst/Rs is not integer. Using OvTx = round(fst/Rs).');
  end

  if abs(cfg.fsr / cfg.Rs - OvRx) > 1e-10
    warning('fsr/Rs is not integer. Using OvRx = round(fsr/Rs).');
  end

  % Tx frame generation: exact N samples
  tx = tx_generate(cfg, OvTx);

  if length(tx.waveform) ~= cfg.N
    error('Tx waveform length is not equal to cfg.N.');
  end

  % synchronization reference
  % 同期用にはデータ部のみを使用
  txRefSyncTxRate = tx.waveform(1:tx.dataLen);
  txRefSyncRxRate = resample_to_fs(txRefSyncTxRate, cfg.fsr, cfg.fst);

  % full frame for extraction length
  txFrameFullRxRate = resample_to_fs(tx.waveform, cfg.fsr, cfg.fst);

  fprintf('len(tx.waveform)     = %d\n', length(tx.waveform));
  fprintf('len(txRefSync)       = %d\n', length(txRefSyncRxRate));
  fprintf('len(txFrameFullRef)  = %d\n', length(txFrameFullRxRate));

  % repeated Tx waveform into channel
  txRepeated = repmat(tx.waveform, cfg.repeatCount, 1);

  switch lower(cfg.mode)
    case 'simulation'
      ch = channel_simulation(cfg, txRepeated, cfg.fst);
      rxAll = resample_to_fs(ch.rxWaveform, cfg.fsr, cfg.fst);

      if cfg.saveTxFile
        tx_to_channel = txRepeated;
        fst = cfg.fst;
        fsr = cfg.fsr;
        Rs = cfg.Rs;
        modType = cfg.modType;
        M = cfg.M;
        repeatCount = cfg.repeatCount;
        N = cfg.N;
        save(cfg.txFile, 'tx_to_channel', 'fst', 'fsr', 'Rs', ...
             'modType', 'M', 'repeatCount', 'N');
        % for AWG (FPGA+DAC) 14bit, N=2^20
        tmp = tx.waveform;
        tmp = resample_signal(tmp, fst, 2.5e9, 2^20);
        tmp = (tmp - min(tmp))/(max(tmp) - min(tmp));
        tmp = floor(tmp * (-1 + 2^14));
        tx_data = tmp;
        save('-mat', cfg.tx_dataFile, 'tx_data');
      end

    case 'experiment'
      rxAll = load_rx_waveform(cfg.rxFile);

    otherwise
      error('Unknown cfg.mode: %s', cfg.mode);
  end

  rxAll = rxAll(:);
  fprintf('len(rxAll)           = %d\n', length(rxAll));

  %------------------------------------------------------------
  % Synchronization
  % - selectedStartIndex already includes cfg.syncSampleOffset
  % - rxAllCorrected is polarity-corrected if needed
  %------------------------------------------------------------
  sync = synchronization_xcorr(rxAll, txRefSyncRxRate, cfg);

  rxAllCorr = sync.rxAllCorrected;
  startIdxEst = sync.selectedStartIndexEst;
  startIdx = sync.selectedStartIndex;
  frameLen = length(txFrameFullRxRate);

  fprintf('Estimated start index = %d\n', startIdxEst);
  fprintf('Manual sync offset    = %+d samples\n', cfg.syncSampleOffset);
  fprintf('Corrected start index = %d\n', startIdx);
  fprintf('Selected peak value   = %+.6e\n', real(sync.selectedPeakValueRaw));

  if sync.polarityFlipped
    fprintf('Rx polarity           = flipped (rxAll <- -rxAll)\n');
  else
    fprintf('Rx polarity           = unchanged\n');
  end

  if startIdx < 1 || (startIdx + frameLen - 1) > length(rxAllCorr)
    error('Selected frame exceeds rxAll bounds.');
  end

  rxFrame = rxAllCorr(startIdx : startIdx + frameLen - 1);

  rx = rx_process(cfg, tx, rxFrame, OvRx);

  [nBitErr, nTotalBits, ber] = metrics_ber(rx.txSymIdxAligned, rx.decSymIdx, cfg.M);
  [evmPercent, evmRms, evmDb] = metrics_evm(rx.txSymbolsAligned, rx.eqSymbols);

  result = struct();
  result.cfg = cfg;
  result.tx = tx;
  result.sync = sync;
  result.rx = rx;
  result.metrics = struct();
  result.metrics.total_symbols = length(rx.txSymIdxAligned);
  result.metrics.total_bits = nTotalBits;
  result.metrics.bit_errors = nBitErr;
  result.metrics.ber = ber;
  result.metrics.evm_percent_rms = evmPercent;
  result.metrics.evm_rms = evmRms;
  result.metrics.evm_db = evmDb;
  result.metrics.skip_initial_symbols = rx.skipInitialSymbols;
  result.metrics.selected_start_index = startIdx;
  result.metrics.selected_peak_order = sync.selectedPeakOrder;
  result.metrics.selected_peak_lag = sync.selectedLag;
  result.metrics.selected_peak_value_raw = sync.selectedPeakValueRaw;
  result.metrics.polarity_flipped = sync.polarityFlipped;

  bitsPerSymbol = log2(cfg.M);
  bitRate = cfg.Rs * bitsPerSymbol;
  result.bitsPerSymbol = bitsPerSymbol;
  result.bitRate = bitRate;

  fprintf('\n============================================================\n');
  fprintf('Mode                  : %s\n', cfg.mode);
  fprintf('Modulation            : %s, M = %d\n', upper(cfg.modType), cfg.M);
  fprintf('Tx sampling freq fst  : %.6e Hz\n', cfg.fst);
  fprintf('Rx sampling freq fsr  : %.6e Hz\n', cfg.fsr);
  fprintf('Symbol rate Rs        : %.6e sym/s\n', cfg.Rs);
  fprintf('Bits/symbol           : %.2f\n', result.bitsPerSymbol);
  fprintf('Bit rate              : %.3f Mbps\n', result.bitRate/1e6);
  fprintf('OvTx                  : %d samples/symbol\n', OvTx);
  fprintf('OvRx                  : %d samples/symbol\n', OvRx);
  fprintf('Symbols/frame         : %d\n', tx.nSym);
  fprintf('Zero padding samples  : %d\n', tx.zeroPadLen);
  fprintf('Repeat count          : %d\n', cfg.repeatCount);
  fprintf('Skip initial symbols  : %d\n', rx.skipInitialSymbols);
  fprintf('Selected peak order   : %d\n', sync.selectedPeakOrder);
  fprintf('Estimated start index : %d\n', startIdxEst);
  fprintf('Manual sync offset    : %+d\n', cfg.syncSampleOffset);
  fprintf('Corrected start index : %d\n', startIdx);
  fprintf('Polarity flipped      : %d\n', sync.polarityFlipped);
  fprintf('Eval symbols          : %d\n', result.metrics.total_symbols);
  fprintf('Total bits            : %d\n', nTotalBits);
  fprintf('Bit errors            : %d\n', nBitErr);
  fprintf('BER                   : %.6e\n', ber);
  fprintf('RMS EVM               : %.4f %% (%.4f dB)\n', evmPercent, evmDb);
  fprintf('============================================================\n\n');

  if cfg.enablePlots
    plot_results(cfg, tx, txFrameFullRxRate, rxAllCorr, rxFrame, rx, OvRx, result.metrics);
  end

end

