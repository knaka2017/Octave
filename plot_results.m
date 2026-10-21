function plot_results(cfg, tx, txFrameRefRxRate, rxAll, rxFrame, rx, OvRx, metrics)

  nPlotSym = cfg.plotNumSymbols;

  nSkipSym = 5;
  if isfield(metrics, 'skip_initial_symbols')
    nSkipSym = metrics.skip_initial_symbols;
  end

  % ------------------------------------------------------------
  % Plot ranges after skipping initial symbols
  % ------------------------------------------------------------
  txStart = min(length(tx.waveform), nSkipSym * tx.OvTx + 1);
  txEnd   = min(length(tx.waveform), txStart + nPlotSym * tx.OvTx - 1);

  rxStart = min(length(rxFrame), nSkipSym * OvRx + 1);
  rxEnd   = min(length(rxFrame), rxStart + nPlotSym * OvRx - 1);

  if txEnd < txStart
    txStart = 1;
    txEnd = min(length(tx.waveform), nPlotSym * tx.OvTx);
  end

  if rxEnd < rxStart
    rxStart = 1;
    rxEnd = min(length(rxFrame), nPlotSym * OvRx);
  end

  % ------------------------------------------------------------
  % Time axes [ns]
  % ------------------------------------------------------------
  tTx_ns = ((txStart:txEnd) - 1) / cfg.fst * 1e6;
  tRx_ns = ((rxStart:rxEnd) - 1) / cfg.fsr * 1e6;

  % ------------------------------------------------------------
  % Tx waveform (physical waveform layer, normalized waveform)
  % ------------------------------------------------------------
  figure(1);
  plot(tTx_ns-min(tTx_ns) -0.5*1e6/cfg.Rs, real(tx.waveform(txStart:txEnd)), 'LineWidth', 1.2);
  grid on;
  set(gca, 'fontsize', 14);
  xlabel('Time [us]');
  ylabel('Voltage [V]');
  title(sprintf('Ideal signal'));
  xlim([0 max(tTx_ns-min(tTx_ns))]);

  % ------------------------------------------------------------
  % Extracted Rx waveform
  % ------------------------------------------------------------
  figure(2);
  plot(tRx_ns-min(tRx_ns)-0.5*1e6/cfg.Rs, real(rxFrame(rxStart:rxEnd)), 'LineWidth', 1.2);
  grid on;
  set(gca, 'fontsize', 14);
  xlabel('Time [us]');
  ylabel('Voltage');
  title(sprintf('Received signal', nSkipSym));
  xlim([0 max(tRx_ns-min(tRx_ns))]);

  % ------------------------------------------------------------
  % Tx/Rx waveform comparison
  % Both waveforms are on Rx-rate time axis
  % ------------------------------------------------------------
  figure(3);
  nShow = min(txEnd - txStart + 1, rxEnd - rxStart + 1);

  txCompStart = txStart;
  txCompEnd   = txStart + nShow - 1;
  rxCompStart = rxStart;
  rxCompEnd   = rxStart + nShow - 1;

  tComp_ns = ((0:nShow-1) + (rxCompStart - 1)) / cfg.fsr * 1e6;

  tmp = real(txFrameRefRxRate(txCompStart:txCompEnd));
  plot(tComp_ns-min(tComp_ns)-0.5*1e6/cfg.Rs, real(txFrameRefRxRate(txCompStart:txCompEnd))/max(abs(tmp)), 'LineWidth', 1.2);
  xlim([0 max(tRx_ns-min(tRx_ns))]);
  hold on;
  tmp = real(rxFrame(rxCompStart:rxCompEnd));
  plot(tComp_ns-min(tComp_ns)-0.5*1e6/cfg.Rs, real(rxFrame(rxCompStart:rxCompEnd))/max(abs(tmp)), 'LineWidth', 1.0);
  xlim([0 max(tRx_ns-min(tRx_ns))]);
  grid on;
  set(gca, 'fontsize', 14);
  xlabel('Time [us]');
  ylabel('Norm. Amplitude');
  legend('Ideal', 'Received');
  title(sprintf('Waveform comparison'));

  % ------------------------------------------------------------
  % Power Spectral Density (PSD)  (one-sided spectrum)
  % Tx: use txFrameRefRxRate (already on Rx sampling rate)
  % Rx: use extracted rxFrame
  % ------------------------------------------------------------

  txPsdSig = txFrameRefRxRate(:);
  rxPsdSig = rxFrame(:);

  nFftTx = 2^nextpow2(length(txPsdSig));
  nFftRx = 2^nextpow2(length(rxPsdSig));

  TX = fft(txPsdSig, nFftTx);
  RX = fft(rxPsdSig, nFftRx);

  psdTx = abs(TX).^2 / max(length(txPsdSig),1);
  psdRx = abs(RX).^2 / max(length(rxPsdSig),1);

  % one-sided
  psdTx = psdTx(1:nFftTx/2+1);
  psdRx = psdRx(1:nFftRx/2+1);

  % scale for one-sided spectrum (except DC and Nyquist)
  if length(psdTx) > 2
    psdTx(2:end-1) = 2 * psdTx(2:end-1);
  end

  if length(psdRx) > 2
    psdRx(2:end-1) = 2 * psdRx(2:end-1);
  end

  fTx_MHz = (0:nFftTx/2) / nFftTx * cfg.fsr / 1e6;
  fRx_MHz = (0:nFftRx/2) / nFftRx * cfg.fsr / 1e6;

  figure(4);
  tmp = 10*log10(psdTx + eps);
  plot(fTx_MHz, tmp - max(tmp), 'LineWidth', 1.2);
  grid on;
  set(gca, 'fontsize', 14);
  xlabel('Frequency [MHz]');
  ylabel('Norm. Power Spectral Density [dB]');
  title('Ideal signal PSD (one-sided)');

  figure(5);
  tmp = 10*log10(psdRx + eps);
  plot(fRx_MHz, tmp - max(tmp), 'LineWidth', 1.2);
  grid on;
  set(gca, 'fontsize', 14);
  xlabel('Frequency [MHz]');
  ylabel('Norm. Power Spectral Density [dB]');
  title('Received signal PSD (one-sided)');

  % ------------------------------------------------------------
  % Constellation / symbol scatter
  % tx.txSymbolsAligned and rx.eqSymbols are both on ideal symbol scale
  % ------------------------------------------------------------
##  figure(6);
##  plot(real(rx.txSymbolsAligned), imag(rx.txSymbolsAligned), 'o');
##  hold on;
##  plot(real(rx.eqSymbols), imag(rx.eqSymbols), '.');
##  plot(real(rx.decSymbols), imag(rx.decSymbols), 'x');
##  grid on;
##  set(gca, 'fontsize', 14);
##  axis equal;
##  xlabel('In-phase');
##  ylabel('Quadrature');
##  legend('Tx symbols (ideal)', 'Equalized Rx symbols', 'Decision symbols');
##
##  if isfield(metrics, 'evm_db')
##    title(sprintf('Constellation / BER = %.3e, EVM = %.2f %% (%.2f dB), skip = %d sym', ...
##          metrics.ber, metrics.evm_percent_rms, metrics.evm_db, nSkipSym));
##  else
##    title(sprintf('Constellation / BER = %.3e, EVM = %.2f %%, skip = %d sym', ...
##          metrics.ber, metrics.evm_percent_rms, nSkipSym));
##  end

  % ------------------------------------------------------------
  % Symbol sequence plots
  % These should be shown on ideal constellation scale.
  % ------------------------------------------------------------
  nStem = min(100, length(rx.txSymbolsAligned));

  figure(7);
  subplot(2,1,1);
  stem(real(rx.txSymbolsAligned(1:nStem)), 'filled');
  grid on;
  set(gca, 'fontsize', 14);
  xlabel('Symbol index');
  ylabel('Amplitude');
  title(sprintf('Ideal symbols', nStem));

  subplot(2,1,2);
  stem(real(rx.eqSymbols(1:nStem)), 'filled');
  grid on;
  set(gca, 'fontsize', 14);
  xlabel('Symbol index');
  ylabel('Amplitude');
  title(sprintf('Detected symbols @ Rx', nStem));

  % ------------------------------------------------------------
  % Eye diagram
  % ------------------------------------------------------------
##  figure;
##  plot_eye_simple(real(rxFrame(rxStart:end)), OvRx, cfg.eyeNumTraces);
##  set(gca, 'fontsize', 14);
##  title(sprintf('Eye diagram (real part, after skipping first %d symbols)', nSkipSym));

  % ------------------------------------------------------------
  % Received whole waveform and selected frame start
  % Octave-compatible vertical line drawing
  % ------------------------------------------------------------
##  figure;
##  tAll_ns = (0:length(rxAll)-1) / cfg.fsr * 1e9;
##  plot(tAll_ns, real(rxAll), 'LineWidth', 1.0);
##  hold on;
##  grid on;
##  set(gca, 'fontsize', 14);
##  xlabel('Time [ns]');
##  ylabel('Amplitude');
##  title('Received whole waveform and selected frame start (real part)');
##
##  yl = ylim;
##  xs = (metrics.selected_start_index - 1) / cfg.fsr * 1e9;
##  plot([xs xs], [yl(1) yl(2)], '--r', 'LineWidth', 1.5);
##
##  legend('Rx all', 'Selected frame start');

end

