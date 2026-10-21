function sync = synchronization_xcorr(rxAll, txRef, cfg)

  rxAll = rxAll(:);
  txRef = txRef(:);
  L_rx = length(rxAll);
  L_tx = length(txRef);

  if isempty(rxAll) || isempty(txRef)
    error('rxAll and txRef must be non-empty.');
  end

  if L_rx < L_tx
    error('length(rxAll) must be >= length(txRef).');
  end

  % ------------------------------------------------------------
  % 1. Fast Cross-correlation using FFT
  % ------------------------------------------------------------
  txEnergy = sum(abs(txRef).^2);
  if txEnergy < eps
    error('txRef energy is too small.');
  end

  Nfft = 2^nextpow2(L_rx + L_tx - 1);
  Rx_fft = fft(rxAll, Nfft);
  Tx_fft = fft(conj(flipud(txRef)), Nfft);
  y = ifft(Rx_fft .* Tx_fft);

  cValidRaw = y(L_tx : end);
  validLen = L_rx - L_tx + 1;
  cValidRaw = cValidRaw(1:validLen);

  cValid = cValidRaw / txEnergy;
  magc = abs(cValid);

  % ------------------------------------------------------------
  % 2. Peak detection
  % ------------------------------------------------------------
  minPeakDistance = max(1, round(cfg.minPeakDistanceFactor * L_tx));

  if isfield(cfg, 'peakQualityThreshold')
    threshold = cfg.peakQualityThreshold;
  else
    threshold = max(magc) * 0.1;
  end

  peaksAbs = [];
  peaksRaw = [];
  startIndices = [];
  magc_search = magc;

  max_peaks = ceil(L_rx / minPeakDistance) + 5;

  for k = 1:max_peaks
    [max_val, max_loc] = max(magc_search);

    if max_val < threshold || max_val == 0
      break;
    end

    peaksAbs(end+1) = max_val;
    peaksRaw(end+1) = cValid(max_loc);
    startIndices(end+1) = max_loc;

    idx_start = max(1, max_loc - minPeakDistance);
    idx_end   = min(length(magc_search), max_loc + minPeakDistance);
    magc_search(idx_start:idx_end) = 0;
  end

  if isempty(startIndices)
    error('No valid peaks found matching the criteria.');
  end

  peakLags = startIndices - 1;

  % ------------------------------------------------------------
  % Sort in time order
  % ------------------------------------------------------------
  [startIndicesSorted, order] = sort(startIndices, 'ascend');
  peakValsAbsSorted = peaksAbs(order);
  peakValsRawSorted = peaksRaw(order);
  peakLagsSorted = peakLags(order);

  % ------------------------------------------------------------
  % Peak selection
  % ------------------------------------------------------------
  switch lower(cfg.peakSelectMode)
    case 'largest'
      [~, idxSel] = max(peakValsAbsSorted);

    case 'index'
      if cfg.peakSelectIndex < 1 || cfg.peakSelectIndex > length(startIndicesSorted)
        error('cfg.peakSelectIndex is out of range.');
      end
      idxSel = cfg.peakSelectIndex;

    otherwise
      error('Unsupported peakSelectMode: %s', cfg.peakSelectMode);
  end

  selectedStartIndexEst = startIndicesSorted(idxSel);
  selectedLag = peakLagsSorted(idxSel);
  selectedPeakValueRaw = peakValsRawSorted(idxSel);
  selectedPeakValueAbs = peakValsAbsSorted(idxSel);

  % ------------------------------------------------------------
  % If selected correlation peak is negative, flip rx polarity
  % ------------------------------------------------------------
  polarityFlipped = false;
  rxAllCorrected = rxAll;

  if real(selectedPeakValueRaw) < 0
    rxAllCorrected = -rxAll;
    polarityFlipped = true;
  end

  % ------------------------------------------------------------
  % Manual synchronization offset is added here
  % ------------------------------------------------------------
  selectedStartIndex = selectedStartIndexEst + cfg.syncSampleOffset;

  % ------------------------------------------------------------
  % Simulation-mode debug
  % ------------------------------------------------------------
  if isfield(cfg, 'mode') && strcmpi(cfg.mode, 'simulation')
    trueLagFirst = round(cfg.delaySym * cfg.fsr / cfg.Rs);
    frameLenRx = round(cfg.N * cfg.fsr / cfg.fst);
    trueLags = trueLagFirst + (0:cfg.repeatCount-1) * frameLenRx;
    trueLags = trueLags(:);

    [~, idxNearest] = min(abs(trueLags - selectedLag));
    nearestTrueLag = trueLags(idxNearest);
    lagError = selectedLag - nearestTrueLag;

    fprintf('\n--- synchronization_xcorr debug (simulation mode) ---\n');
    fprintf('Detected peak lags (time order):\n');
    disp(peakLagsSorted.');

    fprintf('Detected normalized peak abs values:\n');
    disp(peakValsAbsSorted.');

    fprintf('True lag candidates:\n');
    disp(trueLags.');

    fprintf('Selected lag          = %d\n', selectedLag);
    fprintf('Nearest true lag      = %d\n', nearestTrueLag);
    fprintf('Lag error             = %d samples\n', lagError);
    fprintf('----------------------------------------------------\n\n');
  end

  % ------------------------------------------------------------
  % Output
  % ------------------------------------------------------------
  sync = struct();
  sync.selectedStartIndexEst = selectedStartIndexEst;
  sync.selectedStartIndex    = selectedStartIndex;
  sync.polarityFlipped       = polarityFlipped;
  sync.rxAllCorrected        = rxAllCorrected;
  sync.selectedPeakOrder     = idxSel;
  sync.selectedLag           = selectedLag;
  sync.selectedPeakValueRaw  = selectedPeakValueRaw;

end

