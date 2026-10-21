function [evmPercent, evmRms, evmDb] = metrics_evm(ref, rx)

  % ------------------------------------------------------------
  % EVM on symbol layer
  % ref and rx must already be on the ideal constellation scale.
  %
  % Output:
  %   evmPercent : RMS EVM in percent
  %   evmRms     : RMS EVM (ratio)
  %   evmDb      : 20*log10(evmRms)
  % ------------------------------------------------------------

  ref = ref(:);
  rx  = rx(:);

  nUse = min(length(ref), length(rx));
  ref = ref(1:nUse);
  rx  = rx(1:nUse);

  if isempty(ref)
    error('metrics_evm: empty input.');
  end

  refPower = mean(abs(ref).^2);
  if refPower < eps
    error('metrics_evm: reference power is too small.');
  end

  err = rx - ref;
  errPower = mean(abs(err).^2);

  evmRms = sqrt(errPower / refPower);
  evmPercent = 100 * evmRms;
  evmDb = 20 * log10(max(evmRms, eps));

end

