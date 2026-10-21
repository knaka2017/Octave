function [nErr, nTotal, ber] = metrics_ber(txSymIdx, rxSymIdx, M)

  txSymIdx = txSymIdx(:);
  rxSymIdx = rxSymIdx(:);

  nUse = min(length(txSymIdx), length(rxSymIdx));
  txSymIdx = txSymIdx(1:nUse);
  rxSymIdx = rxSymIdx(1:nUse);

  k = round(log2(M));

  txBits = de2bi(txSymIdx, k, 'left-msb');
  rxBits = de2bi(rxSymIdx, k, 'left-msb');

  diffBits = (txBits ~= rxBits);

  nErr = sum(diffBits(:));
  nTotal = numel(txBits);
  ber = nErr / max(nTotal, 1);

end

