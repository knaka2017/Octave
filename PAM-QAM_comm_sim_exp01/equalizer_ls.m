function [yEq, w, txAligned, txIdxAligned] = equalizer_ls(rxSym, txSym, txIdx, L, D, reg)

  rxSym = rxSym(:);
  txSym = txSym(:);
  txIdx = txIdx(:);

  if L < 1
    yEq = rxSym;
    w = 1;
    txAligned = txSym;
    txIdxAligned = txIdx;
    return;
  end

  N = length(rxSym);

  if N <= L + D
    warning('Too few symbols for equalizer. Bypass equalizer.');
    yEq = rxSym;
    w = 1;
    txAligned = txSym;
    txIdxAligned = txIdx;
    return;
  end

  nStart = D + 1;
  nEnd   = N - L + 1;

  numRows = nEnd - nStart + 1;
  U = zeros(numRows, L);
  d = zeros(numRows, 1);
  idxKeep = zeros(numRows, 1);

  row = 1;
  for n = nStart:nEnd
    U(row, :) = rxSym(n-D : n-D+L-1).';
    d(row) = txSym(n);
    idxKeep(row) = n;
    row = row + 1;
  end

  w = (U' * U + reg * eye(L)) \ (U' * d);
  yEq = U * w;

  txAligned = txSym(idxKeep);
  txIdxAligned = txIdx(idxKeep);

end
