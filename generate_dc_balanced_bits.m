function bits = generate_dc_balanced_bits(nBits, runBlockLen)

  %------------------------------------------------------------
  % Generate non-reversible DC-balanced test bit sequence
  %
  % nBits        : total number of bits
  % runBlockLen  : block length used for run suppression
  %                e.g. 2 means 2-bit blocks such as 00 / 11
  %
  % This is NOT reversible coding.
  % It is intended only for test waveform generation.
  %------------------------------------------------------------

  if nargin < 2
    runBlockLen = 2;
  end

  if nBits < 1
    bits = [];
    return;
  end

  bits = zeros(nBits, 1);

  % first block: random
  firstLen = min(runBlockLen, nBits);
  bits(1:firstLen) = randi([0 1], firstLen, 1);

  pos = firstLen + 1;

  while pos <= nBits
    blkLen = min(runBlockLen, nBits - pos + 1);

    prevBlk = bits(pos - runBlockLen : pos - 1);
    prevMean = mean(prevBlk);

    % random candidate block
    cand = randi([0 1], blkLen, 1);

    % if previous block is all 1, bias next block toward 0
    % if previous block is all 0, bias next block toward 1
    if all(prevBlk == 1)
      if rand < 0.8
        cand(:) = 0;
      end
    elseif all(prevBlk == 0)
      if rand < 0.8
        cand(:) = 1;
      end
    end

    % additional run suppression:
    % if immediate previous bit run is too long, force inversion
    maxRun = runBlockLen;
    runVal = bits(pos - 1);
    runCount = 1;
    k = pos - 2;
    while k >= 1 && bits(k) == runVal
      runCount = runCount + 1;
      k = k - 1;
    end

    if runCount >= maxRun
      cand(:) = 1 - runVal;
    end

    bits(pos : pos + blkLen - 1) = cand;
    pos = pos + blkLen;
  end

end

