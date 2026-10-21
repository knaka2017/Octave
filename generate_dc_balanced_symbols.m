function [symIdx, symbols] = generate_dc_balanced_symbols(M, nSym, varargin)
  %------------------------------------------------------------
  % Generate non-reversible DC-balanced symbol sequence for PAM
  %
  % M      : modulation order (e.g. 4 or 8)
  % nSym   : number of symbols
  %
  % Optional name-value pairs:
  %   'WindowLen'       : local average control window
  %   'NumCandidates'   : number of random candidates per symbol
  %   'MaxRun'          : max allowed consecutive identical symbols
  %
  % This function is intended for PAM only.
  % It generates a symbol sequence that:
  %   - suppresses long runs
  %   - suppresses local DC bias
  %   - avoids overly periodic patterns
  %------------------------------------------------------------

  p = inputParser();
  addParameter(p, 'WindowLen', 16);
  addParameter(p, 'NumCandidates', 8);
  addParameter(p, 'MaxRun', 2);
  parse(p, varargin{:});

  windowLen = p.Results.WindowLen;
  numCandidates = p.Results.NumCandidates;
  maxRun = p.Results.MaxRun;

  if M < 2 || mod(log2(M), 1) ~= 0
    error('M must be a power of 2.');
  end

  % Ideal PAM levels: e.g. M=4 -> [-3 -1 +1 +3]
  levels = (-(M-1):2:(M-1)).';

  symIdx = zeros(nSym, 1);
  symbols = zeros(nSym, 1);

  for n = 1:nSym
    bestCost = Inf;
    bestIdx = randi([0, M-1]);
    bestSym = levels(bestIdx + 1);

    for k = 1:numCandidates
      candIdx = randi([0, M-1]);
      candSym = levels(candIdx + 1);

      %--------------------------------------------------------
      % Cost 1: local mean suppression
      %--------------------------------------------------------
      i1 = max(1, n - windowLen + 1);
      prev = symbols(i1:n-1);
      tmp = [prev; candSym];
      costMean = abs(mean(tmp));

      %--------------------------------------------------------
      % Cost 2: identical symbol run suppression
      %--------------------------------------------------------
      runCount = 0;
      t = n - 1;
      while t >= 1 && symbols(t) == candSym
        runCount = runCount + 1;
        t = t - 1;
      end

      if runCount >= maxRun
        costRun = 1000;   % strongly penalize
      else
        costRun = runCount;
      end

      %--------------------------------------------------------
      % Cost 3: weak penalty for too-regular sign alternation
      % to avoid strict periodic patterns like + - + - + -
      %--------------------------------------------------------
      costAlt = 0;
      if n >= 3
        s1 = symbols(n-1);
        s2 = symbols(n-2);
        if sign(candSym) == sign(s2) && sign(s1) ~= sign(s2)
          costAlt = 0.3;
        end
      end

      %--------------------------------------------------------
      % Cost 4: weak amplitude-balance penalty
      %--------------------------------------------------------
      if n >= 2
        i2 = max(1, n - windowLen + 1);
        recent = symbols(i2:n-1);
        costAmp = abs(mean(abs([recent; candSym])) - mean(abs(levels)));
      else
        costAmp = 0;
      end

      cost = 1.0 * costMean + 1.5 * costRun + 1.0 * costAlt + 0.2 * costAmp;

      if cost < bestCost
        bestCost = cost;
        bestIdx = candIdx;
        bestSym = candSym;
      end
    end

    symIdx(n) = bestIdx;
    symbols(n) = bestSym;
  end
end
