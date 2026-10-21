function w = local_tukeywin(N, alpha)
% local_tukeywin  - Tukey window (alpha fraction taper)
% w = local_tukeywin(N, alpha)
%   N     : window length
%   alpha : 0 = rectangular, 1 = Hann
  if nargin < 2, alpha = 0.5; end
  if alpha <= 0
      w = ones(N,1);
      return;
  elseif alpha >= 1
      n = (0:N-1).';
      w = 0.5*(1 - cos(2*pi*n/(N-1))); % Hann
      return;
  end
  w = ones(N,1);
  n = (0:N-1).';
  L = floor(alpha*(N-1)/2);
  if L > 0
      nl = (0:L).';
      w(1:L+1) = 0.5*(1 + cos(pi*(2*nl/(alpha*(N-1)) - 1)));
      nr = (0:L).';
      w(N-L:N) = 0.5*(1 + cos(pi*(2*nr/(alpha*(N-1)) + 1)));
  end
end
