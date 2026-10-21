function h = cosdesign_octave(beta, span, sps, shape)
  % Octave-compatible RC / RRC pulse generator
  % shape: 'rc' or 'rrc'

  if nargin < 4
    shape = 'rc';
  end

  N = span * sps;
  t = (-N/2:N/2).' / sps;   % in symbol durations

  h = zeros(size(t));

  switch lower(shape)
    case 'rc'
      for k = 1:length(t)
        tk = t(k);

        if abs(tk) < 1e-12
          h(k) = 1.0;

        elseif beta > 0 && abs(abs(2*beta*tk) - 1) < 1e-9
          h(k) = (pi/4) * sinc_scalar(1/(2*beta));

        else
          num = sinc_scalar(tk) * cos(pi*beta*tk);
          den = 1 - (2*beta*tk)^2;
          h(k) = num / den;
        end
      end

    case 'rrc'
      for k = 1:length(t)
        tk = t(k);

        if abs(tk) < 1e-12
          h(k) = 1 - beta + 4*beta/pi;

        elseif beta > 0 && abs(abs(4*beta*tk) - 1) < 1e-9
          h(k) = (beta/sqrt(2)) * ...
                 ((1 + 2/pi) * sin(pi/(4*beta)) + ...
                  (1 - 2/pi) * cos(pi/(4*beta)));

        else
          num = sin(pi*tk*(1-beta)) + 4*beta*tk*cos(pi*tk*(1+beta));
          den = pi*tk*(1 - (4*beta*tk)^2);
          h(k) = num / den;
        end
      end

    otherwise
      error('Unsupported shape: %s', shape);
  end
end

