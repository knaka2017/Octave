function [pulse, groupDelay] = generate_pulse_shape(cfg, Ov)

  pulseType = lower(cfg.pulseShape);

  switch pulseType
    case 'rect'
      pulse = ones(Ov, 1);
      groupDelay = floor((length(pulse)-1)/2);

    case 'sine'
      % One-symbol pulse with integer cycles per symbol
      n = (0:Ov-1).';
      fmult = cfg.sineCyclesPerSymbol;
      pulse = sin(2*pi*fmult*n/Ov);
      groupDelay = floor((length(pulse)-1)/2);

    case 'rc'
      % Raised cosine pulse
      beta = cfg.rcRollOff;
      span = cfg.rcSpanSymbols;   % in symbols
      pulse = cosdesign_octave(beta, span, Ov, 'rc');
      groupDelay = floor((length(pulse)-1)/2);

    case 'rrc'
      % Root raised cosine pulse
      beta = cfg.rcRollOff;
      disp(beta);
      span = cfg.rcSpanSymbols;   % in symbols
      pulse = cosdesign_octave(beta, span, Ov, 'rrc');
      groupDelay = floor((length(pulse)-1)/2);

    case 'rz'
      % Return-to-zero rectangular pulse within one symbol
      duty = cfg.rzDuty;
      duty = min(max(duty, 1/Ov), 1.0);

      nOn = max(1, round(duty * Ov));
      pulse = [ones(nOn,1); zeros(Ov - nOn,1)];
      groupDelay = floor((length(pulse)-1)/2);

    case 'manchester'
      % Manchester-shaped pulse:
      % first half positive, second half negative
      % zero-mean over one symbol
      n1 = floor(Ov/2);
      n2 = Ov - n1;
      pulse = [ones(n1,1); -ones(n2,1)];
%      pulse = zeros(Ov,1); pulse(1:10)=1; pulse(11:20)=-1;
      groupDelay = floor((length(pulse)-1)/2);

    case 'supergaussian'
      % Super-Gaussian pulse over multiple symbols
      span = cfg.sgSpanSymbols;
      ord = cfg.sgOrder;
      BT = cfg.sgBT;

      n = (-span*Ov/2 : span*Ov/2).';
      t = n / Ov;   % in symbol durations

      sigma = max(1e-6, 0.5 / max(BT,1e-6));
      pulse = exp(-0.5 * (abs(t)/sigma).^(2*ord));
      groupDelay = floor((length(pulse)-1)/2);

    otherwise
      error('Unsupported pulseShape: %s', cfg.pulseShape);
  end

  pulse = pulse(:);

  % Normalize pulse energy
  e = sqrt(sum(abs(pulse).^2));
  if e < eps
    error('Pulse energy is too small.');
  end
  pulse = pulse / e;

end

