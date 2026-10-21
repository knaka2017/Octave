function y = sinc_scalar(x)
  if abs(x) < 1e-12
    y = 1.0;
  else
    y = sin(pi*x) / (pi*x);
  end
end
