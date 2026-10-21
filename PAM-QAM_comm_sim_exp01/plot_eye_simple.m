function plot_eye_simple(x, sps, nTraces)

  x = x(:);

  span = 3 * sps;
  nAvail = floor(length(x) / sps) - 2;
  nTr = min(nTraces, max(0, nAvail));

  hold on;
  for k = 1:nTr
    idx1 = (k - 1) * sps + 1;
    idx2 = idx1 + span - 1;
    if idx2 <= length(x)
      plot(0:span-1, x(idx1:idx2), 'b');
    end
  end
  grid on;
  xlabel('Sample within 3 symbols');
  ylabel('Amplitude');

end

