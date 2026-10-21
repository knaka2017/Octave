function y = resample_to_fs(x, fs_out, fs_in)

  x = x(:);

  if abs(fs_out - fs_in) < 1e-12
    y = x;
    return;
  end

  [p, q] = rat(fs_out / fs_in, 1e-12);
  y = resample(x, p, q);

end

