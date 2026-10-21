function [yCorr, a, b] = affine_symbol_fit(rxSym, txSym, mode)

  rxSym = rxSym(:);
  txSym = txSym(:);

  nUse = min(length(rxSym), length(txSym));
  rxSym = rxSym(1:nUse);
  txSym = txSym(1:nUse);

  switch lower(mode)
    case 'pam'
      % Real affine fit: tx ≈ a * rx + b
      xr = real(rxSym);
      tr = real(txSym);

      U = [xr, ones(length(xr),1)];
      p = U \ tr;

      a = p(1);
      b = p(2);

      yCorr = a * xr + b;

    case 'qam'
      % Complex affine fit: tx ≈ a * rx + b
      U = [rxSym, ones(length(rxSym),1)];
      p = U \ txSym;

      a = p(1);
      b = p(2);

      yCorr = a * rxSym + b;

    otherwise
      error('Unsupported mode: %s', mode);
  end

end
