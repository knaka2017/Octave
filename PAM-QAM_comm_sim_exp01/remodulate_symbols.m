function s = remodulate_symbols(symIdx, modType, M)

  symIdx = symIdx(:);

  switch lower(modType)
    case 'qam'
      s = qammod(symIdx, M);
    case 'pam'
      s = pammod(symIdx, M);
    otherwise
      error('Unsupported modType: %s', modType);
  end

end

