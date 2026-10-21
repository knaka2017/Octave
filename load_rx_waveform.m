function rx = load_rx_waveform(filename)

  if ~exist(filename, 'file')
    error('Rx file not found: %s', filename);
  end

  s = load(filename);
  names = fieldnames(s);

  if isempty(names)
    error('No variable found in Rx file: %s', filename);
  end

  if isfield(s, 'rx_from_channel')
    rx = s.rx_from_channel;
    return;
  end

  if isfield(s, 'rx_waveform')
    rx = s.rx_waveform;
    return;
  end

  rx = s.(names{1});

end

