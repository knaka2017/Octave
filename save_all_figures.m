function save_all_figures(fmt, dpi)
  % save_all_figures(fmt, dpi)
  % Save all currently opened figures with sequential filenames
  % If any target file already exists, retry from keyword input
  %
  % Usage:
  %   save_all_figures()              % PNG, 600 dpi
  %   save_all_figures('png', 300)    % PNG, 300 dpi
  %   save_all_figures('emf')         % Try EMF w/o dpi (environment dependent)

  if nargin < 1 || isempty(fmt)
    fmt = 'png';
  end
  if nargin < 2 || isempty(dpi)
    dpi = 600;
  end

  % Export font size (used only during PNG export)
  export_fontsize = 10;

  % ----------------------------
  % Get all figures
  % ----------------------------
  figs = findall(0, 'type', 'figure');

  if isempty(figs)
    fprintf('No figures found to save.\n');
    return;
  end

  % Get and sort figure numbers
  fig_nums = zeros(size(figs));
  for k = 1:numel(figs)
    fig_nums(k) = get(figs(k), 'number');
  end
  [fig_nums_sorted, idx] = sort(fig_nums);
  figs = figs(idx);

  % Number of digits for filename numbering
  max_fig_num = max(fig_nums_sorted);
  ndigits = max(2, length(num2str(max_fig_num)));

  % ----------------------------
  % Keyword input loop
  % ----------------------------
  while true

    keyword = '';
    try
      answer = inputdlg('Enter filename keyword (leave empty to use timestamp):', ...
                        'Save Figures', 1, {''});
      if isempty(answer)
        keyword = '';
      else
        keyword = strtrim(answer{1});
      end
    catch
      keyword = strtrim(input('Enter filename keyword (leave empty to use timestamp): ', 's'));
    end

    % Use timestamp if empty
    if isempty(keyword)
      keyword = strftime('%Y%m%d%H%M', localtime(time()));
      fprintf('No keyword provided. Using timestamp: %s\n', keyword);
    end

    % Generate filenames
    filenames = cell(numel(figs), 1);
    for k = 1:numel(figs)
      fig_num = fig_nums_sorted(k);
      numstr = sprintf(['%0', num2str(ndigits), 'd'], fig_num);
      filenames{k} = [keyword, numstr, '.', lower(fmt)];
    end

    % Check existing files
    exist_flag = false;
    for k = 1:numel(filenames)
      if exist(filenames{k}, 'file')
        exist_flag = true;
        break;
      end
    end

    if exist_flag
      fprintf('\nWARNING: One or more target files already exist.\n');
      fprintf('Please enter a different keyword.\n\n');
      continue;
    else
      break;
    end
  end

  fprintf('Saving %d figures as %s files...\n', numel(figs), upper(fmt));

  % ----------------------------
  % Save figures
  % ----------------------------
  for k = 1:numel(figs)
    fig = figs(k);
    filename = filenames{k};

    figure(fig);
    drawnow();

    switch lower(fmt)
      case 'png'

        % --- Save original font sizes ---
        objs = findall(fig, '-property', 'FontSize');
        original_sizes = zeros(size(objs));

        for i = 1:numel(objs)
          original_sizes(i) = get(objs(i), 'FontSize');
          set(objs(i), 'FontSize', export_fontsize);
        end

        % --- Save ---
        try
          print(fig, filename, '-dpng', ['-r', num2str(dpi)]);
        catch err
          % Restore on error
          for i = 1:numel(objs)
            set(objs(i), 'FontSize', original_sizes(i));
          end
          rethrow(err);
        end

        % --- Restore original font sizes ---
        for i = 1:numel(objs)
          set(objs(i), 'FontSize', original_sizes(i));
        end

      case 'emf'
        try
          print(fig, filename, '-dmeta');
        catch
          error(['Failed to save as EMF (environment dependent). ', ...
                 'Please use PNG instead.']);
        end

      otherwise
        error('Unsupported format. Use ''png'' or ''emf''.');
    end

    fprintf('Saved: %s\n', filename);
  end

  fprintf('All figures (%d) have been saved successfully.\n', numel(figs));
end

