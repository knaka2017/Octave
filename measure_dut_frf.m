function frf = measure_dut_frf(z, x, fs_z, fs_x, varargin)
% measure_dut_frf_fast
%   必須は z, x, fs_z, fs_x のみ。
%   その他は既定値を自動設定し、estimate_dut_frf_fast を呼ぶ。
%
% Octave 互換版（MATLAB string を使用しない）

  % --- 必須チェック ---
  if nargin < 4
    error('Arguments required: z, x, fs_z, fs_x');
  end
  if isempty(z) || isempty(x)
    error('z and x must be non-empty vectors.');
  end
  if fs_z <= 0 || fs_x <= 0
    error('fs_z and fs_x must be > 0.');
  end

  % --- 既定値 ---
  params.fl           = 100e3;
  params.fh           = min(0.98e9, 0.98*(fs_z/2));   % ガード
  params.use_window   = true;
  params.max_cycles   = 4;
  params.decim_coarse = 32;
  params.fine_win     = 200000;
  params.loc_win      = 5000;
  params.use_single   = false;

  % --- Name-Value の上書き（Octave互換） ---
  if mod(numel(varargin), 2) ~= 0
      error('Name-Value pairs are malformed.');
  end

  for k = 1:2:numel(varargin)
      name = lower(char(varargin{k}));   % Octave 互換
      val  = varargin{k+1};

      switch name
          case 'fl'
              params.fl = val;
          case 'fh'
              params.fh = val;
          case 'use_window'
              params.use_window = logical(val);
          case 'max_cycles'
              params.max_cycles = val;
          case 'decim_coarse'
              params.decim_coarse = val;
          case 'fine_win'
              params.fine_win = val;
          case 'loc_win'
              params.loc_win = val;
          case 'use_single'
              params.use_single = logical(val);
          otherwise
              error('Unknown parameter name: %s', name);
      end
  end

  % --- 最終ガード（ナイキスト 0.98×） ---
  params.fh = min(params.fh, 0.98*(fs_z/2));

  fprintf('[simple] fl=%.3f kHz, fh=%.3f MHz, cycles=%d\n',...
          params.fl/1e3, params.fh/1e6, params.max_cycles);

  % --- 本体呼び出し ---
  frf = estimate_dut_frf( ...
          z, x, fs_z, fs_x, ...
          params.fl, params.fh, params.use_window, params.max_cycles, ...
          params.decim_coarse, params.fine_win, params.loc_win, params.use_single);
end
