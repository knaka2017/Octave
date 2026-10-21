function y = resample_signal(x, Fsi, Fso, Nout)
% =========================================================
% resample_signal.m
% サンプリングレート変換関数
%
% 使い方:
%   y = resample_signal(x, Fsi, Fso, Nout)
%
% 引数:
%   x    - 入力信号ベクトル（行ベクトルまたは列ベクトル、長さ Nin）
%   Fsi  - 入力サンプリング周波数 [Hz]
%   Fso  - 出力サンプリング周波数 [Hz]
%   Nout - 出力信号ベクトルの長さ（サンプル数）
%          Nout > リサンプル結果長 → 末尾をゼロパディング
%          Nout < リサンプル結果長 → 末尾を切り捨て
%
% 戻り値:
%   y    - 出力信号ベクトル（長さ Nout）
% =========================================================

  pkg load signal;

  % --- 入力チェック ---
  if ~isvector(x)
    error('x はベクトルでなければなりません。');
  end
  if Fsi <= 0 || Fso <= 0
    error('サンプリング周波数は正の値でなければなりません。');
  end
  if ~isscalar(Nout) || Nout < 1 || floor(Nout) ~= Nout
    error('Nout は正の整数でなければなりません。');
  end

  is_col = iscolumn(x);   % 列ベクトルかどうか記憶
  x      = x(:)';         % 行ベクトルに統一して処理

  % --- サンプリング周波数が同じ場合 ---
  if Fsi == Fso
    y_raw = x;
  else
    % --- 有理数比 p/q を求める（最大公約数で約分）---
    %   resample() は整数比 p/q を要求する
    %   浮動小数点誤差を避けるため、十分な精度で整数化してから gcd を取る
    precision = 1e6;
    p_raw = round(Fso * precision);
    q_raw = round(Fsi * precision);
    g     = gcd(p_raw, q_raw);
    p     = p_raw / g;   % アップサンプル係数
    q     = q_raw / g;   % ダウンサンプル係数

    % --- リサンプリング実行 ---
    %   signal パッケージの resample(x, p, q) を使用
    %   内部でポリフェーズFIRフィルタによるアンチエイリアス処理を行う
    y_raw = resample(x, p, q);
  end

  % --- 出力長を Nout に揃える ---
  M = length(y_raw);

  if M >= Nout
    % 長すぎる場合は末尾を切り捨て
    y = y_raw(1:Nout);
  else
    % 短すぎる場合は末尾をゼロパディング
    y = [y_raw, zeros(1, Nout - M)];
  end

  % 元が列ベクトルなら列ベクトルに戻す
  if is_col
    y = y(:);
  end

end



