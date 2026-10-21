function y = extract_one_period(x)
  % 周期的な信号波形ベクトルxから、１周期分yを取り出す
  % xには白色雑音も含まれる
  % xは縦ベクトル、横ベクトルいずれの場合もありえる
  % yの波形値およびyの長さ（1周期分）いずれも未知
  % xの長さ > yの長さ

  % 入力ベクトルが横ベクトルの場合は縦ベクトルに変換
  if (size(x, 2) > size(x, 1))
    x = x';
  endif

  % 1. xの自己相関関数を計算
  correlation = xcorr(x);

  % 2. 自己相関関数のピークを検出
  % 中央のピーク（ラグ0）は除外するため、中央付近のインデックスを計算
  n = length(correlation);
  center_index = floor(n / 2) + 1;

  % ピーク検出の範囲を設定（中央のピーク付近は除く）
  search_start = 1;
  search_end = n;

  % 極大値を探索する（簡易的なピーク検出）
  peaks_indices = [];
  for i = 2:(n - 1)
    if (correlation(i) > correlation(i - 1) && correlation(i) > correlation(i + 1))
      if (i < center_index - 5 || i > center_index + 5) % 中央のピーク付近は除く
        peaks_indices = [peaks_indices, i];
      endif
    endif
  end

  if (isempty(peaks_indices) || length(peaks_indices) < 2)
    warning("自己相関関数で十分な数のピークが検出できませんでした。");
    y = [];
    return;
  endif

  % 3. 隣り合うピークの時間間隔を計算し、周期の候補とする
  period_candidates = diff(sort(peaks_indices));

  if (isempty(period_candidates))
    warning("隣り合うピークが見つかりませんでした。");
    y = [];
    return;
  endif

  % 最も頻繁に出現する時間間隔を周期の長さの推定値とする
  period_length_estimate = mode(period_candidates);

  if (isnan(period_length_estimate) || period_length_estimate <= 0 || period_length_estimate >= length(x))
    warning("推定された周期長が不正です。");
    y = [];
    return;
  endif

  % 4. xから1周期分の波形を抽出
  % 自己相関関数の最初の有意なピークの位置を探す（中央のピーク以降）
  first_peak_after_center_index = -1;
  for i = center_index + 1:n
    if (correlation(i) > 0) % 自己相関が正の値になった最初の点
      first_peak_after_center_index = i;
      break;
    endif
  end

  if (first_peak_after_center_index == -1)
    warning("自己相関関数の有意なピークが見つかりませんでした。");
    y = [];
    return;
  endif

  % 推定された周期長でxから1周期分の波形を抽出
  % 開始位置は、自己相関の最初の有意なピークに対応するラグの後にする
  start_index = first_peak_after_center_index - center_index + 1;

  % 開始位置がxの範囲外の場合は調整
  if (start_index <= 0)
    start_index = 1;
  endif

  % 抽出する終端位置を計算
  end_index = start_index + period_length_estimate - 1;

  % 終端位置がxの範囲を超える場合は調整
  if (end_index > length(x))
    end_index = length(x);
    if (end_index - start_index + 1 < period_length_estimate / 2) % 極端に短い場合は抽出を諦める
      warning("抽出範囲が短すぎます。");
      y = [];
      return;
    endif
  endif

  % 波形を抽出
  y = x(start_index:end_index);

  % 結果が横ベクトルの場合は縦ベクトルに戻す
  if (size(y, 2) > size(y, 1))
    y = y';
  endif

endfunction
