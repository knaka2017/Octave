function y = extract_periodic_signal(x)
    % 周期的な信号波形ベクトルxから1周期分を取り出す関数
    % 入力:
    %   x - 測定された信号ベクトル（縦または横ベクトル）
    % 出力:
    %   y - 抽出された1周期分の信号ベクトル

    % ベクトルを縦ベクトルに変換
    x = x(:);

    % 自己相関を計算
    [autocorr, lags] = xcorr(x, 'coeff');

    % 正のラグ部分に絞る（最初の自己相関のピークを見つける）
    pos_autocorr = abs(autocorr(length(x):end));

    % 最初のピークを検出（0.9以上の相関係数のピーク）
    [pks, peak_indices] = findpeaks(pos_autocorr, 'MinPeakHeight', 0.3);

    % 最初のピークが見つからない場合のフォールバック
    if isempty(peak_indices)
        % ピークの低い閾値で再試行
        [pks, peak_indices] = findpeaks(pos_autocorr, 'MinPeakHeight', 0.2);
    end

    % それでもピークが見つからない場合
    if isempty(peak_indices)
        error('周期を検出できませんでした');
    end

    % 最初のピークの位置が周期
    period = peak_indices(1);

    % 入力信号が周期より短い場合はエラー
    if length(x) < period
        error('入力ベクトルが1周期より短いです');
    end

    % 最適な開始点を見つける
    best_start = 1;
    best_total_correlation = -Inf;

    % 可能な開始点を探索
    for start = 1:(length(x) - period + 1)
        current_cycle = x(start:start+period-1);
        total_correlation = 0;
        valid_comparisons = 0;

        % 可能な限り多くの周期と比較
        for j = 1 : floor((length(x) - start) / period)
            comparison_start = start + j * period;

            % サイクルの長さが一致することを確認
            if comparison_start + period - 1 <= length(x)
                comparison_cycle = x(comparison_start : comparison_start + period - 1);

                % 相関を計算
                try
                    correlation = corr(current_cycle, comparison_cycle);
                    total_correlation = total_correlation + correlation;
                    valid_comparisons = valid_comparisons + 1;
                catch
                    % 何らかの理由で相関計算に失敗した場合はスキップ
                    continue;
                end
            end
        end

        % 平均相関を計算
        if valid_comparisons > 0
            avg_correlation = total_correlation / valid_comparisons;

            % より良い相関が見つかった場合に更新
            if avg_correlation > best_total_correlation
                best_total_correlation = avg_correlation;
                best_start = start;
            end
        end
    end

    % 1周期分のベクトルを抽出
    y = x(best_start : best_start + period - 1);
end

