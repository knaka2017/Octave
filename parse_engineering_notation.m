% ========================================================================
% 補助関数: エンジニアリング表記を解析
% ========================================================================
function value = parse_engineering_notation(str)
% PARSE_ENGINEERING_NOTATION エンジニアリング表記を数値に変換
%
% 入力例: '16G' -> 16e9, '800M' -> 800e6, '1.5K' -> 1500
%
% 対応プレフィックス: K (kilo), M (mega), G (giga), T (tera)

    str = upper(strtrim(str));

    prefix_mult = 1;
    last_char = str(end);

    % プレフィックスをチェック
    switch last_char
        case 'K'
            prefix_mult = 1e3;
            numeric_str = str(1:end-1);
        case 'M'
            prefix_mult = 1e6;
            numeric_str = str(1:end-1);
        case 'G'
            prefix_mult = 1e9;
            numeric_str = str(1:end-1);
        case 'T'
            prefix_mult = 1e12;
            numeric_str = str(1:end-1);
        otherwise
            numeric_str = str; % プレフィックスなし
    endswitch

    numeric_part = str2double(numeric_str);

    % str2doubleの失敗をisnanでチェック
    if isnan(numeric_part)
        value = 0; % 変換失敗時はゼロを返す
        warning('parse_engineering_notation:ConversionFailed', ...
                '文字列 "%s" を数値に変換できませんでした。', str);
    else
        value = numeric_part * prefix_mult;
    endif
endfunction
