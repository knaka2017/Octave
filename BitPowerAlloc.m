function [zb, zp] = BitPowerAlloc(x, y, b, p)
%% Bit and Power loading Allocation
%%
%% [zb, zp] = BitPowerAlloc(x, y, b, p)
%% zb : Bit loading (2^zb-QAM)
%% zp : Power loading [dB]
%%
%% x : sub-carrier index
%% y : SNR [dB] @ x
%% b : set of Bit loading coef (2^b-QAM)
%% p : set of threshold Power loading coef [dB] (same size of b)
%%
%% ex.) [zb, zp] = BitPowerAlloc(fig8x, fig8y, [8, 6, 4, 2], [28, 22, 16, 10]);

  % 縦ベクトルに変換
  [lx, kx] = size(x);
  if (lx <= kx)
    x = x.';
  end
  [ly, ky] = size(y);
  if (ly <= ky)
    y = y.';
  end;
  % 直線近似
  X = [ones(size(x)), x];
  A = X\y;
  Y = X*A;
  z = Y;
  % b, pのソート(大→小)
  [B, iB] = sort(b, 'descend');
  P = p(iB);
  % 階段の各ステップの間隔（変化分）
  dP = [inf, abs(diff(P)), 6];
  % パワー割り当て＆ビット割り当て
  zb = [];
  zp = [];
  for k=1:length(b)
    temp = z < (P(k) + (dP(k)/2)) & z >= (P(k)-(dP(k+1)/2));
    zb(temp) = B(k);
    zp(temp) = P(k);
  end;
  temp  = z < (P(k)-(dP(k+1)/2));
  zb(temp) = 0;
  zp(temp) = -inf;
  if (ly <= ky)
    Y = Y.';
  end;
  zp = zp - Y;
  % 確認用グラフ
  subplot(3,1,1);
    plot(x, y, '.', x, Y, 'b-', x, zp+Y, 'r-');
    xlabel("sub-carrier index");  ylabel("SNR [dB]");
  subplot(3,1,2);
    plot(x, zb, 'o');
    xlabel("sub-carrier index");  ylabel("number of bit");
  subplot(3,1,3);
    plot(x, zp, '-');
    xlabel("sub-carrier index");  ylabel("Power coef");
  save("-mat", "zbzp.mat", "zb", "zp");
endfunction
