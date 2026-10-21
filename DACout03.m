function [uuout uusout uusout_impulse LSB] = DACout03(uu, Nq, fs, fdac, fl_shape)
% uu: signal (floating)
% Nq: quantiz level
% fs: sampling freq. [Hz]
% fdac: DAC speed [Sa/s]
% fl_shape: 1:rect, 2:sinc, otherwise:rect

% Sf = 0.9;		% Scale factor
Sf = 1;

uumin = min([min(real(uu)), min(imag(uu))]);
uumax = max([max(real(uu)), max(imag(uu))]);
LSB = (uumax-uumin)/(2^Nq-1);

uu = uu * Sf;

[uuindex_r, uu_r] = quantiz(real(uu)-(LSB/2), linspace(uumin,uumax,2^Nq), linspace(uumin,uumax,2^Nq));
[uuindex_i, uu_i] = quantiz(imag(uu)-(LSB/2), linspace(uumin,uumax,2^Nq), linspace(uumin,uumax,2^Nq));
uuout = complex(uu_r, uu_i);
%uu = complex(uuindex_r, uuindex_i);
%uuout = complex(uuindex_r, uuindex_i);


% インパルス列
uusout_impulse = upsample(uuout, floor(fs/fdac));

%%% 出力パルス整形フィルタ
switch fl_shape
	case 1
		% 矩形フィルタ
		b = ones(floor(fs/fdac),1);  a = 1;
	case 2
		% sinc
		b = sinc(-4:0.25:4);  a = 1;
	otherwise
		% 矩形フィルタ
		b = ones(floor(fs/fdac),1);  a = 1;
endswitch
uusout = filter(b,a, uusout_impulse); % 出力パルス整形フィルタ

%uus = resample(uuout, fs, fdac);

