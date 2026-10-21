function x = qamdemodM(sh, Ms)

x = [];
[kx, ky] = size(sh);
for k=1:kx
	x = [x; qamdemod(sh(k,:), Ms(k))];
end;
