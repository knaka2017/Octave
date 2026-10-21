function x = qammodM(sh, Ms)

x = [];
[kx, ky] = size(sh);
for k=1:kx
	x = [x; qammod(sh(k,:), Ms(k))];
end;
