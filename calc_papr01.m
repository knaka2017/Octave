tic;
N = 64;
M = 16;
sb = randint(N/2-1, 1, [0, M-1]);
s = qammod(sb, M);
x = [0; s; 0; conj(flipud(s))];

bufY = [0];
err = [];
maxx = [0 0 0 0]';
minx = [100 -100 0 0]';
count = 1;
for a=-100:1:100
  for b=-10:0.1:10
    x(1+N/2) = a + j*b;
    temp = real(ifft(x)); tempx = temp;
    y = max(temp) - min(temp);
    if (y > max(bufY) || count == 1)
      maxx = x;
    end
    if (y < min(bufY) || count == 1)
      minx = x;
    end
    bufY(count) = y;
    count++;
    temp = fft(temp);
    err = [err mean(abs(temp(2:(-1+N/2))-x(2:(-1+N/2))))];
  end
  plot(tempx); title([a y min(bufY)]); ylim([-2 2]);
  pause(0.1);
end

disp([maxx, real(ifft(maxx))]);
disp([minx, real(ifft(minx))]);

toc;
