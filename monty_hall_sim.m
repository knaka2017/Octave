M = 100;
N = 100000;

Nwin_n = zeros(1,N);
Nwin_s = zeros(1,N);
for k = 1:N
  car = randi(M);
  choice = 1;

  if (car == choice)
    Nwin_n(k) = 1;
  end

  hs = [];
  for d = 1:M
    if (d != choice && d != car)
      hs(end+1) = d;
    end
  end

  if (choice != car)
    hsc = car;
  else
    hsc = hs(randi(length(hs)));
  end

  if (car == hsc)
    Nwin_s(k) = 1;
  end

%  disp(sprintf("%d %f %f", k, sum(Nwin_n)/k, sum(Nwin_s)/k));
end

disp(sprintf("%f %f", sum(Nwin_n)/N, sum(Nwin_s)/N));
