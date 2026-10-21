function sb = randintM(N, BLEN, Ms, SEEDs)

sb = [];
for k=1:N
	sb = [sb; randint(1,BLEN, [0, Ms(k)-1], SEEDs+k)];
end;
