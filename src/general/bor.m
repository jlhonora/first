N = 256;
M = 2;
Nacq = 4;
if 1==1,
	x = randn([1 M*N*2+2*N]);
	kx = randn([Nacq N]);
	kf = randn([Nacq N]);
	MM = complex(randn([Nacq N]),randn([Nacq N]));
	fo = randn([1 M]);
	xpos = (1:N)-N/2;
    phix = randn([1 N])+1;
else
	mw = zeros([1 N]); mf = mw;
	mw(1:8) = [1 2 3 4 3 2 1 0];
	mf(1:8) = [1 1 1 1 1 1 1 1];
	x = [real(mw) imag(mw) real(mf) imag(mf) zeros([1 N]) zeros([1 N])];
end;

%xpos = -0.4:0.1:0.3;
% kx = [xpos; xpos+10; xpos+20; xpos+30];

%MM = fftshift(fft(fftshift(mw)));
%MM = complex(real(MM),imag(MM));

%kx = -0.8:0.2:0.6;

RR = 2;
tc = zeros([1 RR+1]);
tm = zeros([1 RR+1]);
TT = 2;
for rr = 1:RR,
	tic
	for ii=1:TT,
		%[fc,gc] = gfR2b(x,kx,kf,MM,fo,xpos,phix);
        %[fc,gc] = gf5(x,kx,kf,MM,fo,xpos);
        [fc,gc] = gf7(x,kx,kf,MM,fo,xpos);
	end;
	tc(rr) = toc;
end;

for rr = 1:RR,
	tic
	for ii = 1:TT,
		%[fm,gm] = GradientFunctionR2v2(x,kx,kf,MM,fo,xpos,phix.');
        %[fm,gm] = GradientFunction5(x,kx,kf,MM,fo,xpos);
        [fm,gm] = GradientFunction7v2(x,kx,kf,MM,fo,xpos);
	end;
	tm(rr) = toc;
end;

fprintf('fc = %.1f fm = %.1f (diff = %.1f)\n',fc,fm,fc-fm);
fprintf('gc(10) = %f gm(10) = %f\n',gc(10),gm(10));

tc(RR+1) = mean(tc(1:RR))
tm(RR+1) = mean(tm(1:RR))

ratios = tc./tm*100
fprintf(sprintf('Average speedup: x%.4f\n',1/mean(ratios)*100))
