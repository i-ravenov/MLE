function res = GetfBm (H, q, T)
N = 2^q - 1;
lambda = Lambda(H,N);
fGnsample = FGN(lambda);
res = cumsum(fGnsample);
res = res * (T/N)^H; 
%x = 0:T/(N-2):T;
%plot(x,res);