function res = GetXt (H1, H2, q, theta, T)
N = 2^q - 1;
x = 0:T/(N-2):T;
fBmH1 = GetfBm(H1, q, T);
fBmH2 = GetfBm(H2, q, T);
x = x.^(2-2*H1);
drift = theta * beta(3/2-H1,3/2-H1) * x;
res = drift + fBmH1 + fBmH2;
end

