function res = estimationForTheta(H1, H2, q, theta, T)
%
%   DON'T FORGET TO RECOMPILE kappa0.cpp
%   FOR NEW VALUES OF H1 and H2
%
sz = 2^q-2;
constH1 = sqrt(gamma(3/2-H1) / (2*H1*gamma(H1+1/2)*gamma(3-2*H1)));

disp('start solving equation')
sol = lin_interp(T,H1,H2);
soln = sol.x';
t = sol.s;
disp('end solving equation')
N = size(soln,2);

disp('start calculating estimation')
res = 0;
numberOfPaths = 1000;
denom = zeros(1,numberOfPaths);
numer = zeros(1,numberOfPaths);
for j = 1:numberOfPaths
    Xt = GetXt(H1, H2, q, theta, T)';
    tx = linspace(0,T,sz)';
    Xti = interp1q(tx,Xt,t);
    denominator = 0;
    for k = 2:N
        denominator = denominator + ( soln(k) ) * ( (t(k))^(2-2*H1) - (t(k-1))^(2-2*H1) );
    end
    denominator = constH1 * denominator;
    denom(j) = denominator;
    numerator = 0;

    for m = 2:N 
        numerator = numerator + soln(m-1) * ( Xti(m) - Xti(m-1) );
    end
    numer(j) = numerator;
    est = numerator/denominator;
    res = res + est;
end

res = res/numberOfPaths;
end
