function ai = alg_integrals(b,alpha)
% For k = 0,1, evaluate the integrals
%                     1
%       I_k = Integral  u^k / |u+c|^alpha du
%                     0 
% The exponent alpha satisfies 0 < alpha < 1. 

I_0 = ( 1 / (1 - alpha) ) * ( (1+b)^(1-alpha) - b^(1-alpha) );
I_1 = ( 1 / (2 - 3*alpha + alpha^2) ) * ( b^(2-alpha) - (1 + b)^(1-alpha) * (alpha + b - 1) );

ai = [I_0,I_1];

end % alg_integrals
