function [sol] = lin_interp(Time, H1, H2)
%
%   DON'T FORGET TO RECOMPILE kappa0.cpp
%   FOR NEW VALUES OF H1 and H2
%
alpha_ = 1 + 2*H1 - 2*H2;
beta_ = 2*H1 - 1;
gammaH1 = sqrt( 2 * H1 * (2-2*H1) * ((gamma(3/2-H1))^3) * gamma(H1+1/2) / (gamma(3-2*H1)) );
lower_bound = 0;
upper_bound = Time;

m = 8;
n = 2^m;
disp(n);
[t,soln] = mod_simp(n);
sol.s = t;
sol.x = soln;
sol.H1 = H1;
sol.H2 = H2;
sol.T = Time;

%===Nested function========================================================

function [t,soln] = mod_simp(n)
    h = (upper_bound - lower_bound)/n;
    h_power1 = h^(1-alpha_);
    h_power2 = h^(1-alpha_ - beta_);

    t = linspace(lower_bound,upper_bound,n+1)';
    [S,T] = meshgrid(t,t);
    Kvals = kappa0(S,T);
    Msize = n + 1;
    rhs = ones(Msize, 1);
    kermat = zeros(Msize);           
    
    for j = 1:n
        for i = 1:n+1
            if ( t(j) >= t(i) )
                b = (t(j) - t(i))/h;
                ints = alg_integrals(b,alpha_);
                I1 = h_power1 * (ints(1) - ints(2)); % right
                I2 = h_power1 * ints(2); 
                kermat(i,j) = kermat(i,j) + I1;
                kermat(i,j+1) = kermat(i,j+1) + I2;
            else 
                a = t(j) / h;
                b = (t(i) - t(j)) / h;
                ints = alg_integrals_doublep(a,b,alpha_,beta_);
                I1 = h_power2 * (t(i))^beta_ * (ints(1) - ints(2));
                I2 = h_power2 * (t(i))^beta_ * ints(2);
                
                kermat(i,j) = kermat(i,j) + I1;
                kermat(i,j+1) = kermat(i,j+1) + I2;                
            end
        end
    end
    kermat = 1/(gammaH1^2)*kermat;
    kermat = Kvals.*kermat;    
    for j = 1:Msize
        kermat(j,j) = kermat(j,j) + 1;
    end
    soln = kermat \ (rhs);  
end % mod_simp

%==========================================================================
end 