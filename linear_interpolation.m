function [t,soln] = linear_interpolation(Time, H1, H2)
%
%   DON'T FORGET TO RECOMPILE kappa0.cpp
%   FOR NEW VALUES OF H1 and H2
%
atol = 1/10^6;
rtol = 1/10^3;
alpha_ = 1 + 2*H1 - 2*H2;
beta_ = 2*H1 - 1;
gammaH1 = sqrt( (gamma(3/2-H1)) / (2*H1*gamma(H1+1/2)*gamma(2-2*H1)) );
%deltaH1 = (2 - 2*H1) * beta(3/2 - H1,3/2 - H1) / gammaH1;
lower_bound = 0;
upper_bound = Time;
n = 2^3;
[told,solnold] = mod_simp(n);
disp('outside the main loop');
max_m = 9;
for m = 4:max_m
    n = 2^m;
    disp(n);
    [t,soln] = mod_simp(n); 
    temp = soln(1:2:n+1);
    delta = norm(temp - solnold, inf);
    if m > 4
        rate = min(0.5, max(delta/deltaold, 0.0625));
        errest = (rate /(1-rate)) * delta;
        
        disp(errest);
        disp(atol + rtol*norm(soln,inf));
    end
        
    deltaold = delta;
    solnold = soln;
end
    
%===Nested functions=======================================================

function [t,soln] = mod_simp(n)
    h = (upper_bound - lower_bound)/n;
    h_power1 = h^(1-alpha_);
    h_power2 = h^(1-alpha_ - beta_);

    t = linspace(lower_bound,upper_bound,n+1)';
    [S,T] = meshgrid(t,t);
    disp('start computing kappa0');
    Kvals = kappa0(S,T)';
    disp('end computing kappa0');
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