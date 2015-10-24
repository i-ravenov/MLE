function xint = ntrpFie(sol,sint)
% Evaluate the Nystroem interpolant for a solution structure sol at the
% points SINT.  
    iscol = (size(sint,2) == 1);
    if iscol, sint = sint'; end
    t = sol.s; x = sol.x;            % local variables
    rhs = ones(size(sint));
    [S,T] = meshgrid(sint,t);
    kermat = kappa0(S,T)';
    n = length(t) - 1;   % The number of subdivisions.

    xint = zeros(size(sint)); 
    alpha = sol.alpha;            
    for j = 1:n
        h = t(j) - t(j-1);
        h_power = 0.5*h^(1-alpha);
            for i = 1:length(sint)
                    tau = sint(i);
                    
                    beta = (tau - t(j))/h;
                    ints = alg_integrals(beta,alpha);
                    I1 =   h_power*(ints(3) - ints(2));
                    I2 = 2*h_power*(ints(1) - ints(3));
                    xint(i) = xint(i) + I1*x(j-1)*kermat(i,j-1)...
                           + I2*x(j)*kermat(i,j) + I3*x(j+1)*kermat(i,j+1);   
            end
    end
    xint = (rhs - xint);
end % ntrpFie
