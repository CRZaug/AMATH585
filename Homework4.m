%% Homework 4

%% Problem 1
%
% We use the `chebfun` package to solve a second order BVP with homogeneous
% boundary conditions.

%%
d = domain(0,1);
L = chebop(@(x,u) -1*diff((1+x^2).*diff(u,1),1),d,0,0);
f = chebfun(@(x) 2*(3*x.^2-x+1),d);
u = chebfun(@(x) x.*(1-x),d);

approx = L\f;

infNorm = norm(u-approx,inf)
L2Norm = norm(u-approx)

figure(1)
plot(approx)

%%
% The infinity-norm and the L2-norm are both near machine precision
% (very small).

%% Problem 2
% 
% We solve Poisson's equation in the unit square with nonhomogeneous
% boundary conditions.
%%

m = 99;

uexact = PoissonSolver5(m);

[X_fine,Y_fine] = meshgrid([0,linspace(1/(m+1),1-1/(m+1),m),1]);

figure(2)
surf(uexact)


%%
error = [];
hVector = [];

for n = 1:4
    
    m = 10*n-1;
    u = PoissonSolver5(m);
       
    [X_coarse,Y_coarse] = meshgrid([0,linspace(1/(m+1),1-1/(m+1),m),1]);
    
    unew = interp2(X_coarse, Y_coarse, u, X_fine, Y_fine);

    e = norm(unew-uexact);

    error = [error,e];
    hVector = [hVector, 1/(m+1)];
    
end

T = table(hVector',error');
T.Properties.VariableNames = {'h' 'Error'};
T

%%
% From the table, the error does not appear to have second-order convergence.
% This is likely due to an error comparing values across grids.


%% Problem 3
% We repeat the last problem with a 9-point stencil for the Laplacian.

m = 99;
uexact = PoissonSolver9(m);

[X_fine,Y_fine] = meshgrid([0,linspace(1/(m+1),1-1/(m+1),m),1]);

figure(3)
surf(uexact)


%%

error = [];
hVector = [];
for n = 1:5
    
    m = 10*n-1;
    u = PoissonSolver9(m);
    
    [X_coarse,Y_coarse] = meshgrid([0,linspace(1/(m+1),1-1/(m+1),m),1]);
    
    unew = interp2(X_coarse, Y_coarse, u, X_fine, Y_fine);
        
    e = norm(unew-uexact);
    
    error = [error,e];
    hVector = [hVector, 1/(m+1)];
    
end

T = table(hVector',error');
T.Properties.VariableNames = {'h' 'Error'};
T
%%
% From the table, we see that the error still doesn't have the expected convergence (fourth-order).
% This is again due to probable issues comparing values between different
% grid sizes. 
%
% The reason (in general) that we might expect not to see 4th order convergence mentioned in
% the text is roundoff error. The errors here (even though not converging rapidly) are small, and as the grid
% size decreases, at fourth-order, the error would rapidly near machine precision at fourth order.

