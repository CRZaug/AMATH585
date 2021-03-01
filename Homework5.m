%% AMATH 585 Homework 5
%
% Camille Zaug
%
% 3/1/2021

%% Problem 2
% We use different iterative methods to solve a two-dimensional boundary
% problem (discretized using finite differences). We study the convergence
% of these methods as a function of the number of iterations.

%%
% Start with 10 meshpoints and create a matrix.

n=10;
h = 1/n;
N = (n-1)^2;
[A,b] = createMatrix(n);


%% 
% 
% The SOR method has a parameter $\omega$ that should be optimized. Here we
% try a few different values and plot the number of iterations it takes for
% the method to converge.
%

%%

testspace = linspace(1,2,100);
iterations = [];
for o = 1:length(testspace)
    omega = testspace(o);
    [iters_SOR, res_SOR, u_SOR] = SOR(A,b,ones(N,1),1e-6,omega,100000,"false");
    iterations = [iterations, iters_SOR(end)];
end

semilogy(testspace,iterations)
title("Finding an optimal \omega for SOR")
ylabel("Number of iterations until convergence")
xlabel("\omega")

%% 
% In the plot, the optimal value of $\omega$ occurs around 1.5. However,
% this changes depending on the number of grid points. In the end, we
% settle on $$ \omega =2-2\pi h$ (where $h$ is the mesh width) as our
% optimal value. This value is smaller than 2 and the same as the value
% described for the Poisson problem in the text.
%
% Next, we test the Jacobi, Gauss-Siedel, and SOR methods on this matrix.
% We plot their convergence below.

%% 

omega_opt =  2-2*pi*h;
[iters_J, res_J, u_Jacobi] = Jacobi(A,b,ones(N,1),1e-6,100000);
[iters_GS, res_GS, u_GS] = GaussSeidel(A,b,ones(N,1),1e-6,10000);
[iters_SOR, res_SOR, u_SOR] = SOR(A,b,ones(N,1),1e-6,omega_opt,100000,"true");

figure(1)
loglog(iters_J, res_J)
hold on;
loglog(iters_GS, res_GS)
loglog(iters_SOR, res_SOR)
title("Convergence of Iterative Methods")
legend('Jacobi','Gauss-Siedel',"SOR")
xlabel("Iterations")
ylabel("Relative Residual Norm ")

%%
% From the plot, Jacobi converges the slowest, Gauss-Siedel converges
% faster than Jacobi, and SOR converges the fastest. This is predicted in
% the text. 
%
%%
% Now we solve the system with the conjugate gradient method with and without
% preconditioners (we use the incomplete Cholesky decomposition as a
% preconditioner).

%%

% CG without preconditioner
[u_CG,~,~,iter_CG,resvec_CG] = pcg(A,b,1e-6,300);

% Incomplete Cholesky decomposition
L = ichol(A);

% CG with preconditioner
[u_CGP,~,~,iter_CGP,resvec_CGP] = pcg(A,b,1e-6,300,L,L');

figure(2)
loglog(0:iter_CG,resvec_CG)
hold on;
loglog(0:iter_CGP,resvec_CGP)
title("Convergence of Iterative Methods")
legend('CG', 'PCG')
xlabel("Iterations")
ylabel("Relative Residual Norm ")
%%
% From the plot, using the incomplete Cholesky as a precondition 
% clearly causes faster convergence of the CG method. Without a
% preconditioner, CG converges about as fast as SOR.
%%
% Finally, we repeat with different grid sizes to see how the convergence
% changes.
%
%

%%
n = 50;
h = 1/n;
N = (n-1)^2;
[A,b] = createMatrix(n);

% Optimal omega
omega_opt =  2-2*pi*h;

% Incomplete Cholesky decomposition
L = ichol(A);

[iters_J, res_J, u_Jacobi] = Jacobi(A,b,ones(N,1),1e-6,100000);
[iters_GS, res_GS, u_GS] = GaussSeidel(A,b,ones(N,1),1e-6,10000);
[iters_SOR, res_SOR, u_SOR] = SOR(A,b,ones(N,1),1e-6,omega_opt,100000,"true");
[u_CG,~,~,iter_CG,resvec_CG] = pcg(A,b,1e-6,300);
[u_CGP,~,~,iter_CGP,resvec_CGP] = pcg(A,b,1e-6,300,L,L');

figure(3)
loglog(iters_J, res_J)
hold on;
loglog(iters_GS, res_GS)
loglog(iters_SOR, res_SOR)
loglog(0:iter_CG,resvec_CG)
loglog(0:iter_CGP,resvec_CGP)
title("Convergence of Iterative Methods")
legend('Jacobi','Gauss-Siedel',"SOR",'CG', 'PCG')
xlabel("Iterations")
ylabel("Relative Residual Norm ")

%% 
n = 100;
h = 1/n;
N = (n-1)^2;
[A,b] = createMatrix(n);

% Optimal omega
omega_opt =  2-2*pi*h;

% Incomplete Cholesky decomposition
L = ichol(A);

[iters_J, res_J, u_Jacobi] = Jacobi(A,b,ones(N,1),1e-6,100000);
[iters_GS, res_GS, u_GS] = GaussSeidel(A,b,ones(N,1),1e-6,100000);
[iters_SOR, res_SOR, u_SOR] = SOR(A,b,ones(N,1),1e-6,omega_opt,100000,"true");
[u_CG,~,~,iter_CG,resvec_CG] = pcg(A,b,1e-6,300);
[u_CGP,~,~,iter_CGP,resvec_CGP] = pcg(A,b,1e-6,300,L,L');

figure(4)
loglog(iters_J, res_J)
hold on;
loglog(iters_GS, res_GS)
loglog(iters_SOR, res_SOR)
loglog(0:iter_CG,resvec_CG)
loglog(0:iter_CGP,resvec_CGP)
title("Convergence of Iterative Methods")
legend('Jacobi','Gauss-Siedel',"SOR",'CG', 'PCG')
xlabel("Iterations")
ylabel("Relative Residual Norm ")

%% 
% From the plots, it appears that the relative behavior of each of the
% iterative methods is the same regardless of the number of mesh
% points. Preconditioned CG converges fastest and Jacobi converges slowest. However, the number of iterations to converge increases overall.
% When $n=50$, no method took longer than $10^4$ iterations to converge.
% When $n=100$, however, Jacobi and Gauss-Siedel took over $10^4$ iterations
% to converge. As adding more mesh points increases the size of the
% problem, it makes sense that convergence takes longer.