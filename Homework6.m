%% Homework 6

%% Problem 2
% We repeat the experiment in the text to see if GS and CG make good
% smoothers for a multigrid method. We first create the described system
% below.

%% 

n = 255;

h = 1/(n+1);
l = 0;
r=1;
a = 0.5;
ul = 1;
ur = 3;

phi = @(x) 20 * pi * x.^3;
phi1 =  @(x) 60 * pi * x.^2;
phi2 = @(x) 120 * pi * x;

actual = @(x) 1+12*x-10*x.^2+a* sin(phi(x));

x = linspace(l+r/(n+1),r-r/(n+1),n);

A = 1/h^2*spdiags(ones(n,1)*[-1 2 -1],-1:1,n,n);

f = -1*(-20 + a*phi2(x).*cos(phi(x))-a*phi1(x).^2.*sin(phi(x)));
f(n) = f(n)+(1/h^2)*ur;
f(1) = f(1) +1/h^2*ul;
f = f';


u0 = (1+2*x)';
%%
% First we solve with GS.
%%
GS5 = GaussSeidelH6(A,f,u0,5);
GS10 = GaussSeidelH6(A,f,u0,10);
GS20 = GaussSeidelH6(A,f,u0,10000);


figure(1)
subplot(4,2,1);
hold on;
plot(x,u0)
plot(x,actual(x))
title("GS after 0 steps")

subplot(4,2,2);
hold on;
plot(x,u0'-actual(x))
title("Error after 0 steps")

subplot(4,2,3);
hold on;
plot(x,GS5)
plot(x,actual(x))
title("GS after 5 steps")

subplot(4,2,4);
hold on;
plot(x,GS5'-actual(x))
title("Error after 5 steps")

subplot(4,2,5);
hold on;
plot(x,GS10)
plot(x,actual(x))
title("GS after 10 steps")

subplot(4,2,6);
hold on;
plot(x,GS10'-actual(x))
title("Error after 10 steps")

subplot(4,2,7);
hold on;
plot(x,GS20)
plot(x,actual(x))
title("GS after 20 steps")

subplot(4,2,8);
hold on;
plot(x,GS20'-actual(x))
title("Error after 0 steps")

%%
% Next, we solve with CG (preconditioned with the incomplete Cholesky
% decomposition).
%%
L = ichol(A);

[CG0,~,~,~] = pcg(A,ones(255,1),1e-10,0,L,L',u0);
[CG5,~,~,~] = pcg(A,f,1e-10,5,L,L',u0);
[CG10,~,~,~] = pcg(A,f,1e-10,10,L,L',u0);
[CG20,~,~,~] = pcg(A,f,1e-10,20,L,L',u0);

figure(2)
subplot(4,2,1);
hold on;
plot(x,CG0)
plot(x,actual(x))
title("CG after 0 steps")

subplot(4,2,2);
hold on;
plot(x,CG0'-actual(x))
title("Error after 0 steps")

subplot(4,2,3);
hold on;
plot(x,CG5)
plot(x,actual(x))
title("CG after 5 steps")

subplot(4,2,4);
hold on;
plot(x,CG5'-actual(x))
title("Error after 5 steps")

subplot(4,2,5);
hold on;
plot(x,CG10)
plot(x,actual(x))
title("CG after 10 steps")

subplot(4,2,6);
hold on;
plot(x,CG10'-actual(x))
title("Error after 10 steps")

subplot(4,2,7);
hold on;
plot(x,CG20)
plot(x,actual(x))
title("CG after 20 steps")

subplot(4,2,8);
hold on;
plot(x,CG20'-actual(x))
title("Error after 20 steps")

%%
% From the plots, we can see that CG reduces the error (generally) much more rapidly
% than GS. However, it does a bad job of reducing the high frequency
% components of the error. This means it would not be a good smoother. Even though GS
% converges more slowly overall, it reduces the high frequency error
% components and so would be an effective smoother.

%% Problem 3
%
% We implement a two grid method using Gauss-Seidel and Weighted Jacobi as
% smoothers. We then test their convergence as a function of different grid
% sizes. We test this out on the problem $u''(x) = -6x-2$, which satisfies
% the boundary conditions $u(0)=u(1)=0$. The actual solution is $u(x) =
% -x^3-x^2$.

%%  
% Gauss Seidel

n = 100;

actual = @(x) -(x.^3-x.^2);
x = linspace(0+1/(n+1),1-1/(n+1),n);


[u,iConvGS] = twoGridGS(n,10000,1e-5);
         
figure(3)
hold on;
plot(x,u)
plot(x,actual(x),'.')
legend("2grid","actual")
title("2 Grid Method with GS Smoother")
         
%% 
% Weighted Jacobi


n = 100;

actual = @(x) -(x.^3-x.^2);
x = linspace(0+1/(n+1),1-1/(n+1),n);

[u,iConvWJ] =twoGridWJ(n,30,1e-5);

figure(4)
hold on;
plot(x,u) 
plot(x,actual(x),'.')
legend("2grid","actual")
title("2 Grid Method with WJ Smoother")
 
%% 
% We see that both GS and WJ as smoothers cause the 2-grid method to
% converge rather quickly, within 20 iterations. Now we wish to show that
% this is independent of grid size.

%%

GSvect = [];
WJvect = [];

for j = 1:10
    n = 2^j;
    [~,iConvGS] = twoGridGS(n,10000,1e-5);
    [~,iConvWJ] =twoGridWJ(n,30,1e-5);
    
    GSvect  = [GSvect,iConvGS];
    WJvect  = [WJvect,iConvWJ];
end

figure(5)
plot(GSvect)
hold on;
plot(WJvect)
title("Convergence of the 2-grid method")
legend("GS", "WJ")
xlabel("Grid size log(n)")
ylabel("Iterations to converge")

%% 
% From the plots, the size of the grid has little effect on the convergence
% of the 2-grid method. Even with a grid size of $2^{10}$, the method takes
% less than iterations to converge.