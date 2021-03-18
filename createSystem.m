function [Aj,fj] = createSystem(n,m)


h = 1/(m+1);
l = 0;
r = 1;
ul = 0;
ur = 0;

x = linspace(l+r/(m+1),r-r/(m+1),m);

Aj = 1/h^2*(diag(ones(m-1,1),-1)-2*diag(ones(m,1))+diag(ones(m-1,1),1));

f = -(6*x-2);
fj = f';



end