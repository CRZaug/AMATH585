function [iterationVect, residualVect,x] = GaussSeidel(A,b,x,tol,iter)

M = tril(A);

x_new = ones(length(x),1);

residualVect = [];

for i = 1:iter
    
    r = b-A*x;

    z = M\r;
    x_new = x+z;
    
    e = sum(abs(x-x_new));

    
    if e < tol
        fprintf("Gauss Seidel converged within tolerance %f after %i iterations\n", tol, i)
        break
    end
    
    residualVect = [residualVect,norm(b-A*x)/norm(b)];

    x = x_new;
    
end

iterationVect = 1:(i-1);





end
