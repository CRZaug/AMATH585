function  [iterationVect, residualVect,x] = SOR(A,b,x,tol,omega,iter,pri)



M = 1/omega*diag(diag(A))+tril(A,-1);


x_new = ones(length(x),1);

residualVect =[];

for i = 1:iter
    
    r = b-A*x;
    z = M\r;
    x_new = x+z;
    
    e = sum(abs(x-x_new));
    
     if e < tol
         if pri == "true"
            fprintf("SOR converged within tolerance %f after %i iterations\n", tol, i)
         end
        break
     end
    
    
    residualVect = [residualVect,norm(b-A*x)/norm(b)];
     
    x = x_new;
    
end

iterationVect = 1:(i-1);

end
