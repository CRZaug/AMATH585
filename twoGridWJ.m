function [u,k] = twoGridWJ(n,iter,tol)


actual = @(x) -(x.^3-x.^2);
x = linspace(0+1/(n+1),1-1/(n+1),n);

u_act = actual(x);

[A,f] = createSystem(n,n);
u=zeros(n,1);

hold on;

    % Create the projection matrix
    if mod(n,2) ==0
        m = n/2; 
    else
        m = floor(n/2)+1;
    end
 

    % Create the proj/interp matrices
    Ip = zeros(m,n);
    for i = 1:m-1
            Ip(i,2*i-1) = 0.5;
            Ip(i,2*i) = 1;
            Ip(i,2*i+1) = 0.5;
    end

    Ip(m,n-1)=0.5;
    Ip(m,n)= 1;
    Ii = 0.5*Ip';

   
    % Begin the cycle
     for k = 1:iter
         
         r = f-A*u;
         
         % Project residual onto grid level
         fj = Ip*r;
         [Aj,~] = createSystem(n,m);
         
         % Directly solve on the coarsest grid
         d = Aj\fj;
        
         
         % Interpolate
         unew = u+Ii*d;
        
         % Perform WJ
         
         M = 3/2*diag(diag(A));
         unew = unew+M\(f-A*unew);
         
         e = sum(abs(u-unew));
         
         if e < tol
            fprintf("2grid converged within tolerance %f after %i iterations\n", tol, k)
            break
         end

        u = unew;
        
     end
end