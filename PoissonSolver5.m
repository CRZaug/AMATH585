function u = PoissonSolver5(m)

    h = 1/(m+1);

    x= linspace(1/(m+1),1-1/(m+1),m);
    y= linspace(1/(m+1),1-1/(m+1),m);
    
    f = @(x,y) x.^2+y.^2;

    T = diag(-4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);

    A = 1/h^2*full(blktridiag(T,eye(m),eye(m),m)) ;

    
    F = [];
    for i = 1:m

        col = f(x(i),y);
        col(1) = col(1)-1/h^2;
        col(m) = col(m)-1/h^2;

        if i == 1 || i == m
            F = [F, col-1/h^2];
        else
            F = [F, col];
        end


    end

    u = reshape(A\F',[m,m]);
    
    % Add ones to the solution
    u = [ones(m,1),u];
    u = [u,ones(m,1)];
    
    u = u';
    u = [ones(m+2,1),u];
    u = [u,ones(m+2,1)];
    
 
end