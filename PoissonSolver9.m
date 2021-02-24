function u = PoissonSolver9(m)
    f = @(x,y) x.^2+y.^2;


    h = 1/(m+1);

    x= linspace(1/(m+1),1-1/(m+1),m);
    y= linspace(1/(m+1),1-1/(m+1),m);

    T1 = diag(-20*ones(1,m)) + diag(4*ones(1,m-1),1) + diag(4*ones(1,m-1),-1);
    T2 = diag(4*ones(1,m)) + diag(1*ones(1,m-1),1) + diag(1*ones(1,m-1),-1);

    A = 1/(6*h^2)*full(blktridiag(T1,T2,T2,m)) ;

    F = [];
    for i = 1:m

        col = f(x(i),y);
        col(1) = col(1)-1/h^2;
        col(m) = col(m)-1/h^2;

        if i == 1 || i == m
            col(2:m-1) = col(2:m-1)-1/h^2;
            col(1) = col(1)-5/(6*h^2);
            col(m) = col(m)-5/(6*h^2);
            F = [F, col];
        else
            F = [F, col];
        end

    end

    F = F+h^2/3;

    u = reshape(A\F',[m,m]);
    
    % Add ones to the solution
    u = [ones(m,1),u];
    u = [u,ones(m,1)];
    
    u = u';
    u = [ones(m+2,1),u];
    u = [u,ones(m+2,1)];
    
 
end