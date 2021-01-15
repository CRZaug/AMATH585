function u = heatEquation(a,b,h,func)
    % Solve the heat equation with conductivity 
    % a: The left endpoint
    % b: The right endpoint
    % h: Mesh size
    % f: RHS function
    
    ua = 1; % The value of u(a)

    % Conductivity
    c = @(x) (1+x^2);

    m = (b-a)/h; % m out of m+1 grid points

    % Initialize the matrix and f
    uMatrix = zeros(m);
    fVector = zeros(m,1);

    % Fill in the first row
    x1 = a+h;
    x1_m = x1-h/2;
    x1_p = x1+h/2;
    uMatrix(1,1) = -1/h^2*(c(x1_m)+c(x1_p));
    uMatrix(1,2) = 1/h^2*c(x1_p);

    fVector(1) = func(x1)-ua*(1/h^2*c(x1_m));


    % Fill in all rows except the last
 
    for j = 2:m-1

        xj = a+j*h;
        xj_m = xj-h/2;
        xj_p = xj+h/2;

        uMatrix(j,j-1) = c(xj_m)/h^2;
        uMatrix(j,j) = -1*(c(xj_m)+c(xj_p))/h^2;
        uMatrix(j,j+1) = c(xj_p)/h^2;

        fVector(j) = func(xj);
        
    end

    % Fill in the last row

    uMatrix(m,m-2) = 1/h*1/2;
    uMatrix(m,m-1) = -2/h;
    uMatrix(m,m) = 1/h*3/2;

    fVector(m) = 0;

    % Finally, solve uMatrix u = fVector
    uMatrix;
    u = uMatrix\fVector;
end