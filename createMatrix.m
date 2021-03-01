function [A,b] = createMatrix(n)

    h = 1/n;
    N = (n-1)^2;

    %  Form block tridiagonal finite difference matrix A and right-hand side 
    %  vector b.

    A=sparse(zeros(N,N));
    b = ones(N,1);         % Use right-hand side vector of all 1's.

    %  Loop over grid points in y direction.

    for j=1:n-1
      yj = j*h;
      yjph = yj+h/2;  yjmh = yj-h/2;

    %    Loop over grid points in x direction.

      for i=1:n-1
        xi = i*h;
        xiph = xi+h/2;  ximh = xi-h/2;
        aiphj = 1 + xiph^2 + yj^2;
        aimhj = 1 + ximh^2 + yj^2;
        aijph = 1 + xi^2 + yjph^2;
        aijmh = 1 + xi^2 + yjmh^2;
        k = (j-1)*(n-1) + i;
        A(k,k) = aiphj+aimhj+aijph+aijmh;
        if i > 1, A(k,k-1) = -aimhj; end
        if i < n-1, A(k,k+1) = -aiphj; end
        if j > 1, A(k,k-(n-1)) = -aijmh; end
        if j < n-1, A(k,k+(n-1)) = -aijph; end
      end
    end
    A = (1/h^2)*A;   % Remember to multiply A by (1/h^2).
end