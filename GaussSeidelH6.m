function x = GaussSeidelH6(A,b,x,iter)


M = tril(A);

x_new = ones(length(x),1);

for i = 1:iter
    
    r = b-A*x;
    
    z = M\r;
    x_new = x+z;
    
    e = sum(abs(x-x_new));

    x = x_new;
    
end

end
