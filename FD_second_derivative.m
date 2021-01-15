function result = FD_second_derivative(f, x, h)
% Take a second derivative using finite differences
% f: the function to differentiate
% x: the point at which to take the derivative
% h: the step size

    result = 1/h^2*(f(x+h)+f(x-h)-2*f(x));
end
