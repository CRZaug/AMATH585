function result = richardsonExtrapolation1(f,x,h)
    % One step of Richardson extrapolation
    result = 1/3*(4*FD_second_derivative(f,x,h/2)-FD_second_derivative(f,x,h));

end