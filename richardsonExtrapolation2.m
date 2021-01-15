function result = richardsonExtrapolation2(f,x,h)
    % Two steps of Richardson extrapolation
    result = 1/15*(16*richardsonExtrapolation1(f,x,h/2)-richardsonExtrapolation1(f,x,h));

end