%% AMATH 585 Homework 1
% Jan 15, 2020
%
% Camille Zaug
%
% Problems (3,4,5,6a) that did not require MATLAB are attached.
% Functions are also attached.

%% Problem 1.
% We compare a numerical implementation of the second derivative of
% $\sin(x)$ at pi/3 to the actual second derivative.

actual = -sin(pi/6);

hValues = [];
numericalResults = [];
error = [];

for exponent = 1:16

     h = 10^(-1*exponent);

    FD = FD_second_derivative(@sin, pi/6,h);
    
    hValues = [hValues, h];
    numericalResults = [numericalResults,FD];
    error = [error, actual-FD];
    
end


T = table(hValues',numericalResults',error');
T.Properties.VariableNames = {'h' 'NumericalResults' 'Error'};

%%
% The actual value of $-\sin(pi/6)$ is -0.5.
%
% The numerical results and error for different values of h are given by
% the table below.

T

%%
% From this table, we notice that as we decrease $h$ by an order of
% 1, the error appears to go down by an order of 2 (which makes sense, as
% the approximation was in fact described in the problem as second order
% accurate). This only occurs for
% $h=0.1, 0.01$, and $0.001$, however. After this point, the pattern by which
% the error decreases is less predictable. This is due to roundoff error.
% When $h$ is smaller than this, $h^2$ starts to become a signifcant fraction
% of $10e-16$ (machine precision), producing roundoff error in computations. 
% This is why calculations with very small
% values of $h$ produce such variable errors. This is especially evident
% when $h$ is extremely small. When this happens, the computed value of the
% second derivative varies so wildly that it jumps back and forth from 0 to
% 1e16!

%% Problem 2. 
% We repeat the goal of the previous problem: Compute $u''(x)$ where $u(x)
% = \sin(x)$ at $x=\pi/6$. Now, however, we use Richardson extrapolation to
% increase the order of accuracy. We first take the $h=0.2$, so $h/2 =
% 0.1$ and $h/4 = 0.05$.
%
% Below, we compute one step of Richardson extrapolation using two sets of
% $h$ values.

actual = -sin(pi/6);

RE1_hValues = [.2,.1];
RE1_numericalResults = [];
RE1_error = [];

for h = RE1_hValues
    
    FD = richardsonExtrapolation1(@sin, pi/6,h);

    RE1_numericalResults = [RE1_numericalResults,FD];
    RE1_error = [RE1_error, actual-FD];

end

RE1_T = table(RE1_hValues',RE1_numericalResults',RE1_error');
RE1_T.Properties.VariableNames = {'h' 'NumericalResults' 'Error'};

RE1_T
%%
% One step of Richardson extrapolation has more accuracy than our previous
% method of computing the second derivative. We can see from the table that
% reducing $h$ by a factor of 2 reduced the error by an order of magnitude!
% Furthermore, with the second order accurate method, when $h=0.1$, the
% error was on the order of $10e-4$. With one step of Richardson
% extrapolation and $h =0.1$ , the error is on the order of $10e-8$, which
% is twice the accuracy. It appears that one step of Richardson
% extrapolation is 4th order accurate, which is what we predict from studying the
% error of this method via Taylor series. This is a factor of 2 improvement from the
% previous method. 
%
% Now we'll repeat with 2 steps of Richardson extrapolation, which uses the
% results of the previous part (this is encoded into the function called
% below).

actual = -sin(pi/6);
h = .2;
  
FD = richardsonExtrapolation2(@sin, pi/6, h);

RE2_T = table(h,FD,actual-FD);
RE2_T.Properties.VariableNames = {'h' 'NumericalResults' 'Error'};

RE2_T
%% 
% With $h=0.2$, the error is on the order of $10e-11$ for two-step
% Richardson Extrapolation. This is four orders of improvement from one-step Richardson 
% extrapolation, which had an error of $10e-11$. Overall, it appears that two-step Richardson
% extrapolation is 6th order accurate, which is what we'd expect given a
% Taylor series analysis of this method.
%
% In summary of the what we've seen here, Richardson extrapolation is an excellent tool
% for reducing error. With this method, we avoid the issue of producing values of $h$ too
% near machine precision to produce meaningful values.

%% Problem 6.
% b,c) Code to solve the difference equations for the heat equation problem
% described in part a). Choosing $u = (1-x)^2$ to be the actual solution,
% so $f(x) = 2(3x^2-2x+1)$.
%
% The heat equation is solved for four different mesh sizes.

error = [];
hVect = [];
for k = 1:4
    n=10^k; % number of paritions
    a = 0; % left endpoint
    b = 1; % right endpoint
    h = 1/n; % Size of mesh

    % The RHS function
    f = @(x) 2*(3*x.^2-2*x+1); 
    x = linspace(a,b,n+1);
    x = x(2:end);

    % Numerically solve the heat equation
    result = heatEquation(a,b,h,f);

    % The true solution
    u = @(x) (1-x).^2; 
    actual = u(x);

    % Calculate the error
    
    e = sqrt(h*sum((actual'-result).^2));
    
    error = [error, e];
    hVect = [hVect, h];
end

heatT = table(hVect',error');
heatT.Properties.VariableNames = {'h' 'Error'};

heatT

%% 
% From the above table, it appears that this method is second-order
% accurate, which is what we expected. Decreasing $h$ by an order of 1  
% decreases the error by an order of 2. This was expected, as we used
% second-order approximations of all derivatives. 
