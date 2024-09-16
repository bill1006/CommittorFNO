function sol = param_well_gen(a,beta)
close all
% solves u'' - beta*(4*x.*(x.^2 - 1)).*u' = 0
% u(-1) = 0, u(1) = 1
%%
ubdry = @(x) 0.5*(1 + x);
x_eval = -1:0.01:1; % points for solution evaluation

N = 750; 

[D,x] = cheb(N); % the Chebyshev differentiation matrix and the Chebyshev grid
D2 = D^2; % the differentiation matrix for the 2nd derivative
D = D(2:N,2:N); % remove boundary rows and columns
D2 = D2(2:N,2:N); % remove boundary rows and columns
L = D2 - beta*(arrayfun(@(x) pieceVder(x, a),x(2:N)))*ones(1,N-1).*D; % the differential operator
f = beta*0.5*arrayfun(@(x) pieceVder(x, a) ,x(2:N)); % -L applied to the boundary function ubdry = 0.5*(1+x) computed exactly.
u = L\f; % solve for the interior points
u = [0;u;0] + ubdry(x); % add the boundary points 
    % evaluate the solution at points x_eval     
%    uu = polyval(polyfit(x,u,N),x_eval); 
 uu = ChebClenshaw(x,u,x_eval); 
 sol = uu;
   


% takes input 

