function well_2params(a,beta)
close all
% solves u'' - beta*(4*x.*(x.^2 - 1)).*u' = 0
% u(-1) = 0, u(1) = 1
%%
ubdry = @(x) 0.5*(1 + x);
x_eval = -1:0.01:1; % points for solution evaluation

Ns = 20 : 50 : 820;
NN = length(Ns);
sol = zeros(length(x_eval),NN);
for k = 1 : NN
    N = Ns(k);
    [D,x] = cheb(N); % the Chebyshev differentiation matrix and the Chebyshev grid
    D2 = D^2; % the differentiation matrix for the 2nd derivative
    D = D(2:N,2:N); % remove boundary rows and columns
    D2 = D2(2:N,2:N); % remove boundary rows and columns
    L = D2 - beta*(arrayfun(@(x) pieceVder(x, a),x(2:N)))*ones(1,N-1).*D % the differential operator
    f = beta*0.5*arrayfun(@(x) pieceVder(x, a) ,x(2:N)) % -L applied to the boundary function ubdry = 0.5*(1+x) computed exactly.
    u = L\f; % solve for the interior points
    u = [0;u;0] + ubdry(x); % add the boundary points 
    % evaluate the solution at points x_eval     
%    uu = polyval(polyfit(x,u,N),x_eval); 
     uu = ChebClenshaw(x,u,x_eval); 
    sol(:,k) = uu;
end
figure;
hold on;
fsz = 20;
for k = 1 : NN-1
    lname = sprintf("N = %d",Ns(k));
    plot(x_eval,abs(sol(:,k) - sol(:,NN)),'linewidth',2,'DisplayName',lname);
end
ylabel('Error','FontSize',fsz);
xlabel('x','Fontsize',fsz);
legend;
set(gca,'Fontsize',fsz,'YScale','log');
%
figure;
hold on;
fsz = 20;
for k = 1 : NN
    lname = sprintf("N = %d",Ns(k));
    plot(x_eval,sol(:,k),'linewidth',2,'DisplayName',lname);
end
ylabel('q','FontSize',fsz);
xlabel('x','Fontsize',fsz);
legend;
set(gca,'Fontsize',fsz);
end


% takes input 

