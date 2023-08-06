function pp = ChebClenshaw(x,fx,xx)
% Compute Chebyshev coefficients and evaluate Chebyshev sum using
% Clenshaw's method
n = length(x) - 1;
%% compute Chebyshev coefficients
theta = pi*(0:n)'/n; % angles theta
c = zeros(n+1,1);
for k = 0 : n
    c(k + 1) = sum(fx.*cos(k*theta))-0.5*(fx(1)*cos(k*theta(1))+fx(end)*cos(k*theta(end)));
end
c = 2*c/n;
fprintf('c = [\n');
fprintf('%d\n',c);
fprintf(']\n');
%% evaluate Chebyshev sum at points t using Clenshaw's method
e = ones(n+1,1);
pp = zeros(size(xx));
for j = 1 : length(xx)
    A = spdiags([e,-xx(j)*2*e,e],-2:0,n+1,n+1);
    bb = A'\c;
    pp(j) = 0.5*(bb(1) - bb(3)) - 0.5*c(end)*cos(n*acos(xx(j)));
end
end
