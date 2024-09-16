% beta 0 - K
K = 10;

% increment
i = .01;

Vpot = @(x)(x.^2 - 1).^2;
 
Vder = @(x)4*x.*(x.^2 - 1);


inputstorage = zeros( (1/i) * K, 201);
betas = zeros( (1/i) * K, 201);

committorstorage = zeros((1/i) * K , 201);

total =  (1/i) * K;

% param*i is our beta value

%201 mesh

% each row is data

for param = 1 : total
    committorstorage(param, :) = committor_2Well_in_1D(param * i);
    inputstorage(param, :) = exp(Vpot( -1:.01:1) * param * i * - 1);
    betas(param, :) = zeros(1, 201) + param * i;
    tempo = param;

end



writematrix(inputstorage, 'thetainputs.csv')
writematrix(betas, "betasdoublewell.csv")
writematrix(committorstorage, 'thetacommittors.csv')




