% beta 0 - K (100)
K = 0:.1:10; 

%a range (200)
A = -.9:.1:.9;

inputstorage = zeros(length(K)*length(A), 201);
betas = zeros(length(K)*length(A), 201);
outputstorage = zeros(length(K)*length(A), 201);

count = 1;

for i = 1:length(K)
    b = K(i); 
    for j = 1:length(A)
        a = A(j);
        outputstorage(count, :) = param_well_gen(a, b); 
        inputstorage(count, :) = exp(-1 * arrayfun(@(x) pieceVpot(x, a), -1:.01:1)*b); 
        betas(count, :) = zeros(1, 201) + b; 
        count = count + 1
    end
end 

writematrix(inputstorage, '2paraminputs.csv')
writematrix(betas, "2parambetas.csv")
writematrix(outputstorage, '2paramcommittors.csv')




