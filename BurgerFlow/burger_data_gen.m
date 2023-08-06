N = 1400


inputstorage = zeros(N, 1600)

outputstorage = zeros(N, 1600)


for i = 1:N
    temp = Burgers_spectral(i+100) ;
    inputstorage(i, :) = temp(:, 1) ;
    outputstorage(i, :) = temp(:, 2) ;
    tempo = i

end

writematrix(inputstorage, "backwardsvisinputstept1sigmoid.csv")
writematrix(outputstorage, "backwardsvisoutputstept1sigmoid.csv")
