function deriv = pieceVpot(x, a)
if x >= a
    deriv = (-2/(a-1)^3)*(x-1)^2*(x-((3*a-1)/2)); 
else
    deriv =  (-2/(a+1)^3)*(x+1)^2*(x-((3*a+1)/2)); 
end 
end

