%x is array?
function deriv = pieceVder(x, a)
if x >= a
    deriv = (-4/(a-1)^3)*(x-1).*(x-(3*a-1)/2) + (-2/(a-1)^3)*(x-1).^2;
else
    deriv = (-4/(a+1)^3)*(x+1).*(x-(3*a+1)/2) + (-2/(a+1)^3)*(x+1).^2;
end 
end

