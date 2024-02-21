% Copyright: Nestor Cardozo 2009. All rights reserved

function y=suppequ(gama,ramp)

% Equation 1 of Hardy (1995)

y = sin(2*gama)/(2*(cos(gama))^2+1) - tan(ramp);

