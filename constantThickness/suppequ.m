% Copyright: Nestor Cardozo 2009. All rights reserved

function y=suppequ(gamstar,ramp)

% Equation 1 of Hardy (1997)

y = (1.+2.*cos(gamstar)*cos(gamstar))/sin(2.*gamstar) + (cos(ramp)-2.)/sin(ramp);