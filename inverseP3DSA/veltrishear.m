% Copyright: Nestor Cardozo 2009. All rights reserved

function [vx, vy] = veltrishear(xx,yy,sinc,m)

% 2D
% simplest velocity of trishear model
% algorithm from Zehnder and Allmendinger (2000)

if xx <0.0
    if yy >=0.0
        vx = sinc;
        vy = 0.0;
    elseif yy<0.0
        vx=0.0;
        vy=0.0;
    end
elseif xx>=0.0
    if yy>=xx*m 
        vx=sinc;
        vy=0.0;
    elseif yy<=-xx*m
        vx=0.0;
        vy=0.0;
    else
        % EQUATION 6 OF ZEHNDER AND ALLMENDINGER (2000)
        vx=(sinc/2.0)*(yy/(xx*m)+1.0);
        vy=(sinc/2.0)*(m/2)*((yy/(xx*m))^2.0-1.0);
    end
end