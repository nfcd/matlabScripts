% Copyright: Nestor Cardozo 2018

function [vx, vy] = veltrishear(xx,yy,m,sinc)

% 2D simplest velocity of trishear model
% algorithm from Zehnder and Allmendinger (2000)
% Journal of Structural Geology 22, 1009-1014

if xx <0.
    if yy >=0.
        vx = sinc;
        vy = 0.;
    elseif yy<0.
        vx=0.;
        vy=0.;
    end
elseif xx>=0.
    if yy>=xx*m 
        vx=sinc;
        vy=0.;
    elseif yy<=-xx*m
        vx=0.;
        vy=0.;
    else
        % EQUATION 6 OF ZEHNDER AND ALLMENDINGER (2000)
        vx=(sinc/2.)*(yy/(xx*m)+1.);
        vy=(sinc/2.)*(m/2.)*((yy/(xx*m))^2.-1.);
    end
end