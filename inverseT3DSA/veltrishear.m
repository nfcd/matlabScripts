% Copyright: Nestor Cardozo 2009. All rights reserved

function [vx,vy,vz] = veltrishear(xx,yy,vox,voz,Azp,At,htra,m)

% TRUE-3D
% BOUNDARY CONDITIONS, Eq. 5 of Cristallini et al. 
% GSA Bulletin v. 116, pp. 938-952 (2004)
if xx<0.
    if yy>=0.
		vx = vox;
        vy = 0.;
        vz = voz;
    elseif yy<0.
        vx = 0.;
        vy = 0.;
        vz = 0.;
    end
elseif xx>=0.
    if yy>=xx*m 
    	vx = vox;
        vy = 0.;
        vz = voz;
    elseif yy<=-xx*m
        vx = 0.;
        vy = 0.;
        vz = 0.;
    else
    	% Equation 11 of Cristallini et al. 2004
	 	Az = voz/2.0;
	 	msq = m*m;
	 	xxsq = xx*xx;
	 	yysq = yy*yy;
		% Equation 16 of Cristallini et al. 2004
		Jz = At/(msq*cos(htra)*cos(htra));
		% Equation 6 of Cristallini et al. 2004
		vx = (vox/2.0)*(yy/(xx * m) + 1.0);
		% Equation 7 of Cristallini et al. 2004
		vz = (voz/2.0)*(yy/(xx * m) + 1.0);
        % Lets drop the Azp*yy term
		Constante = -(vox*m)/4.0 + (Azp*xx*m)/2.0 - (Az*Jz*xx*msq)/2.0;
		vy = (vox*yysq)/(4.0*xxsq*m) - (Azp*yysq)/(2.0*xx*m) + (Az*Jz*yysq)/(2.0*xx) + Constante;
    end
end