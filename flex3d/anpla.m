% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function anpla

%solves for the deformation of an infinite plate on an elastic foundation
%based on equation 179, pp 267 of Timoshenko and Woinowski 1959

global canvas
global flagload
global flagelas
global nnode
global xc 
global yc 
global ff
global emodule poisson t
global disp
global densup grav
global flagsol


if flagload==0
   msgbox('Loads are not specified');
elseif flagelas==0
   msgbox('Elastic Parameters are not specified');
else

tic;

set(canvas,'Pointer','watch');

disp=zeros(nnode,1);

D = (emodule*t^3)/(12*(1-poisson^2));

support = densup*grav;

alfa = (D/support)^(1/4);

for id = 1:nnode
   if ff(id)~=0.0
      for jd = 1:nnode
         X = abs(xc(id)-xc(jd));
         Y = abs(yc(id)-yc(jd));
         R = sqrt(X^2 + Y^2);
         XI = R/alfa;
         kei = pi/2*imag(bessely(0,XI*exp(1i*pi/4))) - pi/2*real(besselj(0,XI*exp(1i*pi/4)));
			disp(jd)= disp(jd)- ((ff(id)*alfa^2)/(2*pi*D))*kei;
         
      end
   end
end

set(canvas,'Pointer','arrow');

toc;

flagsol=1;

end