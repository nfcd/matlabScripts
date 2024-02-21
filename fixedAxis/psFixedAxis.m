% Copyright: Nestor Cardozo 2009. All rights reserved

% Program psFixedAxis: plots ps for fixed
% fault propagation folds

% Algorithm from Hardy and Poblet (1995)
% Marine and Petroleum Geology 12, 165-176

% ramp values
rampd=10:50;
ramp=rampd.*pi/180;

lbrat=zeros(size(ramp));

for i=1:size(ramp,2)
    % See Figure 1 B of Hardy and Poblet (1995)
    % Eq. 7 of Hardy and Poblet (1995)
    gam1=(pi-ramp(i))/2.;
    % Eq. 5 of Hardy and Poblet (1995). Assume no excess shear Sp = 0
    gamestar = acot((-1.*(2.*cos(ramp(i))-3.)/sin(ramp(i)))/2.);
    % Eq. 6 of Hardy and Poblet (1995)
    gamistar=gam1-gamestar;
    % Eq. 4 of Hardy and Poblet (1995)
    game=acot(cot(gamestar)-2.*cot(gam1));
    % Eq. 3 of Hardy and Poblet (1995)
    gami = asin((sin(gamistar)*sin(game))/sin(gamestar));
    % Eq. 8 of Hardy and Poblet (1995)
    % Ratio of backlimb length to total slip
    % This is the same as the propagation to slip ratio
    a1=cot(gamestar)-cot(gam1);
    a2=1./sin(ramp(i))-(sin(gami)/sin(game))/sin(game+gami-ramp(i));
    a3=sin(gam1+ramp(i))/sin(gam1);
    lbrat(i)=a1/a2 + a3;
end

plot(rampd,lbrat,'b-');
axis ([10 50 0 6]);
xlabel('Ramp angle in degrees');
ylabel('P by S');
grid on;



