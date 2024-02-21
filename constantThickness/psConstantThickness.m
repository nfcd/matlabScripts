% Copyright: Nestor Cardozo 2009. All rights reserved

% Program psConstantThickness: plots ps for constant thickness
% fault propagation folds

% Algorithm from Hardy (1997)
% Journal of Structural Geology 19, 893-896

% ramp values
rampd=10:50;
ramp=rampd.*pi/180;

lbrat=zeros(size(ramp));

for i=1:size(ramp,2)
    % SOLVE MAIN PARAMETERS
    % See figure 1 of Hardy (1997)
    % Solve eq. 1 of Hardy (1997)
    options=optimset('display','off');
    gamstar = fzero('suppequ',0.5,options,ramp(i));
    % Eq. 2 of Hardy (1997)
    gam1 = pi/2. - ramp(i)/2.;
    % Eq. 3 of Hardy (1997)
    gam = pi/2.+gamstar-gam1;
    % Eq. 8 of Hardy (1997)
    % Propagation to slip ratio
    lbrat(i) = 1./(1.-sin(ramp(i))/sin(2.*gam-ramp(i)));
end

plot(rampd,lbrat,'b-');
axis ([10 50 0 6]);
xlabel('Ramp angle in degrees');
ylabel('P by S');
grid on;