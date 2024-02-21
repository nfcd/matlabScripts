% Copyright Nestor Cardozo 2009

% Test 4: Variable elastic thickness
% Elastic thickness is 30 km for x <= 250 km
% and x >= 750 km, and decreases linearly to
% 10 km at x = 500 km

% Load extends along y from 350 km to 650 km,
% and along x from 400 to 600 km. Load height
% varies along x, from 3.0 km at the center
% of the plate, to 0 at 100 km from the center 
% of the plate

% Spacing between points = 5.0 km

% number of points in x and y
pointsx = 201;
pointsy = 201;
% grid size in km
delta = 5.0;

% elastic thickness in meters
fid = fopen('te.txt','wt');
for i=1:pointsy
    for j=1:pointsx
        te = 30000.0;
        xval = (j-1)*delta;
        if xval > 250 && xval < 750
            te= 10000.0 + abs(xval-500)/250*20000.0;
        end
        fprintf(fid,'%f\n',te);
    end
end        
fclose(fid);

% density of load in kg/m^3
rho_load = 2700.0;
fid = fopen('loads.txt','wt');
for i=1:pointsy
    for j=1:pointsx
        hei_load = 0.0;
        yval = (i-1)*delta;
        xval = (j-1)*delta;
        if yval > 350. && yval < 650.
        		if xval > 400. && xval < 600.
        			% load height in meters
        			hei_load = 3000.0 * (1 - abs(xval-500)/100);
        		end
        end
        fprintf(fid,'%f %f\n',hei_load,rho_load);
    end
end        
fclose(fid);