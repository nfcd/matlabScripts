% Copyright Nestor Cardozo 2009

% Test 3: 40 km constant elastic thickness

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
te = 40000.0;
fid = fopen('te.txt','wt');
for i=1:pointsy
    for j=1:pointsx
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