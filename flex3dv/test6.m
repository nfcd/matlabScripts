% Copyright Nestor Cardozo 2009

% Test 6: Variable elastic thickness
% At the center of the plate the elastic thickness
% is 10 km. From the center, the elastic thickness 
% increases linearly to 30 km at a distance of 250 km. 
% At farther distances the elastic thickness
% is constant and is 30 km.

% 3000 m high load at the center of the plate
% Then it decreases linearly to 0 at
% 100 km from the center of the plate

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
        yval = (i-1)*delta;
        xval = (j-1)*delta;
        radius = sqrt((xval-500)^2+(yval-500)^2);
        if radius < 250.0
            te = 10000.0 + radius/250*20000.0;
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
        radius = sqrt((xval-500)^2+(yval-500)^2);
        if radius < 100
            hei_load = 3000.0 * (1.0 - radius/100);
        end
        fprintf(fid,'%f %f\n',hei_load,rho_load);
    end
end        
fclose(fid);