% Test 7: More complex case from loads and elastic thickness distribution
% as shown in figures loads.jpg and te.jpg. Digitized load and elastic 
% thickness contours (loadcontours.txt and tecontours.txt) are read into
% the program and then the data is interpolated in a grid using the
% gridfit function by John D'Errico. After that, the files necessary for
% flex3d are written

% Loads

load loadcontours.txt;
x=loadcontours(:,1)';
y=loadcontours(:,2)';
z=loadcontours(:,3)';

xg=0:5:500;
yg=0:5:500;

[ZGL,XG,YG] = gridfit(x,y,z,xg,yg);

% Don't allow negative loads
for i=1:size(yg,2)
    for j=1:size(xg,2)
        if ZGL(i,j) < 0.0
            ZGL(i,j) = 0.0;
        end
    end
end

% visualize interpolated loads
figure('Name','Interpolated loads','NumberTitle','off');
mesh(XG,YG,ZGL);
hold on;
% circles on plot are input data
plot3(x,y,z,'o');
xlabel('East km');
ylabel('North km');
zlabel('load km');
hold off;

% elastic thickness
load tecontours.txt;
xe=tecontours(:,1)';
ye=tecontours(:,2)';
ze=tecontours(:,3)';

% Interpolate elastic thickness to a larger plate: Add 250 km around map
% This is to avoid boundary effects
xg=-250:5:750;
yg=-250:5:750;

[ZGE,XG,YG] = gridfit(xe,ye,ze,xg,yg);

% Don't allow elastic thicknesses greater than 35 km and lower than
% 5 km
for i=1:size(yg,2)
    for j=1:size(xg,2)
        if ZGE(i,j) > 35.0
            ZGE(i,j) = 35.0;
        end
        if ZGE(i,j) < 5.0
            ZGE(i,j) = 5.0;
        end
    end
end

% visualize interpolated elastic thickness
figure('Name','Interpolated elastic thickness','NumberTitle','off');
mesh(XG,YG,ZGE);
hold on;
% circles on plot are input data
plot3(xe,ye,ze,'o');
xlabel('East km');
ylabel('North km');
zlabel('te km');
hold off;

% Write files for flex3d

% Loads in meters
% Notice that loads are located in the central part of the plate
% and that between the edge of the loads and the plate boundary
% there are 250 km. This is to avoid boundary effects.
fid = fopen('loads.txt','wt');

for i=1:size(yg,2)
    for j=1:size(xg,2)
        height = 0.0;
        density = 2700.0;
        if i > 50 && i < 152
            if j > 50 && j < 152
                height = ZGL(i-50,j-50)*1e3;
                % Make loads lower than 1 km, sedimentary loads with 
                % density 2400. Loads higher than 1 km are tectonic 
                % loads with density 2700
                if height < 1000.0
                    density = 2400.0;
                else
                    density = 2700.0;
                end
            end
        end
        fprintf(fid,'%f %f\n',height, density);
    end
end
fclose(fid);

% Elastic thickness in meters
fid = fopen('te.txt','wt');
for i=1:size(yg,2)
    for j=1:size(xg,2)
        elasticthickness = ZGE(i,j)*1e3;
        fprintf(fid,'%f\n',elasticthickness);
    end
end
fclose(fid);