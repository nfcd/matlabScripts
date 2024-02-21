% Copyright: Nestor Cardozo 2009. All rights reserved

% program runmodel

% 3D
% Run trishear model

% INITIAL GEOMETRY
% Extent in X and Y, cell size, and origin
prompt = {'Extent of section in the X (East) direction','Extent of section in the Y (North) direction','Cell size','X Origin','Y Origin'};
def = {'500.0','500.0','5.0','0.0','0.0'};
title = 'X-Y EXTENT, CELL SIZE, AND ORIGIN';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
extentx = str2double(answer{1});
extenty = str2double(answer{2});
cellsize = str2double(answer{3});
origx = str2double(answer{4});
origy = str2double(answer{5});

% Datums (Z increases upwards)
prompt = {'Datums, Z increases upwards'};
def = {'0.0 30.0 60.0 90.0 120.0'};
title = 'DATUMS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
datums = str2num(answer{1});

% Make beds (surfaces)
[XP,YP,ZP]=meshgrid(origx:cellsize:extentx+origx,origy:cellsize:extenty+origy,datums);

% INPUT PARAMETERS OF THE TRISHEAR MODEL

% Strike and dip of fault: Right hand rule is followed (Fault plane dips to the right of its strike direction)
% Only planar faults can be modeled with this algorithm
prompt = {'Strike','Ramp Angle'};
def = {'180.0','30.0'};
title = 'STRIKE AND DIP OF FAULT';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
strike = str2double(answer{1})*pi/180.0;
ramp = str2double(answer{2})*pi/180.0;

% FAULT TIP LINE
% fault tip 1 and length of fault along strike are given, and fault tip 2 is computed from % fault tip 1 along the strike direction
prompt = {'X (east)fault tip 1','Y (north) fault tip 1','Z (up) fault tip 1','Length of fault along strike'};
def = {'200.0','500.0','0.0','500.0'};
title = 'INITIAL FAULT TIP LOCATION';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
xts = str2double(answer{1});
yts = str2double(answer{2});
zts = str2double(answer{3});
faultLength = str2double(answer{4});
xtn = xts+sin(strike)*faultLength;
ytn = yts + cos(strike)*faultLength;
ztn = zts;

% Propagation to Slip ratio (P/S)
prompt = {'P/S at tip 1','P/S at tip 2'};
def = {'1.5','1.5'};
title = 'PROPAGATION TO SLIP RATIO';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
pss = str2double(answer{1});
psn = str2double(answer{2});

% Trishear angle (angle of triangular zone of deformation)
prompt = {'Trishear angle at tip 1','Trishear angle at tip 2'};
def = {'60.0','60.0'};
title = 'TRISHEAR ANGLE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
tras = str2double(answer{1})*pi/180;
tran = str2double(answer{2})*pi/180;

% Fault Slip
prompt = {'Fault slip at tip 1','Fault slip at tip 2','Number of slip increments','Rake of fault slip (0-180)'};
def = {'100.0','100.0','100','90.0'};
title = 'FAULT SLIP, + IF REVERSE, - IF NORMAL';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
sls = str2double(answer{1});
sln = str2double(answer{2});
ninc = str2double(answer{3});
slrake = str2double(answer{4})*pi/180;

forwtrishear;