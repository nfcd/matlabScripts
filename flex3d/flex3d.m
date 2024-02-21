% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

%program flex3d

%the main program for estimating the deflection of a plate on an
%elastic foundation


global canvas
global flagdat
global flagima
global flagscal
global flagdisc
global flagelas
global flagrav
global flagfig
global flagload
global flagsol

canvas = figure('Name','FLEX3D');

flagdat=0;
flagima=0;
flagscal=0;
flagdisc=0;
flagelas=0;
flagrav=0;
flagfig=0;
flagload=0;
flagsol=0;

%setting the graphical user interface

HM_data = uimenu(gcf,'Label','Data');
HM_datload = uimenu(HM_data,'Label','Load Data', 'Callback','lomodres');
HM_datsave = uimenu(HM_data,'Label','Save Data', 'Callback','modresa');

HM_digit = uimenu(gcf,'Label','Image');
HM_digload = uimenu(HM_digit,'Label','Load Image','Callback','isomap');
HM_digscal = uimenu(HM_digit,'Label','Scale Image','Callback','isoscal');

HM_param = uimenu(gcf,'Label','Parameters');
HM_pardis = uimenu(HM_param,'Label','Set Plate Geometry','Callback','discre');
HM_parmec = uimenu(HM_param,'Label','Elastic Properties', 'Callback', 'elas');
HM_parres = uimenu(HM_param,'Label','Gravity and Restoring Force', 'Callback', 'resfor');


HM_load = uimenu(gcf,'Label','Loads');
HM_loin = uimenu(HM_load,'Label','Load Input', 'Callback', 'loadin');
HM_delo= uimenu(HM_load,'Label','Load Removal','Callback','delload');

HM_solve = uimenu(gcf,'Label','Solve');
HM_solan = uimenu(HM_solve,'Label','Classical Theory','Callback','anpla');

HM_plot = uimenu(gcf,'Label','Plot');
HM_plotdis = uimenu(HM_plot,'Label','Plate','Callback','plmesh');
HM_plotimes = uimenu(HM_plot,'Label','Plate and Image','Callback','immesh');
HM_plotin = uimenu(HM_plot,'Label','Undeformed Topography','Callback','plotin');
HM_fiplot = uimenu(HM_plot,'Label','Deformed Topography','Callback','fiplot');
HM_confi = uimenu(HM_plot,'Label','Topographic Contours','Callback','confi');