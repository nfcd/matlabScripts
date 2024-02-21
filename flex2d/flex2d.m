%program flex2d

%defining variables

global flagelas
global flagtyp
global flaggeom
global flaghede
global flagsol
global flagfig


flagelas=0;
flagtyp=0;
flaggeom=0;
flaghede=0;
flagsol=0;
flagfig=0;

%setting the graphical user interface

canvas = figure('Name','FLEX2D');

%------------------------------------------------------------------------

HM_par = uimenu(gcf,'Label','Parameters');
HM_parmec = uimenu(HM_par,'Label','Elastic properties', 'Callback', 'elas');
HM_parana = uimenu(HM_par,'Label','Type of analysis','Callback','tanal');

HM_geom = uimenu(gcf,'Label','Geometry');
HM_geomlo = uimenu(HM_geom,'Label','Load geometry','Callback','geom');
HM_geomde = uimenu(HM_geom,'Label','Heights and Densities of loads','Callback','hdens');
HM_geomfi = uimenu(HM_geom,'Label','Read loads from file','Callback','refile');

HM_sol = uimenu(gcf,'Label','Analytical Solution');
HM_solbe = uimenu(HM_sol,'Label','Solve','Callback','solbeams');
HM_solres = uimenu(HM_sol,'Label','Model Results','Callback','results');

HM_plot = uimenu(gcf,'Label','Plot');
HM_plotin = uimenu(HM_plot,'Label','Undeformed Topography','Callback','undef');
HM_plotout = uimenu(HM_plot,'Label','Deformed Topography','Callback','defor');


HM_help = uimenu(gcf,'Label','Help');
HM_helpab = uimenu(HM_help,'Label','About Flex2d','Callback','boxabout');
HM_helpbox = uimenu(HM_help,'Label','Help','Callback','boxhelp');

%-----------------------------------------------------------------------------

