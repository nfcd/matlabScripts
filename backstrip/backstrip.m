%program backstrip

%defining the variables

global flagfile
global flagunits
global flagage
global flagdens
global flagporo
global flagdecom
global flagsedwe
global flagwater
global canvas

%flags to ensure success flow

flagfile = 0;
flagunits = 0;
flagage = 0;
flagdens = 0;
flagporo = 0;
flagdecom = 0;
flagsedwe = 0;
flagwater = 0;

%setting the graphical user interface

canvas = figure('Name','BACKSTRIP');

%------------------------------------------------------------------------
HM_par = uimenu(gcf,'Label','Parameters');
HM_parunits = uimenu(HM_par,'Label','Geometry','Callback','inunits');
HM_parage = uimenu(HM_par,'Label','Ages','Callback','agein');
HM_pardens = uimenu(HM_par,'Label','Densities', 'Callback', 'densit');
HM_parpor = uimenu(HM_par,'Label','Porosity Parameters', 'Callback', 'porpar');
HM_parfil = uimenu(HM_par,'Label','Read from file','Callback','loparfi');

HM_sol = uimenu(gcf,'Label','Solve');
HM_soldec = uimenu(HM_sol,'Label','Decompact Sediments','Callback','decomp');
HM_soltec = uimenu(HM_sol,'Label','Sediment Load Correction','Callback','sedim');
HM_solwat = uimenu(HM_sol,'Label','Water depth correction','Callback','water');

HM_plot = uimenu(gcf,'Label','Plot');
HM_plotdecd = uimenu(HM_plot,'Label','Decompacted depth vs time','Callback','pdecdep');
HM_plotsub = uimenu(HM_plot,'Label','Basement depth without water correction','Callback','ptecwow');
HM_plotwat = uimenu(HM_plot,'Label','Tectonic Subsidence','Callback','ptecsub');

%----------------------------------------------------------------------------------


