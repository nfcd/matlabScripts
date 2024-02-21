function refile

global flaggeom
global flaghede
global height pload


if flaggeom==0
   warndlg('Geometry of the loads is not defined','WARNING');
elseif flaggeom==1
	[datafile] = uigetfile('*.txt','Choose a txt file');
    A=load(datafile);
    height = A(:,1);
    pload = A(:,2);
    flaghede=1;
end



