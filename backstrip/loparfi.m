function loparfi


global psed pm pw
global zt zb
global agetop agebas
global c surpor
global basset
global flagfile
global flagunits
global flagage
global flagdens
global flagporo


%define main parameters for problem

prompt = {'Water density in kg/m^3','Mantle density in kg/m^3'};
def = {'1000.0','3300.0'};
title = 'MAIN PARAMETERS';
lineNo = 1;
answer=inputdlg(prompt,title,lineNo,def);
pw = str2double(answer{1});
pm = str2double(answer{2});

%read parameters from a txt file

[datafile] = uigetfile('*.txt','Choose a txt file');
A=load(datafile);
zt = A(:,1);
zb = A(:,2);
agetop = A(:,3);
agebas = A(:,4);
psed = A(:,5);
c = A(:,6);
surpor = A(:,7);
basset = A(:,8);

flagfile = 1;
flagunits = 1;
flagage = 1;
flagdens = 1;
flagporo = 1;

   







