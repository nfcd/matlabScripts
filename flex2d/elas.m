function elas

global E v h rigid
global flagelas


prompt = {'Enter the Young''s Modulus in Pa','Enter the Poisson''s ratio ','Enter the elastic thickness of the lithosphere in kilometers'};
def = {'7e10','0.25','30'};
title = 'INPUT ELASTIC PARAMETERS OF THE LITHOSPHERE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
E = str2double(answer{1});
v = str2double(answer{2});
h = str2double(answer{3})*1e3;
rigid = (E*h*h*h)/(12*(1.0-v*v));

flagelas=1;





