% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function loadin

global x_pos y_pos
global x_pix y_pix
global kmw kml
global flagrav
global flagfig
global flagload
global WIDTH LENGTH
global dence
global nel nnel
global xc
global yc
global np
global elcol
global hel
global dens
global ff
global grav

if flagfig==0
   msgbox('No Plate on display');
elseif flagrav==0  
   msgbox('Gravity is not specified');
else
   
prompt = {'Enter 1 if background is image'};
def = {'1'};
title = 'MOUSE LEFT BUTTON = ENTER, RIGHT BUTTON = END';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
licon = str2num(answer{1});
 	 
hold on;

if licon==1
   AW = (WIDTH - kmw)/2;
   AL = (LENGTH - kml)/2;
end

C =1;

XV = zeros(nnel,1);
YV = zeros(nnel,1);
XT = zeros(nnel,1);
YT = zeros(nnel,1);

while C==1;
    
    [A,B,C] = ginput(1);
    
    if C==1
        prompt = {'Enter load height in m','Enter load density in kg/m^3'};
 	 	def = {'1000','2700'};
 	 	title = 'Load Specifications';
    	lineNo=1;
    	answer=inputdlg(prompt,title,lineNo,def);
 	 	height = str2double(answer{1,1})*1e-3;
 	 	density = str2double(answer{2,1});
        if licon==0
    		range=WIDTH/(2*dence);
        elseif licon==1
    		range=WIDTH/(2*dence)*x_pix/kmw ;
        end
        
        for ifor=1:1:nel
            if elcol(ifor)==0
                for jfor=1:1:nnel
                    XV(jfor)= xc(np(ifor,jfor));
                    YV(jfor)= yc(np(ifor,jfor));
                    if licon==1
                        XT(jfor) = (xc(np(ifor,jfor))-AW)*x_pix/kmw + x_pos(1);
                        YT(jfor) = (LENGTH-AL-yc(np(ifor,jfor)))*y_pix/kml + y_pos(1);
                    end
                end
                if licon==0
                    diff = abs(XV(1)-A);
                    diffv = abs(YV(1)-B);
                elseif licon==1
                    diff = abs(XT(1)-A);
                    diffv = abs(YT(1)-B);
                end
                if diff<=range && diffv<=range
                    if licon==0 
                        in = inpolygon(A,B,XV,YV);
                    elseif licon==1
                        in = inpolygon(A,B,XT,YT);
                    end
                    if in==1
                        elcol(ifor)=1;       
                        hel(ifor)=hel(ifor)+height;
                        dens(ifor)=dens(ifor)+density;
                        if licon==0
                            if dens(ifor)>2400
                                fill(XV,YV,'g');
                            elseif dens(ifor)<=2400
                                fill(XV,YV,'y');
                            end
                        elseif licon==1
                            if dens(ifor)>2400
                                fill(XT,YT,'g');
                            elseif dens(ifor)<=2400
                                fill(XT,YT,'y');
                            end
                        end
                        w= abs(XV(1)-XV(2));
                        h= abs(YV(1)-YV(3));
                        load=(density*grav*height*w*h)/4;
                        for jfor=1:1:nnel
                            ff(np(ifor,jfor))= ff(np(ifor,jfor))+load; 
                        end
                    end
                end   
            end
        end
    end
  
end

flagload=1;
hold off;

end