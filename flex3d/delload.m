% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function delload

global x_pos y_pos
global x_pix y_pix
global kmw kml
global flagload
global WIDTH LENGTH
global dence
global nel nnel
global xc 
global yc
global np
global ff
global grav
global hel
global dens
global elcol


if flagload==0
   msgbox('No Loads have been entered');
else

prompt = {'Enter 1 if background is image'};
def = {'1'};
title = 'MOUSE LEFT BUTTON = ENTER, RIGHT BUTTON = END';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
licon = str2num(answer{1});

hold on

if licon==1
   AW = (WIDTH - kmw)/2;
   AL = (LENGTH - kml)/2;
end

C =1;

XV = zeros(nnel,1);
YV = zeros(nnel,1);
XT = zeros(nnel,1);
YT = zeros(nnel,1);

while C==1
    
    [A,B,C] = ginput(1);
    
    if C==1
        
        if licon==0
    		range=WIDTH/(2*dence);
        elseif licon==1
    		range=WIDTH/(2*dence)*x_pix/kmw;
        end
        
        for ifor=1:nel
            if elcol(ifor)==1
                for jfor=1:nnel
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
                        elcol(ifor)=0;
                        if licon==0
                            fill(XV,YV,'w');
                        elseif licon==1
                            fill(XT,YT,'w');
                        end
                        w= abs(XV(1)-XV(2));
                        h= abs(YV(1)-YV(3));
                        load=(dens(ifor)*grav*hel(ifor)*w*h)/4;
                        hel(ifor)=0.0;
                        dens(ifor)=0.0;
       
                        for jfor=1:nnel
                            ff(np(ifor,jfor))=ff(np(ifor,jfor))-load;
                        end
                    end
                end
            end
        end
    end
end

hold off;

end