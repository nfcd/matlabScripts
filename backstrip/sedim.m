function sedim

% SEDIMENT LOAD CORRECTION (Ym) assuming local compensation (Airy isostasy)
% rhos = effective density of the entire sedimentary column at a given time
% S = the total thickness of the entire sedimentary column at a given time

global psed pm pw
global c surpor
global basset
global Y1 Y2
global S Ym
global flagfile
global flagdecom
global flagsedwe


if flagdecom==0
   warndlg('Sediments are not decompacted','WARNING');
else
   
if flagfile==0
   prompt = {'1 if filling a marine basin, 2 if the basin is continental'};
	def = {'1'};
	title = 'BASIN SETTING FOR EACH UNIT';
	lineNo = 10;
	answer=inputdlg(prompt,title,lineNo,def);
	basset  = str2num(answer{1});
end

a = size(psed,1);

%estimate thickness 

thick = Y2 - Y1;
rhos = zeros(a,1);
S = Y2(1,:);

for i = 1:a
    for j = i:-1:1
		phi = porosity(Y2(j,i),Y1(j,i),surpor(j),c(j));
		pb = phi * pw + (1-phi) * psed(j);
		rhos(i) = rhos(i) + pb * thick(j,i) / S(i);
    end
    if basset(i)==1
        Ym(i) = S(i) * (pm - rhos(i)) / (pm - pw);
    elseif basset(i)==2
        Ym(i) = S(i) * (pm - rhos(i)) / (pm);
    end
end

flagsedwe = 1;

end




