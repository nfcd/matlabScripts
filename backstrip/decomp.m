function decomp


%decompact sedimentary column

global zt zb 
global c surpor
global Y1 Y2 
global flagunits
global flagage
global flagporo
global flagdens
global flagdecom
global canvas

if flagunits == 0
   warndlg('Missing tops and bases', 'WARNING');
elseif flagage==0
   warndlg('Missing ages', 'WARNING');
elseif flagporo==0  
   warndlg('Missing porosity parameters', 'WARNING');
elseif flagdens==0   
   warndlg('Missing sediment dry densities', 'WARNING');
else
   
% to avoid division by 0

for i=1:size(c,1)
   if c(i)==0 
      c(i)=0.0001; 
   end
end

% setup initial Y1 and Y2 matrices

for i = 1:size(c,1)
	for j = i:-1:1
        if j == i 
            Y1(j,i) = 0; 
			Y2(j,i) = zb(j) - zt(j);
        else 
            Y1(j,i) = Y2(j+1,i);
			Y2(j,i) = zb(j) - zt(j) + Y1(j,i);
        end
	end
end

% solve equation 8.23 in Allen and Allen

set(canvas,'Pointer','watch');

for i = 1:size(c,1)
	for j = i:-1:1
        if j == i 
            Y1(j,i) = 0; 
        else 
            Y1(j,i) = Y2(j+1,i);
        end
        options=optimset('display','off');
        yy = fzero('allequ',Y2(j,i),options,Y1(j,i),zb(j),zt(j),surpor(j),c(j));
        Y2(j,i) = yy;
	end
end

set(canvas,'Pointer','arrow');

flagdecom = 1;

end



