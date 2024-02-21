% Copyright: Nestor Cardozo 2012

%Print upper four beds of forward model with seven beds
fid = fopen('bed4.txt','wt');
for j=1:size(XP,2)
    fprintf(fid,'%f %f\n',XP(4,j),YP(4,j));
end
fclose(fid);

fid = fopen('bed5.txt','wt');
for j=1:size(XP,2)
    fprintf(fid,'%f %f\n',XP(5,j),YP(5,j));
end
fclose(fid);

fid = fopen('bed6.txt','wt');
for j=1:size(XP,2)
    fprintf(fid,'%f %f\n',XP(6,j),YP(6,j));
end
fclose(fid);

fid = fopen('bed7.txt','wt');
for j=1:size(XP,2)
    fprintf(fid,'%f %f\n',XP(7,j),YP(7,j));
end
fclose(fid);

