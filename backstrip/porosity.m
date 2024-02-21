function y = porosity(y2,y1,phi,c)


% calculates the porosity of a layer between depths y1 and y2 using an
% exponential decay with depth given initial porosity phi and decay constant c
% an if-else statement to avoid the division by zero error in case of an
% unconformity

if y2 == y1
	y = phi/c * (exp(-c*y1) - exp(-c*y2)) / 1;
else
	y = phi/c * (exp(-c*y1) - exp(-c*y2)) / (y2 - y1);
end