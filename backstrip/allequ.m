function y = allequ(Y2,Y1,y2,y1,phi,c)

% equation 8.23 in Allen and Allen

y = Y2 - Y1 - (y2 - y1 - phi/c * (exp(-c*y1) - exp(-c*y2)) + phi/c *(exp(-c*Y1) - exp(-c*Y2)));
