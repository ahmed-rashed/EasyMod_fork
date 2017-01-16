function [a,b] = moindre_carre(x,y)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Best-line fit using the least squares method.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


N = length(x) ; 
Sx = sum(x) ; 
Sy = sum(y) ; 
Sxx = sum(x.*x) ; 
Sxy = sum(x.*y) ; 
b = (Sxx*Sy-Sx*Sxy)/(N*Sxx-Sx^2) ; 
a = (N*Sxy-Sx*Sy)/(N*Sxx-Sx^2) ; 

