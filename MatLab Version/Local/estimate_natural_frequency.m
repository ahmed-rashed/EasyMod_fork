function [fr,theta_r]  =  estimate_natural_frequency(ind_max,delta2,f_delta2,theta_delta2)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Searching the natural frequency using the method of radial deviation.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


for jj = 1:length(ind_max)
    xa = f_delta2(ind_max(jj)-1) ;
    xb = f_delta2(ind_max(jj)) ;
    ya = delta2(ind_max(jj)-1) ;
    yb = delta2(ind_max(jj)) ;
    a = (yb-ya)/(xb-xa) ;
    b = yb-a*xb ;
    P = [a b] ;
    fr(jj) = roots(P) ;
    f1 = f_delta2(ind_max(jj)-1) ;
    f2 = f_delta2(ind_max(jj)) ;
    theta1 = theta_delta2(ind_max(jj)-1) ;
    theta2 = theta_delta2(ind_max(jj)) ;
    a = (theta2-theta1)/(f2-f1) ;
    b = theta2-a*f2 ;
    theta_r(jj) = a*fr(jj)+b ;
end


