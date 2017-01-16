function [ww_delta, delta]  =  estimate_delta(ww_local,H_local,Omega)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Calculation of delta as it is defined in the line-fit method.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


index_Omega = find(ww_local==Omega) ;
H_omega = H_local(index_Omega,1) ;
jj = 0 ;
for index = 1:length(ww_local)
    if index == index_Omega
        continue
    else
    jj = jj+1 ;
    Den(jj,1) = H_local(index,1)-H_omega ;
    Num(jj,1) = ww_local(index,1).^2-Omega^2 ;
    ww_delta(jj,1) = ww_local(index,1) ;
    end
end
delta = Num./Den ;
