function [w_delta_vec, delta_vec] = estimate_Delta(w_local_vec,Receptance_local_vec,index_Omega)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Calculation of delta as it is defined in the line-fit method.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS
% Cleaned by Ahmed Rashed

ind_vec=[1:index_Omega-1,index_Omega+1:length(w_local_vec)];
Den_vec=Receptance_local_vec(ind_vec)-Receptance_local_vec(index_Omega);
Num_vec=w_local_vec(ind_vec).^2-w_local_vec(index_Omega)^2;
w_delta_vec=w_local_vec(ind_vec);
delta_vec=Num_vec./Den_vec;