function [w_delta_vec, delta_vec] = estimate_delta(w_local_vec,Receptance_local_vec,index_w)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Calculation of delta as it is defined in the line-fit method.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


ind_vec=[1:index_w-1,index_w+1:length(w_local_vec)];
jj=1;
for index=1:length(w_local_vec)
    if index == index_w
        continue
    end
    
    Den_vec(jj)=Receptance_local_vec(index,1)-Receptance_local_vec(index_w);
    Num_vec(jj)=w_local_vec(index,1).^2-w_local_vec(index_w)^2;
    w_delta_vec(jj)=w_local_vec(index);

    jj=jj+1;
end
delta_vec=Num_vec./Den_vec;