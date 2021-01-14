function [f_r_new_col,zeta_r_new_col,bDampingFreq_Stabilized_r_new_col]=rec(f_r_col,zeta_r_col,f_r_old_col,zeta_r_old_col,prec_f_r,prec_zeta_r)
% Modified by Ahmed Rashed
% Compare frequency and damping of the current step vs those of the previous step
%
%  Input data:
%  f_r_col: natural frequency vector,
%  zeta_r_col: damping vector,
%  f_max: maximum freqneucy covered by the data (for directly eliminating the frequency out of range),
%  prec_f_r: tolerance in frequency,
%  prec_zeta_r: tolerance in damping.
%
%  Output data:
%  f_r_new_col: updated eigenvalue matrix,
%  zeta_r_new_col: updated damping matrix,
%  bDampingFreq_Stabilized_r_new_vec: test matrix,
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT

f_r_new_temp1_col=nan(size(f_r_old_col));
zeta_r_new_temp1_col=f_r_new_temp1_col;
f_r_new_temp2_col=f_r_new_temp1_col;
zeta_r_new_temp2_col=zeta_r_new_temp1_col;
bDampingFreq_Stabilized_r_new_temp1_col=false(size(f_r_old_col));

N_new=0;
for n_r=1:length(f_r_col)
    [D_f_min,ind_min]=min(abs(f_r_col(n_r)-f_r_old_col));
    if (D_f_min/f_r_old_col(ind_min))<=prec_f_r  % Estimated natural frequency stabilized
        f_r_new_temp1_col(ind_min)=(f_r_col(n_r)+f_r_old_col(ind_min))/2;

        if (abs(zeta_r_col(n_r)-zeta_r_old_col(ind_min))/zeta_r_old_col(ind_min)) < prec_zeta_r  % Estimated zeta stabilized
           zeta_r_new_temp1_col(ind_min)=(zeta_r_col(n_r)+zeta_r_old_col(ind_min))/2;
           bDampingFreq_Stabilized_r_new_temp1_col(ind_min)=true;
        else
           zeta_r_new_temp1_col(ind_min)=zeta_r_col(n_r);
           bDampingFreq_Stabilized_r_new_temp1_col(ind_min)=false;
        end
    else  %f_r_col(n_r) does not match any of the old estimated natural frequencies
        N_new=N_new+1;
        f_r_new_temp2_col(N_new)=f_r_col(n_r);
        zeta_r_new_temp2_col(N_new)=zeta_r_col(n_r);
    end
end

ind_f_r_1_col=find(~isnan(f_r_new_temp1_col));
f_r_new_col=[f_r_new_temp1_col(ind_f_r_1_col);f_r_new_temp2_col(1:N_new)];
zeta_r_new_col=[zeta_r_new_temp1_col(ind_f_r_1_col);zeta_r_new_temp2_col(1:N_new)];
bDampingFreq_Stabilized_r_new_col=[bDampingFreq_Stabilized_r_new_temp1_col(ind_f_r_1_col);false(N_new,1)];