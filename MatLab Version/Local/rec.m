function [f_r_new_col,zeta_r_new_col,f_r_stabilized_col,f_zeta_r_stabilized_col]=rec(f_r_col,zeta_r_col,f_r_old_col,zeta_r_old_col,prec_f_r,prec_zeta_r)
% Modified by Ahmed Rashed
% Compare frequency and damping of the current step vs those of the previous step
%
%  Input data:
%  f_r_col: natural frequency vector,
%  zeta_r_col: damping vector,
%  prec_f_r: tolerance in frequency,
%  prec_zeta_r: tolerance in damping.
%
%  Output data:
%  f_r_new_col: updated eigenvalue matrix,
%  zeta_r_new_col: updated damping matrix,
%  bDampingFreq_Stabilized_r_new_vec: test matrix,
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT

f_r_new_col=f_r_col;
zeta_r_new_col=zeta_r_col;
f_r_stabilized_col=false(size(f_r_col));
f_zeta_r_stabilized_col=f_r_stabilized_col;

for n_r=1:length(f_r_col)
    [D_f_min,ind_min]=min(abs(f_r_col(n_r)-f_r_old_col));
    if (D_f_min/f_r_old_col(ind_min))<=prec_f_r  % Estimated natural-frequency stabilized
        f_r_new_col(n_r)=(f_r_col(n_r)+f_r_old_col(ind_min))/2;
        f_r_stabilized_col(n_r)=true;

        if (abs(zeta_r_col(n_r)-zeta_r_old_col(ind_min))/zeta_r_old_col(ind_min)) < prec_zeta_r  % Estimated natural-frequency & zeta stabilized
           zeta_r_new_col(n_r)=(zeta_r_col(n_r)+zeta_r_old_col(ind_min))/2;
           f_zeta_r_stabilized_col(n_r)=true;
        end
    end
end