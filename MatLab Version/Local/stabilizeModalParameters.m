function [f_r_new_col,zeta_r_new_col,stabilized_f_r_col,stabilized_f_zeta_r_col]=stabilizeModalParameters(f_r_col,zeta_r_col,f_r_old_col,zeta_r_old_col,prec_f_r,prec_zeta_r)
% Created by Ahmed Rashed based on EasyMod "rec.m"

f_r_new_col=f_r_col;
zeta_r_new_col=zeta_r_col;
stabilized_f_r_col=false(size(f_r_col));
stabilized_f_zeta_r_col=stabilized_f_r_col;

for n_r=1:length(f_r_col)
    [D_f_min,ind_min]=min(abs(f_r_col(n_r)-f_r_old_col));
    if (D_f_min/f_r_old_col(ind_min))<=prec_f_r  % Estimated natural-frequency stabilized
        f_r_new_col(n_r)=(f_r_col(n_r)+f_r_old_col(ind_min))/2;
        stabilized_f_r_col(n_r)=true;

        if (abs(zeta_r_col(n_r)-zeta_r_old_col(ind_min))/zeta_r_old_col(ind_min)) < prec_zeta_r  % Estimated natural-frequency & zeta stabilized
           zeta_r_new_col(n_r)=(zeta_r_col(n_r)+zeta_r_old_col(ind_min))/2;
           stabilized_f_zeta_r_col(n_r)=true;
        end
    end
end