function [A_r_phys,w_n_r_phys,zeta_r_phys]=mode_lsce(h_cols,D_t,Z_col,N_modes,N_t)
 
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Extraction of "physical modes", that's to say the modes corresponding to complex conjugate poles.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

if N_t<2*N_modes
    error('Number of samples insufficient for the requested number of modes!')
end

ind=find(Z_col==0,1);
V_r_col=Z_col(1:ind-1);

lambda_r=log(V_r_col)./D_t;
w_d_r=imag(lambda_r); 
delta_r=real(lambda_r); 
w_n_r=sqrt(w_d_r.^2+delta_r.^2); 
zeta_r=-(delta_r./w_n_r);

[w_n_r,i_sort]=sort(w_n_r/2/pi);
V_r_col=V_r_col(i_sort);
zeta_r=zeta_r(i_sort);

% Solving equation V_mat * A_r_temp= h [Maia, eqn 4.14]
V_mat=zeros(2*N_modes,2*N_modes); 
for ii=1:2*N_modes
    V_mat(ii,:)=V_r_col(1:2*N_modes).'.^(ii-1); 
end

[L,N_out]=size(h_cols); 
Ni=floor(L/N_t);
A_r=nan(2*N_modes,N_out);
for ind=1:N_out
    A_r_temp=zeros(2*N_modes,Ni);
    for ind2=1:Ni
        A_r_temp(:,ind2)=V_mat\h_cols(N_t*(ind2-1)+(1:2*N_modes),ind);
    end
    A_r(:,ind)=mean(A_r_temp,2);
    A_r(:,ind)=A_r(:,ind)/A_r(1,ind);   % Eigenvector normalization
end

% Extracting the physical modes (frequencies which appear two times)
[w_n_r_phys,i_w_n_r]=uniquetol(w_n_r,1e-8*2*pi/max(abs(w_n_r)));
zeta_r_phys=zeta_r(i_w_n_r);
A_r_phys=A_r(i_w_n_r,:);