function [A_r,w_r_col,zeta_r_col]=mode_lsce(h_cols,D_t,Z_col,N_modes)
 
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Extraction of "physical modes", that's to say the modes corresponding to complex conjugate poles.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

[N_t,N_outputs]=size(h_cols); 
N_inputs=floor(N_t/N_t); %Bug is here
if N_t<2*N_modes
    error('Number of samples insufficient for the requested number of modes!')
end

lambda_r_col=log(Z_col)/D_t;
w_d_r_col_temp=imag(lambda_r_col); 
delta_r_col_temp=real(lambda_r_col); 
w_r_col_temp=sqrt(w_d_r_col_temp.^2+delta_r_col_temp.^2); 
zeta_r_col_temp=-(delta_r_col_temp./w_r_col_temp);

[w_r_col_temp,i_sort]=sort(w_r_col_temp);
Z_col=Z_col(i_sort);
zeta_r_col_temp=zeta_r_col_temp(i_sort);

% Extracting the physical modes (frequencies which appear in pairs)
[w_r_col,i_w_n_r]=uniquetol(w_r_col_temp,eps*1e3*max(abs(w_r_col_temp)));
zeta_r_col=zeta_r_col_temp(i_w_n_r);

% Solving equation W_V_mat * A_r_temp= h [Maia, eqn (4.48) or (4.14)]
W_V_mat=zeros(2*N_modes,2*N_modes); 
for n_mode=1:2*N_modes
    W_V_mat(n_mode,:)=Z_col(1:2*N_modes).'.^(n_mode-1); %This is bug
end

A_r_temp1=nan(2*N_modes,N_outputs);
for n_out=1:N_outputs
    A_r_temp2=zeros(2*N_modes,N_inputs);
    for n_input=1:N_inputs
        A_r_temp2(:,n_input)=W_V_mat\h_cols(N_t*(n_input-1)+(1:2*N_modes),n_out);  %[Maia, eqn (4.48) or (4.14)]
    end
    A_r_temp1(:,n_out)=mean(A_r_temp2,2);
    A_r_temp1(:,n_out)=A_r_temp1(:,n_out)/A_r_temp1(1,n_out);   % Eigenvector normalization
end
A_r=A_r_temp1(i_w_n_r,:);