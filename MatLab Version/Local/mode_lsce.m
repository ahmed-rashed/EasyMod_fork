function A_r=mode_lsce(h_cols,V_r_col,N_modes,ind2_col)
 
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Extraction of "physical modes", that's to say the modes corresponding to complex conjugate poles.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

[N_t,N_outputs]=size(h_cols); 
N_inputs=floor(N_t/N_t); %Bug is here

% Solving equation W_V_mat * A_r_temp= h [Maia, eqn (4.48) or (4.14)]
W_V_mat=zeros(2*N_modes,2*N_modes); 
for n_mode=1:2*N_modes
    W_V_mat(n_mode,:)=V_r_col(1:2*N_modes).'.^(n_mode-1); %This seems a bug
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
A_r=A_r_temp1(ind2_col,:);