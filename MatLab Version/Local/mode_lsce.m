function A_r=mode_lsce(h_cols,N_inputs,V_r_col,ind_col)
 
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Extraction of "physical modes", that's to say the modes corresponding to complex conjugate poles.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

[N_t,N_outputs]=size(h_cols); 
N_modes=size(V_r_col,1)/2;

% Solving equation W_V_mat*A_r_temp=h [Maia, eqn (4.48) or (4.14)]
W_V_mat=V_r_col.'.^((1:2*N_modes).'-1);

A_r_temp1=nan(2*N_modes,N_outputs);
for n_output=1:N_outputs
    A_r_temp2=zeros(2*N_modes,N_inputs);
    for n_input=1:N_inputs
        A_r_temp2(:,n_input)=W_V_mat\h_cols(N_t*(n_input-1)+(1:2*N_modes),n_output);  %Bug here if N_inputs>1!!!! [Maia, eqn (4.48) or (4.14)]
    end
    A_r_temp1(:,n_output)=mean(A_r_temp2,2);
end
A_r=A_r_temp1(ind_col,:)./A_r_temp1(ind_col(1),:);   % Eigenvector normalization