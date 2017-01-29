function h_mat=gen_resp_impul(H_oneSided_cols,f_col,infoFRF)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Calculation of impulse responses from FRF.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

[n,N_H]=size(H_oneSided_cols);
N_t=2*n-2;

% Initializa the matrix with zeros
N_out=max([infoFRF.response]);
N_inp=max([infoFRF.excitation]);
h_mat=zeros(N_out,max(N_out,N_inp)*N_t);
for ind=1:N_H
      i=infoFRF(ind).response;
      j=infoFRF(ind).excitation;
      h_col=ifft_one_sided(H_oneSided_cols(:,ind))*N_t;
      h_mat(i,(j-1)*N_t+1:j*N_t)=h_col;
end

% Time vector is added and it is identitical for all the responses
D_f=f_col(2)-f_col(1);
f_s=N_t*D_f;
D_t=1/f_s;
t_row=(0:N_t-1)/D_t;
h_mat(N_out+1,1:N_t)=t_row;