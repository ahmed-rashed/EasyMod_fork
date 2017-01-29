function [h_cols,N_inp,N_out,N_t,D_t]=calculate_h(H_oneSided_cols,D_f,infoFRF)

[N_f_max,N_H]=size(H_oneSided_cols);
N_t=2*N_f_max-2;

% Initializa the matrix with zeros
N_out=max([infoFRF.response]);
N_inp=max([infoFRF.excitation]);
h_cols=nan(N_t,N_out);
for ind=1:N_H
    h_col=ifft_one_sided(H_oneSided_cols(:,ind))*N_t;
    i=infoFRF(ind).response;
    h_cols(:,i)=h_col;
end

f_s=N_t*D_f;
D_t=1/f_s;