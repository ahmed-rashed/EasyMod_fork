function h_col=ifft_one_sided(H_one_sided_cols)

N_f_max=size(H_one_sided_cols,1);
N=2*N_f_max-2; %This function assumes N=even. This is a valid assumption since N is usually power of 2. Thus n=N/2+1

% Note that the full H is conjugate symmetric
H_one_sided_unscaled_cols=H_one_sided_cols;
H_one_sided_unscaled_cols(2:N_f_max,:)=H_one_sided_unscaled_cols(2:N_f_max,:)/2;
h_col=ifft(H_one_sided_unscaled_cols,N,1,'symmetric');
