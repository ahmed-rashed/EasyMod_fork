function GG_cols=MatSur2(h_cols,N_inp,p)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function builds the overdetermined matrix from
%    - the reduced matrix resulting from the analyse of the impulse responses matrix,
%    - the number of outputs and inputs (N_out and N_inp),
%    - the last measurement point taken into account (Npmax),
%    - the number of equations by output (Neq),
%    - the number of measurement points (N_t),
%    - the order (p).
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

[N_t,N_out]=size(h_cols);
Neq=floor(N_t/(p+1)); %May be p+1=N_out
% Sub-matrices G assembly
G=zeros(Neq,(p+1)*N_inp,N_out);
for l=1:Neq
    n=(p+1)*l;
    for m=1:p+1
        G(l,(m-1)*N_inp+1:m*N_inp,:)=h_cols(n-(m-1),:);
    end 
end

% Overdetermined matrix GG assembly
GG_cols=zeros(N_out*Neq,(p+1)*N_inp);
for i=1:N_out
   GG_cols((i-1)*Neq+(1:Neq),:)=G(:,:,i); 
end