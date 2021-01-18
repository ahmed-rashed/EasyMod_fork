function h_mat_rows=MatSur2(h_cols,N_inputs,p)
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function builds the overdetermined matrix from
%    - the reduced matrix resulting from the analysis of the impulse responses matrix,
%    - the number of outputs and inputs (N_outputs and N_inputs),
%    - the last measurement point taken into account (Npmax),
%    - the number of equations by output (Neq),
%    - the number of measurement points (N_t),
%    - the order (p).
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

[N_t,N_outputs]=size(h_cols);
Neq=floor(N_t/(p+1)); %May be p+1=N_outputs

h_3D=zeros(Neq,(p+1)*N_inputs,N_outputs);
for n=1:Neq
    for m=1:p+1
        h_3D(n,(m-1)*N_inputs+(1:N_inputs),:)=h_cols((p+1)*n-(m-1),:);    %h matrix of [Maia, eqn. (4.35)]
    end 
end

% Overdetermined matrix GG assembly
h_mat_rows=zeros(Neq*N_outputs,(p+1)*N_inputs);
for n_output=1:N_outputs
   h_mat_rows((n_output-1)*Neq+(1:Neq),:)=h_3D(:,:,n_output); 
end