function [eigVec_up,eigVal_col]=PbValPp(Beta_T_transpose,N_inputs,p)
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function solves the eigenvalue problem in a specific case.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT

% (beta_{0} + z_col*beta_{1} +  ... + z_col^{p}*beta_{p})*L=0

% Building of two matrices p*N_inputs x p*N_inputs    %Maia equation (4.44)
%    J=[beta0   0   0   0]   K=[-beta1 -beta2 -beta3 -beta4]
%      [ 0   I   0   0]        [  I   0   0   0]
%      [ 0   0   I   0]        [  0   I   0   0]
%      [ 0   0   0   I]        [  0   0   I   0]

% By Ahmed Rashed
%This file seems to intentionally confuse readers of the code


J=eye(N_inputs*p);
J(1:N_inputs,1:N_inputs)=Beta_T_transpose((p-1)*N_inputs+(1:N_inputs),:);

if p == 0
   K=eye(N_inputs);
   p=1;
else
   K=diag(ones(N_inputs*(p-1),1),-N_inputs);
   for jj=1:p-1
      K(1:N_inputs,(jj-1)*N_inputs+(1:N_inputs))=-Beta_T_transpose((p-jj-1)*N_inputs+(1:N_inputs),:);
   end
   K(1:N_inputs,(p-1)*N_inputs+(1:N_inputs))=-eye(N_inputs);
end

% Calculating eigenvalues
[eigVec_full,eigVal_mat]=eig(J,K);
eigVal_col=diag(eigVal_mat);

[~,Index]=sort(abs(imag(eigVal_col)));
eigVal_col=eigVal_col(Index);
eigVec_full=eigVec_full(:,Index);

EigValues_prec=eps*N_inputs*p*max(abs(eigVal_col));
if any(abs(eigVal_col)<=EigValues_prec)
    warning('Warning: Rank deficient generalized eigenvalue problem. Eigenvalues are not well determined. Results may be inaccurate.');
end
eigVec_up=eigVec_full(1:N_inputs,:);

%Eigenvector scaling
for k=1:size(eigVec_up,2)
   eigVec_up(:,k)=eigVec_up(:,k)/eigVec_up(1,k);
end
