function [X_up,eigVal_col]=ppolyval(A,B,N_inp,p)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function solves the eigenvalue polynomial problem. It is the
%  modified file polyeig.m adapted to the eigenvalue problem of p-degree.
%  it solve the eigenvalue problem
%           (A_{0} + lambda*A_{1} +  ... + lambda^{p}*A_{p})*x=0.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


% Use the QZ algorithm on the pair of matrices
[alpha,beta,~,~,X]=qz(A,B);

% Extract and check the results
X_up=X(1:N_inp,:);
alpha_diag=diag(alpha);
beta_diag=diag(beta);
atol=100*N_inp*max(abs(alpha_diag))*eps;
btol=100*N_inp*max(abs(beta_diag))*eps;
eigVal_col=zeros(N_inp*p,1);
for j=1:N_inp*p
   if abs(alpha_diag(j)) < atol && abs(beta_diag(j)) < btol
     warning('Warning: Rank deficient generalized eigenvalue problem. Eigenvalues are not well determined.  Results may be inaccurate.');
   end
   eigVal_col(j)=alpha_diag(j)/beta_diag(j);
end
