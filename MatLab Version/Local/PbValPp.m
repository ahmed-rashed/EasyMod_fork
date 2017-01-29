function [L,z_col]=PbValPp(beta_cols,N_inp,p)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function solves the eigenvalue problem in a specific case.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


%  Necessary functions:
%  -----------------------------------------------------------
%  ppolyval.m


% Extraction of coefficients A1,...,A_{p} of dimension (N_inp x N_inp)
%   (A_{p} + z_col*A_{p-1} +  ... + z_col^{p-1}*A_{1}+z_col^{p})*L=0
A=zeros(N_inp,N_inp,p);
for j=1:p
   A(:,:,j)=beta_cols((j-1)*N_inp+(1:N_inp),:);
end

% Transformation for solving the problem
%   (B_{0} + z_col*B_{1} +  ... + z_col^{p}*B_{p})*L=0
B=zeros(N_inp,N_inp,p+1);
for i=1:p
   B(:,:,i)=A(:,:,p-(i-1));
end
B(:,:,p+1)=eye(N_inp);

% Building of two matrices p*N_inp x p*N_inp
%    J=[B0   0   0   0]   K=[-B1 -B2 -B3 -B4]
%      [ 0   I   0   0]     [  I   0   0   0]
%      [ 0   0   I   0]     [  0   I   0   0]
%      [ 0   0   0   I]     [  0   0   I   0]
J=eye(N_inp*p);
J(1:N_inp,1:N_inp)=B(:,:,1);
if p == 0
   K=eye(N_inp);
   p=1;
else
   K=diag(ones(N_inp*(p-1),1),-N_inp);
   for k=1:p
      K(1:N_inp,(k-1)*N_inp+1:k*N_inp)=-B(:,:,k+1);
   end
end

% Using ppolyval function
[L,z_col]=ppolyval(J,K,N_inp,p);
%Eigenvector scaling
for k=1:size(L,2)
   L(:,k)=L(:,k)/L(1,k);
end
