function [L,z] = PbValPp(x,Ni,p)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function solves the eigenvalue problem in a specific case.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


%  Necessary functions:
%  -----------------------------------------------------------
%  ppolyval.m


% Extraction of coefficients A1,...,Ap of dimension (Ni x Ni)
%   (Ap + z*Ap-1 +  ... + z^p-1*A1+z^p)*L = 0
for j=1:p
   A(:,:,j) = x((j-1)*Ni+1:j*Ni,:) ;
end

% Transformation for solving the problem
%   (B0 + z*B1 +  ... + z^p*Bp)*L = 0
for i=1:p
   B(:,:,1+(i-1)) = A(:,:,p-(i-1)) ;
end
B(:,:,p+1) = eye(Ni) ;

% Building of two matrices p*Ni x p*Ni
%    J = [B0   0   0   0]   K = [-B1 -B2 -B3 -B4]
%        [ 0   I   0   0]       [  I   0   0   0]
%        [ 0   0   I   0]       [  0   I   0   0]
%        [ 0   0   0   I]       [  0   0   I   0]
J = eye(Ni*p,Ni*p) ;
J(1:Ni,1:Ni) = B(:,:,1) ;
if p == 0
   K = eye(Ni,Ni) ;
   p = 1 ;
else
   K = diag(ones(Ni*(p-1),1),-Ni) ;
   for k = 1:p
      K(1:Ni,(k-1)*Ni+1:k*Ni) = -B(:,:,k+1) ;
   end
end

% Using ppolyval function
[L,z] = ppolyval(J,K,Ni,p) ;
lz = length(z) ;
for k = 1:lz
   L(:,k) = L(:,k)/L(1,k) ;
end
