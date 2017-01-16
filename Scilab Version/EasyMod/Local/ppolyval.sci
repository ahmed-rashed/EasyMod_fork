function [X,E] = ppolyval(A,B,n,p)

// ------------------   This file is part of EasyMod   ----------------------------
//  Internal function
//
//  This function solves the eigenvalue polynomial problem. It is the
//  modified file polyeig.m adapted to the eigenvalue problem of p-degree.
//  it solve the eigenvalue problem
//           (A0 + lambda*A1 +  ... + lambda^p*Ap)*x = 0.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


// Use the QZ algorithm on the pair of matrices
[alpha,bet,Q,Z] = spec(A,B) ;

// Extract and check the results
X = Z ;
X = X(1:n,:) ;
atol = 100*n*max(abs(alpha))*%eps ;
btol = 100*n*max(abs(bet))*%eps ;
E = zeros(n*p,1) ;
for j = 1:n*p
    if abs(alpha(j)) < atol & abs(bet(j)) < btol
       disp(' ') ;
       disp('Warning: Rank deficient generalized eigenvalue problem.') ;
       disp('Eigenvalues are not well determined.  Results may be inaccurate.') ;
    end
    E(j) = alpha(j)/bet(j) ;
end

endfunction