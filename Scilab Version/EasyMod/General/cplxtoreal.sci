function [MODE_REAL] = cplxtoreal(MODE_CPLX)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function transforms complex mode shapes to real mode shapes.
//
//  Synthax:
//  [MODE_REAL] = cplxtoreal(MODE_CPLX)
//
//  Input data:
//  MODE_CPLX: structure containing the different identified parameters
//                MODE_CPLX.frequencyk = natural frequency
//                MODE_CPLX.etak = loss factor
//                MODE_CPLX.Bijk = modal constant (complex).
//
//  Output data:
//  MODE_REAL: structure containing the different identified parameters
//                MODE_REAL.frequencyk = natural frequency
//                MODE_REAL.etak = loss factor
//                MODE_REAL.Bijk = modal constant (real).
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


function p = polyfit(x,y,n)
//POLYFIT Fit polynomial to data.
//   POLYFIT(X,Y,N) finds the coefficients of a polynomial P(X) of
//   degree N that fits the data, P(X(I))~=Y(I), in a least-squares sense.
//
//   For a cubic P(X), the output vector p=[p1 p2 p3 p4] satisfies
//                3        2
//     P(X) = p1 x  + p2 x  + p3 x  + p4 .

if length(x) ~= length(y)
    error('X and Y vectors must be the same size.')
end

x = x(:);
y = y(:);

// Construct Vandermonde matrix.
V = ones(length(x),n+1);
for j = n:-1:1
    V(:,j) = x.*V(:,j+1);
end

// Solve least squares problem.
[Q,R] = qr(V);
QTy = Q'*y;
p = R(1:n+1,1:n+1)\QTy(1:n+1); // Same as p = V\y;
p = p.';          // Polynomial coefficients are row vectors by convention.
endfunction

psi_cplx = MODE_CPLX.Bijk ;
fk = MODE_CPLX.frequencyk ;
xik = MODE_CPLX.etak ;
[MM NN] = size(psi_cplx) ;
for jj = 1:NN ;
    for ii = 1:MM ;
        X(ii,jj) = real(psi_cplx(ii,jj)) ;
        Y(ii,jj) = imag(psi_cplx(ii,jj)) ;
    end  
    x = X(:,jj) ;
    y = Y(:,jj) ;
    P = polyfit(x,y,1) ;
    coef(jj) = P(1) ;
    
end
THETA = atan(coef) ;
for ind = 1:NN
    psi_real(:,ind) = real(psi_cplx(:,ind)*exp(-%i*THETA(ind))) ;
end
MODE_REAL.frequencyk = fk ;
MODE_REAL.etak = xik ;
MODE_REAL.Bijk = psi_real ;

endfunction