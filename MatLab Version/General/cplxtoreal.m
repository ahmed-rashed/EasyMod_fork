function [MODE_REAL] = cplxtoreal(MODE_CPLX)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function transforms complex mode shapes to real mode shapes.
%
%  Synthax:
%  [MODE_REAL] = cplxtoreal(MODE_CPLX)
%
%  Input data:
%  MODE_CPLX: structure containing the different identified parameters
%                MODE_CPLX.frequencyk = natural frequency
%                MODE_CPLX.etak = loss factor
%                MODE_CPLX.Bijk = modal constant (complex).
%
%  Output data:
%  MODE_REAL: structure containing the different identified parameters
%                MODE_REAL.frequencyk = natural frequency
%                MODE_REAL.etak = loss factor
%                MODE_REAL.Bijk = modal constant (real).
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


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
    [P,S] = polyfit(x,y,1) ;
    coeff(jj) = P(1) ;
    
end
THETA = atan(coeff) ;
for ind = 1:NN
    psi_real(:,ind) = real(psi_cplx(:,ind)*exp(-1i*THETA(ind))) ;
end
MODE_REAL.frequencyk = fk ;
MODE_REAL.etak = xik ;
MODE_REAL.Bijk = psi_real ;

