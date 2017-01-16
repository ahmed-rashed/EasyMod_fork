function modegauss(infoMODE)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function displays the mode shapes in the Gaussian plane in order to
//  verify the collinearity property.
//
//  Synthax:
//  modegauss(infoMODE)
//
//  Input data:
//  infoMODE: structure containing the set of modal parameters
//                infoMODE.frequencyk = natural frequency
//                infoMODE.etak = loss factor
//                infoMODE.Bijk = modal constant.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS



function compass(U,V)
//COMPASS Compass plot.
//   COMPASS(U,V) draws a graph that displays the vectors with
//   components (U,V) as arrows emanating from the origin.
xx = [0 1 .8 1 .8].';
yy = [0 0 .08 0 -.08].';
arrow = xx + yy.*sqrt(-1);
U = U(:);
V = V(:);
z = (U + %i*V).';
a = arrow * z;
th = atan(imag(a),real(a));
r = sqrt(real(a).^2 + imag(a).^2);
h = polarplot(th,r);
endfunction


psi = infoMODE.Bijk ;
fk = infoMODE.frequencyk ;
xik = infoMODE.etak ;
[MM NN] = size(psi) ;
kk = 0 ;
figure ;
for jj = 1:NN ;
    kk = kk+1 ;
    if pmodulo(kk,5) == 0 ;
            kk = 1 ;
    end
    for ii = 1:MM ;
        x(ii) = real(psi(ii,jj)) ;
        y(ii) = imag(psi(ii,jj)) ;     
    end
     subplot(2,2,kk) ;
     compass(x,y) ;
     xlabel('Real part') ;
     ylabel('Imaginary part') ;
     nomode = string(jj) ;
     name = ['mode shape:' nomode] ;
     title(name) ;
end

endfunction