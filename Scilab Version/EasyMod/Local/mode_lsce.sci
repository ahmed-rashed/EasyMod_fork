function [psi_phys,fn_phys,xi_phys] = mode_lsce(H,deltaT,Z,Nmode,Nt) 
 
// ------------------   This file is part of EasyMod   ----------------------------
//  Internal function
//
//  Extraction of "physical modes", that's to say the modes corresponding to
//  complex conjugate poles.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


[MM NN] = size(Z) ;  
for ind = 1:MM
    if  Z(ind,Nmode) == 0
        break
    end
        z(ind,1) = Z(ind,Nmode) ; 
end
lambda = log(z)./deltaT ;
wd = imag(lambda) ; 
delta = real(lambda) ; 
wn = sqrt(wd.^2+delta.^2) ; 
fn = wn/(2*%pi) ; 
xi = -(delta./wn) ; 
[fn_ord,ind_ord] = gsort(fn,'g','i') ;
for index = 1:length(fn_ord)
    z_ord(index) = z(ind_ord(index)) ; 
    xi_ord(index) = xi(ind_ord(index)) ; 
end

// Solving equation Vr=h
VV = zeros(2*Nmode,2*Nmode) ; 
for kk = 1:1:2*Nmode ; 
    for jj = 1:1:2*Nmode ; 
        VV(kk,jj) = z_ord(jj).^(kk-1) ; 
    end
end

[MM,LL] = size(H) ; 
Ni = LL/Nt ; 
No = MM ; 
Phi = [] ; 

if Nt<2*Nmode
    error('Number of samples insufficient for the number of modes under consideration!')
end

for ind = 1:1:No
    rr = zeros(Ni,2*Nmode) ;    
    for ind2 = 1:Ni 
        hh = H(ind,Nt*(ind2-1)+1:1:Nt*(ind2-1)+2*Nmode) ; 
        hh = hh.' ; 
        temp = VV\hh ; 
        rr(ind2,:) = temp.' ;         
    end
    if Ni == 1
        Phi_ord(ind,:) = rr ;
    else
        Phi_ord(ind,:) = sum(rr,1)/Ni ;
    end
    
end

// Extracting the physical modes (frequencies which appear two times)
jj = 0 ; 
phi = [] ; 
for index = 1:length(fn_ord)-1
    if abs(fn_ord(index)-fn_ord(index+1)) < 1e-8 ; 
        jj = jj+1 ; 
        psi_phys(:,jj) = Phi_ord(:,index) ;
        fn_phys(jj) = fn_ord(index) ; 
        xi_phys(jj) = xi_ord(index) ; 
    end
end
[MM,NN] = size(psi_phys) ; 

// Eigenvector normalization
for ind = 1:NN
    if psi_phys(1,ind) == 0 then
        psi_phys(1,ind) = 1e-6 ;
    end
    psi_phys(:,ind) = psi_phys(:,ind)/psi_phys(1,ind) ; 
end

endfunction