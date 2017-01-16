function [recep,mobil,inert] = gen_frf(M,D,K,numin,numout,frq)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function generates a FRF from mass, damping and stiffness matrices,
//  in various format (compliance, mobility, accelerance).
//
//  Synthax:
//  [recep,mobil,inert] = gen_frf(M,D,K,numin,numout,frq)
//
//  Input data:
//  M: mass matrix,
//  D: damping matrix,
//  K: stiffness matrix,
//  numin: excitation DOF,
//  numout: response DOF,
//  frq: frequency vector. 
//
//  Output data:
//  recept: receptance vector,
//  mobil: mobility vector,
//  inert: accelerance vector.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


frq  =  frq(:) ;
omega = frq*2*%pi ;
tfunc1 = omega ;
for i = 1:length(omega)
    MDK = K+%i*omega(i)*D-omega(i)^2*M ;
    MDKi = inv(MDK) ;  
    tfunc1(i) = MDKi(numin,numout) ;
end
tfunc2 = tfunc1.*omega*%i ;
tfunc3 = -tfunc1.*omega.^2 ;
recep = tfunc1 ;
mobil = tfunc2 ;
inert = tfunc3 ;

endfunction