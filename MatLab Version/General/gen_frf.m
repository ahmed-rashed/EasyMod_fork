function [recep,mobil,inert]=gen_frf(M,C,K,numin,numout,freq)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function generates a FRF from mass, damping and stiffness matrices,
%  in various format (compliance, mobility, accelerance).
%
%  Synthax:
%  [recep,mobil,inert]=gen_frf(M,C,K,numin,numout,freq)
%
%  Input data:
%  M: mass matrix,
%  C: damping matrix,
%  K: stiffness matrix,
%  numin: excitation DOF,
%  numout: response DOF,
%  freq: frequency vector. 
%
%  Output data:
%  recept: receptance vector,
%  mobil: mobility vector,
%  inert: accelerance vector.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


freq=freq(:);
omega=freq*2*pi;
recep=nan(size(omega));
for i=1:length(omega)
    H_i=K-omega(i)^2*M+1i*omega(i)*C;
    H_i_inv=inv(H_i);  
    recep(i)=H_i_inv(numin,numout);
end
mobil=recep.*omega*1i;
inert=-recep.*omega.^2;

