function [psi1,freq1,psi2,freq2] = appariement(infoMODE1,infoMODE2,coupled)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function brings together the correlated modes issued from two sets of
%  data.
%
%  Synthax:
%  [psi1,freq1,psi2,freq2] = appariement(infoMODE1,infoMODE2,coupled)
%
%  Input data:
%  infoMODE1: structure containing the first set of parameters
%                infoMODE1.frequencyk = natural frequency
%                infoMODE1.etak = loss factor
%                infoMODE1.Bijk = modal constant,
%  infoMODE2: structure containing the second set of parameters
%                infoMODE2.frequencyk = natural frequency
%                infoMODE2.etak = loss factor
%                infoMODE2.Bijk = modal constant,
%  coupled: (nx2) vector giving the correspondance mode number between
%  infoMODE1 and infoMODE2 data.
%
%  Output data:
%  infoFRF: initial structure where modal parameters are added
%   modaux identifiés
%      FRF(index_FRF).infoMODE(index_mode).frequencyk = natural frequency,
%      FRF(index_FRF).infoMODE(index_mode).etak = loss factor,
%      FRF(index_FRF).infoMODE(index_mode).Bijk = modal constant.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


psi1 = [] ;
psi2 = [] ;
freq1 = [] ;
freq2 = [] ; 
psi1_temp = infoMODE1.Bijk ;
psi2_temp = infoMODE2.Bijk ;
freq1_temp = infoMODE1.frequencyk ;
freq2_temp = infoMODE2.frequencyk ;
[M,N] = size(coupled) ;
for index = 1:M
    psi1 = [psi1 psi1_temp(:,coupled(index,1))] ;
    psi2 = [psi2 psi2_temp(:,coupled(index,2))] ;
    freq1 = [freq1 freq1_temp(1,coupled(index,1))] ;
    freq2 = [freq2 freq2_temp(1,coupled(index,2))] ;
end
