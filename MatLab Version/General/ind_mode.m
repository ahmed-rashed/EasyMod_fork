function [I] = ind_mode(H)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function generates three mode indicators - sum of FRFs, sum of FRFs real
%  part and sum of FRFs imaginary part.
%
%  Synthax:
%  [I] = ind_mode(H)
%
%  Input data:
%  H: FRF matrix containing all the FRFs (column number = FRF number). 
%
%  Output data:
%  I: structure containing the different indicators
%        I.ISUM: sum of included FRFs,
%        I.ISRe: sum of included real part FRFs,
%        I.ISIm: sum of included imaginary part FRFs.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


[m,n] = size(H) ;
for index = 1:m ;
    ISUM(index,1) = sum(abs(H(index,:)))/(n) ;
    ISRe(index,1) = sum(real(H(index,:)))./sum(abs(H(index,:)))/(n) ;   
    ISIm(index,1) = sum(imag(H(index,:)))/(n) ;
end
I = struct('ISUM',ISUM,'ISRe',ISRe,'ISIm',ISIm) ;