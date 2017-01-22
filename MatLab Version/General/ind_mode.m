function I=ind_mode(H_cols)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function generates three mode indicators - sum of FRFs, sum of FRFs real
%  part and sum of FRFs imaginary part.
%
%  Synthax:
%  I=ind_mode(H_cols)
%
%  Input data:
%  H_cols: FRF matrix containing all the FRFs (column number=FRF number). 
%
%  Output data:
%  I: structure containing the different indicators
%        I.ISUM: sum of included abs of FRFs,
%        I.ISRe: sum of included real part of FRFs,
%        I.ISIm: sum of included imaginary part of FRFs.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


n=size(H_cols,2);
ISUM=sum(abs(H_cols),2)/n;
ISRe=sum(real(H_cols),2)/n;
ISIm=sum(imag(H_cols),2)/n;

I=struct('ISUM',ISUM,'ISRe',ISRe,'ISIm',ISIm);