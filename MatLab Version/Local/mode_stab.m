function infoMODE=mode_stab(FRF,psi_phys,fn_phys,xi_phys,RES,prec1,prec2) 
 
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Extraction of "physical modes" corresponding to those stabilized in
%  frequency using the LSCE method.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


M=size(psi_phys,1);
% Research of modes corresponding to the stabilized frequencies
Mres=size(RES,1); 
kk=0;
for ii=1:Mres
    for jj=1:length(fn_phys)
        if (abs(RES(ii,2)-xi_phys(jj))/RES(ii,2) < prec2) && (abs(RES(ii,1)-fn_phys(jj))/RES(ii,1) < prec1)
            kk=kk+1;
            Psi(:,kk)=psi_phys(:,jj);
            Fk(kk)=fn_phys(jj);
            Xik(kk)=xi_phys(jj);
        end
    end
end

psi=zeros(3*M,kk); % Zero initialisation of the modal matrix
for ind=1:kk
    infoMODE.frequencyk(1,ind)=Fk(ind);
    infoMODE.etak(ind)=2*Xik(ind);
end

for index=1:1:M
    direction=FRF(index).dir_response;
    psi(index*3-3+direction,:)=Psi(index,:);     
end
infoMODE.Bijk=psi;
