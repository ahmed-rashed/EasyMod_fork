function infoMODE=mode_stab(infoFRF,A_r_phys,f_n_r_phys,zeta_r_phys,lsce_result,prec1_percent,prec2_percent) 
 % ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Extraction of "physical modes" corresponding to those stabilized in frequency using the LSCE method.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

M=size(A_r_phys,2);
% Research of modes corresponding to the stabilized frequencies
kk=0;
for ii=1:size(lsce_result,1)
    for jj=1:length(f_n_r_phys)
        if (abs(lsce_result(ii,2)-zeta_r_phys(jj))/lsce_result(ii,2) < prec2_percent) && (abs(lsce_result(ii,1)-f_n_r_phys(jj))/lsce_result(ii,1) < prec1_percent)
            kk=kk+1;
            A_r(:,kk)=A_r_phys(jj,:);
            f_n_r(kk)=f_n_r_phys(jj);
            zeta_r(kk)=zeta_r_phys(jj);
        end
    end
end

for n=1:kk
    infoMODE.frequencyk(1,n)=f_n_r(n);
    infoMODE.etak(n)=2*zeta_r(n);
end

infoMODE.Bijk=zeros(3*M,kk);
for n=1:M
    direction=infoFRF(n).dir_response;
    infoMODE.Bijk(3*(n-1)+direction,:)=A_r(n,:);     
end
