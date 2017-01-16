function [DELTA_FREQ,MAC,coupled] = mac(infoMODE1,infoMODE2,macinf,delfreqsup,graph)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function calculates the modal assurance criterion (MAC) of two
%  modal parameters sets and plots, if necessary, the associated chart.
%
%  Synthax:
%  [DELTA_FREQ,MAC,coupled] = mac(infoMODE1,infoMODE2,macinf,delfreqsup,graph)
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
%  macinf: minimal value of MAC for correlation coupling,
%  delfreqsup: maximum frequency gap (%) allowed between two mode shapes,
%  graph: bolean variable for visualizing MAC (1: plot, 0: no plot). 
%
%  Output data:
%  DELTA_FREQ: frequency gap (%),
%  MAC: MAC matrix,
%  coupled: (nx2) vector giving the correspondance mode number between
%  infoMODE1 and infoMODE2 data, using the input criteria.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


coupled = [] ;
compteur = 0 ;
Nexp1 = length(infoMODE1.frequencyk) ;
Nexp2 = length(infoMODE2.frequencyk) ;
psi1 = infoMODE1.Bijk ;
psi2 = infoMODE2.Bijk ;

for ii = 1:Nexp1
    for jj = 1:Nexp2
        DELTA_FREQ(ii,jj) = (abs(infoMODE1.frequencyk(ii)-infoMODE2.frequencyk(jj)))/(infoMODE1.frequencyk(ii)) ;
         MAC(ii,jj) = abs(sum(psi1(:,ii).*conj(psi2(:,jj))))^2/((sum(psi1(:,ii).*conj(psi1(:,ii))))*(sum(psi2(:,jj).*conj(psi2(:,jj))))) ;
        if (MAC(ii,jj) > macinf) && (DELTA_FREQ(ii,jj) <= delfreqsup)
            compteur = compteur+1 ;
            coupled(compteur,1) = ii ;
            coupled(compteur,2) = jj ;
        elseif (MAC(ii,jj) > macinf) && (DELTA_FREQ(ii,jj) > delfreqsup)
               compteur = compteur+1 ;
               coupled(compteur,1) = ii ;
               coupled(compteur,2) = jj ;
        end                 
    end
end

if graph == 1
    [l c] = size(MAC) ;
    MAC = [MAC zeros(l,1)] ;
    MAC = [MAC ; zeros(1,c+1)] ;
    [l c] = size(MAC) ;
    figure ;
    colormap(cool) ;
    pcolor(MAC) ;
    colorbar ;
    title('MAC matrix') ;
    xlabel('Mode shapes of analysis 2') ;
    ylabel('Mode shapes of analysis 1') ;
end
