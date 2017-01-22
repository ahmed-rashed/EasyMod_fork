function [freq_local,H_local,H_gen_local,infoMODE,line_prop]=line_fit(H,freq,fmin,fmax)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Identification based on the line-fit method (SDOF method) giving the
%  natural frequency, the loss factor and the modal constant.
%
%  Synthax :
%  [freq_local,H_local,H_gen_local,infoMODE,line_prop]=line_fit(H,freq,fmin,fmax)
%
%  Input data:
%  H: FRF vector that we want to analyse,
%  freq: frequency vector,
%  fmin and fmax: frequency range embracing the mode that we xant to
%                 analyse.
%
%  Output data:
%  freq_local: studied frequency vector extrated from the overall frequency,
%  H_local: studied FRF vector extrated from the overall FRF,
%  H_gen_local: synthetized FRF from modal parameters,
%  infoMODE: structure containing the different identified parameters
%                infoMODE.frequencyk=natural frequency
%                infoMODE.etak=loss factor
%                infoMODE.Bijk=modal constant,
%  line_prop: structure containing information about the line-fit method.
%          
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


%  Necessary functions:
%  -----------------------------------------------------------
%  estimate_delta.m
%  moindre_carre.m


problem=0;

% Data extraction
temp=find(freq>fmin);
index_low=temp(1,1)-1;
temp=find(freq>fmax);
index_high=temp(1,1)-1;

H_local=H(index_low:index_high,1);
freq_local=freq(index_low:index_high,1);
ww_local=2*pi*freq_local;

% In the case of insufficient number of samples
N_pts=length(freq_local); 
if N_pts < 5
    disp('!!!')
    disp('Insufficient number of samples in the studied frequency range')
    disp(' ')
    problem=1;
end

% Dobson's method application
for ind=1:N_pts
    [ww_delta, delta]=estimate_delta(ww_local,H_local,ww_local(ind));
    Delta(:,ind)=delta;
    ww_Delta(:,ind)=ww_delta;
end

% Best straight line finding
for index=1:N_pts
    [a b]=moindre_carre(ww_Delta(:,index).^2,real(Delta(:,index)));
    tr(index,1)=a;
    cr(index,1)=b;
    [a b]=moindre_carre(ww_Delta(:,index).^2,imag(Delta(:,index)));
    ti(index,1)=a;
ci(index,1)=b;
end
[ur dr]=moindre_carre(ww_local.^2,tr);
[ui di]=moindre_carre(ww_local.^2,ti);
p=ui/ur;
q=di/dr;

% Modal parameters calculation
loss_factor=(q-p)/(1+p*q);
wr=sqrt(dr/((p*loss_factor-1)*ur));
eig_frequency=wr/(2*pi);
ar=wr^2*(p*loss_factor-1)/((1+p^2)*dr);
br=-ar*p;
B=complex(ar,br);

% FRF synthetization from mosal parameters 
H_gen_local=B./(wr^2-ww_local.^2+j*loss_factor*wr^2);

% Data saving
if problem == 0
   line_prop=struct('ur',ur,'dr',dr,'ui',ui,'di',di,'freq_local',freq_local,'tr',tr,'ti',ti,'ww_Delta',ww_Delta,'Delta',Delta);
   infoMODE=struct('frequencyk',eig_frequency,'etak',loss_factor,'Bijk',B);
elseif problem == 1
   infoMODE=struct('frequencyk',NaN,'etak',NaN,'Bijk',0);   
end