function [f_r,eta_r,A_r]=DobsonMethod(f_local_vec,Receptance_local_vec,ShowInternalDetails)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Identification based on the Dobson method (SDOF method) giving the
%  natural frequency, the loss factor and the modal constant.
%
%  Synthax :
%  [line_prop]=DobsonMethod(f_local_vec,Receptance_local_vec)
%
%  Input data:
%  Receptance_local_vec: FRF vector that we want to analyse,
%  f_local_vec: frequency vector,
%
%  Output data:
%  f_local_vec: studied frequency vector extrated from the overall frequency,
%  Receptance_local_vec: studied FRF vector extrated from the overall FRF,
%  H_gen_local: synthetized FRF from modal parameters,
%  infoMODE: structure containing the different identified parameters
%                infoMODE.frequencyk=natural frequency
%                infoMODE.etak=loss factor
%                infoMODE.Bijk=modal constant,
%  line_prop: structure containing information about the Dobson method.
%          
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


%  Necessary functions:
%  -----------------------------------------------------------
%  estimate_delta.m

w_local_vec=2*pi*f_local_vec;

% In the case of insufficient number of samples
N_pts=length(f_local_vec); 
if N_pts < 5
    error('Insufficient number of points in the studied frequency range');
end

% Dobson's method application
w_Delta_cols=nan(N_pts-1,N_pts);
Delta_cols=w_Delta_cols;
for ind=1:N_pts
    [w_Delta_cols(:,ind), Delta_cols(:,ind)]=estimate_delta(w_local_vec,Receptance_local_vec,ind);
end

% Best straight line finding
tr=nan(size(f_local_vec));
ti=tr;
for ind=1:N_pts
    p1=polyfit(w_Delta_cols(:,ind).^2,real(Delta_cols(:,ind)),1);
    tr(ind)=p1(1);
    p2=polyfit(w_Delta_cols(:,ind).^2,imag(Delta_cols(:,ind)),1);
    ti(ind)=p2(1);
end
p1=polyfit(w_local_vec.^2,tr,1);
ur=p1(1);
dr=p1(2);
p2=polyfit(w_local_vec.^2,ti,1);
ui=p2(1);
di=p2(2);

% Modal parameters calculation
p=ui/ur;
q=di/dr;
eta_r=(q-p)/(1+p*q);
w_r=sqrt(dr/((p*eta_r-1)*ur));
f_r=w_r/2/pi;
a_r=w_r^2*(p*eta_r-1)/((1+p^2)*dr);
b_r=-a_r*p;
A_r=complex(a_r,b_r);

if ShowInternalDetails
    % Data saving
    line_prop=struct('ur',ur,'dr',dr,'ui',ui,'di',di,'tr',tr,'ti',ti,'w_Delta',w_Delta_cols,'Delta',Delta_cols);
    plot_line_fit(f_local_vec,Receptance_local_vec,f_r,eta_r,A_r,line_prop);
end
