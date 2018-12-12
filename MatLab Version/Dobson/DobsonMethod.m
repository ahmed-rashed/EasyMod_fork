function [f_r,eta_r,A_r]=DobsonMethod(f_local_vec,Receptance_local_vec,ShowInternalDetails,label_str)

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
Delta_cols=nan(N_pts-1,N_pts);
w_Delta_cols=nan(N_pts-1,N_pts);
for n=1:N_pts
    [w_Delta_cols(:,n), Delta_cols(:,n)]=estimate_Delta(w_local_vec,Receptance_local_vec,n);
end

% Best straight line finding
tr=nan(size(f_local_vec));
ti=nan(size(f_local_vec));
for n=1:N_pts
    p1=polyfit(w_Delta_cols(:,n).^2,real(Delta_cols(:,n)),1);
    tr(n)=p1(1);
    p2=polyfit(w_Delta_cols(:,n).^2,imag(Delta_cols(:,n)),1);
    ti(n)=p2(1);
end
p3=polyfit(w_local_vec.^2,tr,1);
ur=p3(1);
dr=p3(2);
p4=polyfit(w_local_vec.^2,ti,1);
ui=p4(1);
di=p4(2);

% Modal parameters calculation
p=ui/ur;
q=di/dr;
eta_r=(q-p)/(1+p*q);
w_r=sqrt(dr/((p*eta_r-1)*ur));

a_r=w_r^2*(p*eta_r-1)/((1+p^2)*dr);
b_r=-a_r*p;
A_r=complex(a_r,b_r);

f_r=w_r/2/pi;

if ShowInternalDetails
    % Data saving
    line_prop=struct('ur',ur,'dr',dr,'ui',ui,'di',di,'tr',tr,'ti',ti,'w_Delta',w_Delta_cols,'Delta',Delta_cols);
    DobsonInternalDetails(f_local_vec,Receptance_local_vec,f_r,eta_r,A_r,line_prop,label_str);
end