function [H_estimate,freq,infoFRF]=estimate_frf(x,y)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function builds the H1 estimator from time history of input x and
%  output y. 
%
%  Synthax:
%  [H_estimate,freq,infoFRF]=estimate_frf(x,y)
%
%  Input data:
%  x: reference signal in [V],
%  y: response signal in [V].
%
%  Output data:
%  H_estimate: FRF estimated from H1 estimator,
%  freq: frequency vector,
%  infoFRF: structure containing information on FRFs 
%            infoFRF(jj).response=jjth FRF response node
%            infoFRF(jj).dir_response=jjth FRF response direction (1=X, 2=Y,
%             3=Z, 4=RotX, 5=RotY, 6=RotZ)
%            infoFRF(jj).excitation=jjth FRF excitation node
%            infoFRF(jj).dir_excitation=jjth FRF excitation direction
%            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ).
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


% Definition of parameters
prompt={'Sensitivity of channel 1 sensor: ', 'Sensitivity of channel 1 sensor: ','Sampling frequency: ',...
        'Number of sections used for estimating the FRF: ','Overlap [%]: ',...
        'Excitation node: ', 'Excitation direction: ','Response node: ', 'Response direction: '};
        
dlg_title='FRF building';
answer=inputdlg(prompt,dlg_title);
sens1=str2num(char(answer(1,1)));
sens2=str2num(char(answer(2,1)));
fs=str2num(char(answer(3,1)));
section=str2num(char(answer(4,1)));
Overlap=str2num(char(answer(5,1)));
Overlap=Overlap/100;
reference_node=str2num(char(answer(6,1)));
reference_dir=str2num(char(answer(7,1)));
resp_node=str2num(char(answer(8,1)));
resp_dir=str2num(char(answer(9,1)));

% Data loading
x=x*sens1;
y=y*sens2;

% FRF building
nfft=length(x);
n=floor(nfft/section);
nfft=n*section;
[H_estimate,freq]=TFESTIMATE(x,y,hann(n),floor(Overlap*n),n,fs);
fmin=freq(1,1);
Df=freq(2,1)-freq(1,1);
infoFRF=struct('response',resp_node,'dir_response',resp_dir,'excitation',reference_node,'dir_excitation',reference_dir); 



