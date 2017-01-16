function [H_estimate,frq,infoFRF] = estimate_frf(x,y)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function builds the H1 estimator from time history of input x and
//  output y. 
//
//  Synthax:
//  [H_estimate,frq,infoFRF] = estimate_frf(x,y)
//
//  Input data:
//  x: reference signal in [V],
//  y: response signal in [V].
//
//  Output data:
//  H_estimate: FRF estimated from H1 estimator,
//  freq: frequency vector,
//  infoFRF: structure containing information on FRFs 
//            infoFRF(jj).response = jjth FRF response node
//            infoFRF(jj).dir_response = jjth FRF response direction (1=X, 2=Y,
//             3=Z, 4=RotX, 5=RotY, 6=RotZ)
//            infoFRF(jj).excitation = jjth FRF excitation node
//            infoFRF(jj).dir_excitation = jjth FRF excitation direction
//            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ).
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


// Definition of parameters
prompt0 = {'Sensitivity of channel 1 sensor: ', 'Sensitivity of channel 1 sensor: ','Sampling frequency: ',...
        'Number of sections used for estimating the FRF: ','Overlap [%]: ',...
        'Excitation node: ', 'Excitation direction: ','Response node: ', 'Response direction: '};
        
dlg_title = 'FRF building' ;
answer = x_mdialog(dlg_title,prompt0,{'','','','','','','','',''}) ;
sens1 = eval(answer(1,1)) ;
sens2 = eval(answer(2,1)) ;
fs = eval(answer(3,1)) ;
section = eval(answer(4,1)) ;
Overlap = eval(answer(5,1)) ;
Overlap = Overlap/100 ;
reference_node = eval(answer(6,1)) ;
reference_dir = eval(answer(7,1)) ;
resp_node = eval(answer(8,1)) ;
resp_dir = eval(answer(9,1)) ;

// Data loading
x = x*sens1 ;
y = y*sens2 ;

// FRF building
nfft = length(x) ;
n = floor(nfft/section) ;
nfft = n*section ;
[Sxy]=pspect(n*Overlap,n,'hn',x,y) ;
[Sxx]=pspect(n*Overlap,n,'hn',x) ;
H = Sxy./Sxx ;
H_estimate = H(1:1:round(length(H)/2)+1) ;
frq = linspace(0,fs/2,length(H_estimate));
frq = frq' ;
fmin = frq(1,1) ;
Df = frq(2,1)-frq(1,1) ;
infoFRF = struct('response',resp_node,'dir_response',resp_dir,'excitation',reference_node,'dir_excitation',reference_dir) ; 

endfunction