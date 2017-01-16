function plot_circle_fit(freq_local,H_local,H_gen_local,infoMODE,circ_prop)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  Graphical representation of the results provided by the function 'circle_fit'.
//
//  Synthax:
//  plot_circ_prop_fit(freq_local,H_local,H_gen_local,infoMODE,circ_prop)
//
//  Input data:
//  freq_local: studied frequency vector extrated from the overall frequency,
//  H_local: studied FRF vector extrated from the overall FRF,
//  H_gen_local: synthetized FRF from modal parameters,
//  infoMODE: structure containing the difreqerent identified parameters
//                infoMODE.frequencyk = natural frequency
//                infoMODE.etak = loss factor
//                infoMODE.Bijk = modal constant,
//  circ_prop: structure containing information about Nyquist's circle
//          circ_prop.x0 = abscissa of circle center
//          circ_prop.y0 = ordinate of circle center
//          circ_prop.R0 = circ_prop radius
//          circ_prop.D  = damping evolution.
//
//  These data are directly obtained with the function 'circle_fit'.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


if isnan(circ_prop.R0)
    figure
else
if isempty(infoMODE.frequencyk) == %F
x = real(H_local) ;
y = imag(H_local) ;

// Nyquist's circle characteristics
x0 = circ_prop.x0 ;
y0 = circ_prop.y0 ;
R0 = circ_prop.R0 ;
theta_ref = circ_prop.theta_ref ;
eig_theta = circ_prop.eig_theta ;
theta = circ_prop.theta ;
N_pts = length(freq_local) ;

// Modal parameters
eig_frequency = infoMODE.frequencyk ;
loss_factor = infoMODE.etak ;
Z = circ_prop.D ;
B = infoMODE.Bijk ;

figure
// Nyquist's circle visualization
subplot(2,2,1)
plot(x,y,'*','Linewidth',1) ;
set(gca(),"auto_clear","off") ;
t = 0:0.01:2*%pi ;
plot(R0*sin(t)+x0,R0*cos(t)+y0,'color','r','Linewidth',1) ;
set(gca(),"auto_clear","off") ;

// Natural frequency localisation on the curve
ref = 2*%pi-(eig_theta+theta_ref) ;
xa = R0*cos(ref)+x0 ;
ya = R0*sin(ref)+y0 ;
set(gca(),"auto_clear","off") ;
plot([xa x0],[ya y0],'color','green','Linewidth',1) ;
xlabel('Real part FRF') ;
ylabel('Imaginary part FRF') ;
title('Nyquist curve') ;

// FRF building
subplot(2,2,2) ;
plot(freq_local,20*log10(abs(H_gen_local)),'color','r','Linewidth',1) ;
set(gca(),"auto_clear","off") ;
plot(freq_local,20*log10(abs(H_local)),'Linewidth',1) ;
set(gca(),"grid",[1 1]) ;
xlabel('Frequency  [Hz]') ;
ylabel('Gain  [dB]') ;
title('Bode curve') ;
legend('generated FRF','measured FRF',4) ;

// Damping matrix representation
subplot(2,2,4) ;
[M,N] =size(Z) ;
finc = (freq_local(2,1)-freq_local(1,1)) ;
temp = find(theta>eig_theta) ;
k1 = temp(1,1)-1 ;
k2 = temp(1,1) ;
wb = freq_local(k1:-1:1)*2*%pi ; // half-power low frequency
wa = freq_local(k2:N_pts)*2*%pi ; // half-power high frequency
fa = wa/(2*%pi) ;
fb = wb/(2*%pi) ;
fb = fb(1:M) ;
fa = fa(1:N) ;
[X,Y] = meshgrid(fa,fb) ;
if length(fa) == 1 
    plot(fb,Z,'g') ;
    axis([fb(length(fb)) fb(1) Z(1) Z(length(Z))]) ;
    xlabel('After resonance') ;
    ylabel('Loss factor') ;
    title('Loss factor evolution') ;
else
    if length(fb) == 1
        plot(fa,Z,'g') ;
        axis([fa(1,1) fa(length(fa),1) Z(1) Z(length(Z))]) ;
        xlabel('Before resonance') ;
        ylabel('Loss factor') ;
        title('Loss factor evolution') ;
    else
        surf(X,Y,Z) ;
        a = get("current_axes") ;
        a.data_bounds = [fa(1,1),fb(length(fb),1);fa(length(fa),1),fb(1,1)] ;
        xlabel('After resonance') ;
        ylabel('Before resonance') ;
        zlabel('$\eta$') ;
        title('Loss factor evolution') ;
    end
end

// Results displaying
subplot(2,2,3) ;
xstring(0.1,0.7,'Natural frequency [Hz]:') ;
xstring(0.1,0.6,'Damping constant [%]:') ;
xstring(0.1,0.5,'Modal Const MAG:') ;
xstring(0.1,0.4,'Modal Const phase [°]:') ;
Bmod = abs(B) ;
phi = atan(imag(B),real(B))*360/(2*%pi) ;
xstring(0.7,0.7,string(eig_frequency)) ;
xstring(0.7,0.6,string(loss_factor/2*100)) ;
xstring(0.7,0.5,string(Bmod)) ;
xstring(0.7,0.4,string(phi)) ;

elseif isempty(infoMODE.frequencyk) == %T
    figure
end
end

endfunction