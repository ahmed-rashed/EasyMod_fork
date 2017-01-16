% The best way to begin
clear all
close all
clc
% Model parameters:
% mass matrix
M = [1 0 0 ; 0 1 0 ; 0 0 1] ;
% damping matrix
C = [40 0 0 ; 0 40 0 ; 0 0 40] ;
% stiffness matrix
K = [237315 -161000 0 ; -161000 398315 -161000 ; 0 -161000 398315] ;
% FRF storage
Df = 200/400;

freq = [Df:Df:200];
[receptance,mobilite,inertance] = gen_frf(M,C,K,1,1,freq) ;
unv58write(inertance,1,3,1,3,0,Df,'3DL_H11.unv') ;
[receptance,mobilite,inertance]=gen_frf(M,C,K,1,2,freq) ;
unv58write(inertance,1,3,2,3,0,Df,'3DL_H21.unv') ;
[receptance,mobilite,inertance]=gen_frf(M,C,K,1,3,freq) ;
unv58write(inertance,1,3,3,3,0,Df,'3DL_H31.unv') ;

% FRF loading
[H11,freq,infoFRF(1)] = unv58read('3DL_H11.unv') ;
[H21,freq,infoFRF(2)] = unv58read('3DL_H21.unv') ;
[H31,freq,infoFRF(3)] = unv58read('3DL_H31.unv') ;
H = [H11,H21,H31] ;
infoFRF2 = infoFRF ;
infoFRF3 = infoFRF ;

% Visualization on various layouts
figure
subplot(4,3,1)
plot(freq,20*log10(abs(H11)))
xlabel('Frequency [Hz]'), ylabel('H_{11} [dB]')
subplot(4,3,2)
plot(freq,20*log10(abs(H21)))
xlabel('Frequency [Hz]'), ylabel('H_{21} [dB]')
subplot(4,3,3)
plot(freq,20*log10(abs(H31)))
xlabel('Frequency [Hz]'), ylabel('H_{31} [dB]')
subplot(4,3,4)
plot(real(H11),imag(H11))
xlabel('H_{11} [real]'), ylabel('H_{11} [imag]')
subplot(4,3,5)
plot(real(H21),imag(H21))
xlabel('H_{21} [real]'), ylabel('H_{21} [imag]')
subplot(4,3,6)
plot(real(H31),imag(H31))
xlabel('H_{31} [real]'), ylabel('H_{31} [imag]')
subplot(4,3,7)
plot(freq,real(H11))

xlabel('Frequency [Hz]'), ylabel('H_{11} [real]')
subplot(4,3,8)
plot(freq,real(H21))
xlabel('Frequency [Hz]'), ylabel('H_{21} [real]')
subplot(4,3,9)
plot(freq,real(H31))
xlabel('Frequency [Hz]'), ylabel('H_{31} [real]')
subplot(4,3,10)
plot(freq,imag(H11))
xlabel('Frequency [Hz]'), ylabel('H_{11} [imag]')
subplot(4,3,11)
plot(freq,imag(H21))
xlabel('Frequency [Hz]'), ylabel('H_{21} [imag]')
subplot(4,3,12)
plot(freq,imag(H31))
xlabel('Frequency [Hz]'), ylabel('H_{31} [imag]')
pause

[indicators] = ind_mode(H) ;
[cursor1,cursor2] = plot_ind_mode(indicators,freq) ;
figure
subplot(3,1,1)
plot(freq,20*log10(abs(indicators.ISUM)))
xlabel('Frequency [Hz]')
ylabel('Indicator I_{SUM} [dB]')
subplot(3,1,2)
plot(freq,20*log10(abs(indicators.ISRe)))
xlabel('Frequency [Hz]')
ylabel('Indicator I_{S,Re} [dB]')
subplot(3,1,3)
plot(freq,20*log10(abs(indicators.ISIm)))
xlabel('Frequency [Hz]')
ylabel('Indicator I_{S,Im} [dB]')

% Circle-fit
N = size(H,2) ;
bornes_min = [40 75 110] ;
bornes_max = [60 100 130] ;
Nbr_mode = length(bornes_min) ;
for i=1:N
for j=1:Nbr_mode
[freq_local,H_local,H_gen_local,infoMODE,circ_prop] = circle_fit(H(:,i),freq,bornes_min(j),bornes_max(j));
plot_circle_fit(freq_local,H_local,H_gen_local,infoMODE,circ_prop);
infoFRF = add_data...
(i,infoMODE.frequencyk,infoMODE.etak,infoMODE.Bijk,infoFRF,j);
end
end
% Results saving
infoMODE1 = save_result_modal(infoFRF);
unv55write(infoMODE1,'3DL_circle_fit.unv',1);

% Line-fit
N=size(H,2);
bornes_min = [40 75 110];
bornes_max = [60 100 130];
Nbr_mode = length(bornes_min);
for i=1:N
for j=1:Nbr_mode
[freq_local,H_local,H_gen_local,infoMODE,circ_prop] = ...
line_fit(H(:,i),freq,bornes_min(j),bornes_max(j));
plot_line_fit(freq_local,H_local,H_gen_local,infoMODE,circ_prop);
infoFRF2 = add_data...
(i,infoMODE.frequencyk,infoMODE.etak,infoMODE.Bijk,infoFRF2,j);
end
end
% Results saving
infoMODE2 = save_result_modal(infoFRF2);
unv55write(infoMODE2,'3DL_line_fit.unv',1);

% Least-square complex exponential
[RES,infoFRF3,infoMODE3]=lsce(H,freq,infoFRF3);

% Results saving
unv55write(infoMODE3,'3DL_LSCE.unv',1)