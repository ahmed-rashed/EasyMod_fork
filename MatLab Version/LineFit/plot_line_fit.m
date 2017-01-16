function plot_line_fit(freq_local,H_local,H_gen_local,infoMODE,line_prop)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Graphical representation of the results provided by the function 'line_fit'.
%
%  Synthax:
%  plot_line_fit(freq_local,H_local,H_gen_local,infoMODE,line_prop)
% 
%  Input data:
%  freq_local: studied frequency vector extrated from the overall frequency,
%  H_local: studied FRF vector extrated from the overall FRF,
%  H_gen_local: synthetized FRF from modal parameters,
%  infoMODE: structure containing the difreqerent identified parameters
%                infoMODE.frequencyk = natural frequency
%                infoMODE.etak = loss factor
%                infoMODE.Bijk = modal constant,
%  line_prop: structure containing information about the line-fit method.
%          
%  These data are directly obtained with the function 'line_fit'.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


%  Necessary functions: 
%  -----------------------------------------------------------
%  err_fit_circle.m


if isempty(infoMODE.frequencyk) == 0

% Determination of Nyquist's circle using least squares method    
x = real(H_local) ;
y = imag(H_local) ;
[x0,y0,R0] = err_fit_circle(x,y) ;

% Modal properties
freq_local = line_prop.freq_local ;
ww_local = 2*pi*freq_local ;
ff_Delta = 2*pi*line_prop.ww_Delta ;
eig_frequency = infoMODE.frequencyk ;
loss_factor = infoMODE.etak ;
B = infoMODE.Bijk ;

% Lines properties
Delta =line_prop.Delta ;
tr = line_prop.tr ;
ti = line_prop.ti ;
ur = line_prop.ur ;
dr = line_prop.dr ;
ui = line_prop.ui ;
di = line_prop.di ;

figure
% FRF building
subplot(2,2,1) ;
plot(freq_local,20*log10(abs(H_gen_local)),'color','r','Linewidth',1) ;
hold on ;
plot(freq_local,20*log10(abs(H_local)),'Linewidth',1) ;
grid on ;
xlabel('Frequency  [Hz]') ;
ylabel('Gain  [dB]') ;
title('Bode curve') ;
legend('generated FRF','measured FRF',4) ;

% Nyquist's circle visualization
subplot(2,2,2) ;
plot(x,y,'*','Linewidth',1) ;
hold on ;
t = 0:0.01:2*pi ;
plot(R0*sin(t)+x0,R0*cos(t)+y0,'color','r','Linewidth',1) ;
xlabel('Real part FRF') ;
ylabel('Imaginary part FRF') ;
title('Nyquist curve') ;

% Various lines visualization
subplot(4,4,11) ;
plot(freq_local.^2,real(Delta),'color','b') ;
ylabel('Delta') ;
title('Real part') ;
xlim([min(freq_local.^2) max(freq_local.^2)]) ;
subplot(4,4,12) ;
plot(freq_local.^2,imag(Delta),'color','b') ;
title('Imaginary part') ;
xlim([min(freq_local.^2) max(freq_local.^2)]) ;
subplot(4,4,15) ;
plot(freq_local.^2,tr,'*','color','green','MarkerSize',3) ;
hold on ;
plot(freq_local.^2,ur*(ww_local.^2)+dr,'color','r','LineWidth',1) ;
ylabel('Slope') ;
xlabel('Square Frequency  [Hz^2]') ;
xlim([min(freq_local.^2) max(freq_local.^2)]) ;
subplot(4,4,16) ;
plot(freq_local.^2,ti,'*','color','green','MarkerSize',3) ;
hold on ;
plot(freq_local.^2,ui*(ww_local.^2)+di,'color','r','LineWidth',1) ;
xlabel('Square Frequency  [Hz^2]') ;
xlim([min(freq_local.^2) max(freq_local.^2)]) ;

% Results displaying
set(uicontrol,'style','text','FontSize',10,'position',[35 110 150 20],'string','Natural frequency [Hz]:') ;
set(uicontrol,'style','text','FontSize',10,'position',[35 90 150 20],'string','Damping constant [%]:') ;
set(uicontrol,'style','text','FontSize',10,'position',[35 70 150 20],'string','Modal Const MAG:') ;
set(uicontrol,'style','text','FontSize',10,'position',[35 50 150 20],'string','Modal Const phase [°]:') ;
but1=uicontrol('style','text','position',[190 110 95 20]) ;
but2=uicontrol('style','text','position',[190 90 95 20]) ;
but3=uicontrol('style','text','position',[190 70 95 20]) ;
but4=uicontrol('style','text','position',[190 50 95 20]) ;
Bmod = abs(B) ;
phi = atan2(imag(B),real(B))*360/(2*pi) ;
set(but1,'FontSize',10,'string',eig_frequency) ;
set(but2,'FontSize',10,'string',loss_factor/2*100) ;
set(but3,'FontSize',10,'string',Bmod) ;
set(but4,'FontSize',10,'string',phi) ;

elseif isempty(infoMODE.frequencyk) == 1
    figure
end

