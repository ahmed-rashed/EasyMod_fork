function plot_line_fit(f_local_vec,Receptance_local_vec,f_r,eta_r,A_r,line_prop)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Graphical representation of the results provided by the function 'line_fit'.
%
%  Synthax:
%  plot_line_fit(f_local_vec,Receptance_local_vec,f_r,eta_r,A_r,line_prop)
% 
%  Input data:
%  Receptance_local_vec: studied FRF vector extrated from the overall FRF,
%    f_r=natural frequency
%    eta_r=loss factor
%    A_r=modal constant,
%  line_prop: structure containing information about the line-fit method.
%          
%  These data are directly obtained with the function 'line_fit'.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


%  Necessary functions: 
%  -----------------------------------------------------------
%  err_fit_circle.m

if isempty(f_r)
    error('Empty modal parameters')
end

Receptance_local_calculated_vec=A_r./complex((2*pi*f_r)^2-(2*pi*f_local_vec).^2,eta_r*(2*pi*f_r)^2);

% Determination of Nyquist's circle using least squares method
x=real(Receptance_local_vec);
y=imag(Receptance_local_vec);
[x0,y0,R0]=err_fit_circle(x,y);

% Modal properties
eig_frequency=f_r;
loss_factor=eta_r;
B=A_r;

% Lines properties
Delta =line_prop.Delta;
tr=line_prop.tr;
ti=line_prop.ti;
ur=line_prop.ur;
dr=line_prop.dr;
ui=line_prop.ui;
di=line_prop.di;

figure
% FRF building
subplot(2,2,1);
plot(f_local_vec,20*log10(abs(Receptance_local_calculated_vec)),'color','r','Linewidth',1);
hold on;
plot(f_local_vec,20*log10(abs(Receptance_local_vec)),'Linewidth',1);
grid on;
xlabel('Frequency  [Hz]');
ylabel('Gain  [dB]');
title('Bode curve');
legend('generated FRF','measured FRF','Location','SouthEast');

% Nyquist's circle visualization
subplot(2,2,2);
plot(x,y,'*','Linewidth',1);
hold on;
t=0:0.01:2*pi;
plot(R0*sin(t)+x0,R0*cos(t)+y0,'color','r','Linewidth',1);
xlabel('Real part FRF');
ylabel('Imaginary part FRF');
title('Nyquist curve');

% Various lines visualization
subplot(4,4,11);
plot(f_local_vec.^2,real(Delta),'color','b');
ylabel('Delta');
title('Real part');
xlim([min(f_local_vec.^2) max(f_local_vec.^2)]);
subplot(4,4,12);
plot(f_local_vec.^2,imag(Delta),'color','b');
title('Imaginary part');
xlim([min(f_local_vec.^2) max(f_local_vec.^2)]);
subplot(4,4,15);
plot(f_local_vec.^2,tr,'*','color','green','MarkerSize',3);
hold on;
plot(f_local_vec.^2,ur*(2*pi*f_local_vec).^2+dr,'color','r','LineWidth',1);
ylabel('Slope');
xlabel('Square Frequency  [Hz^2]');
xlim([min(f_local_vec.^2) max(f_local_vec.^2)]);
subplot(4,4,16);
plot(f_local_vec.^2,ti,'*','color','green','MarkerSize',3);
hold on;
plot(f_local_vec.^2,ui*(2*pi*f_local_vec).^2+di,'color','r','LineWidth',1);
xlabel('Square Frequency  [Hz^2]');
xlim([min(f_local_vec.^2) max(f_local_vec.^2)]);

% Results displaying
set(uicontrol,'style','text','FontSize',10,'position',[35 110 150 20],'string','Natural frequency [Hz]:');
set(uicontrol,'style','text','FontSize',10,'position',[35 90 150 20],'string','Damping constant [%]:');
set(uicontrol,'style','text','FontSize',10,'position',[35 70 150 20],'string','Modal Const MAG:');
set(uicontrol,'style','text','FontSize',10,'position',[35 50 150 20],'string','Modal Const phase [�]:');
but1=uicontrol('style','text','position',[190 110 95 20]);
but2=uicontrol('style','text','position',[190 90 95 20]);
but3=uicontrol('style','text','position',[190 70 95 20]);
but4=uicontrol('style','text','position',[190 50 95 20]);
Bmod=abs(B);
phi=atan2(imag(B),real(B))*360/(2*pi);
set(but1,'FontSize',10,'string',eig_frequency);
set(but2,'FontSize',10,'string',loss_factor/2*100);
set(but3,'FontSize',10,'string',Bmod);
set(but4,'FontSize',10,'string',phi);