function plot_ind_mode(H_cols,freq)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function displays the three mode indicators - sum of FRFs absolute, sum of FRFs
%  real part and sum of FRFs imaginary part.
%
%  Synthax:
%  plot_ind_mode(H_cols,freq)
%
%  Input data:
%  freq: frequency vector,
%  H_cols: FRF matrix containing all the FRFs (column number=FRF number). 
%
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

n=size(H_cols,2);
ISUM=sum(abs(H_cols),2)/n;
ISRe=sum(real(H_cols),2)/n;
ISIm=sum(imag(H_cols),2)/n;

fig=figure;
subplot(3,1,1)
semilogy(freq,ISUM);
set(gca,'XTickLabel',[])
ylabel('$\sum\left|H\right|$','interpreter', 'latex');

subplot(3,1,2)
semilogy(freq,abs(ISRe))
set(gca,'XTickLabel',[])
ylabel('$\sum\Re\left(H\right)$','interpreter', 'latex')

subplot(3,1,3)
semilogy(freq,abs(ISIm))
xlabel('$f$ (Hz)', 'interpreter', 'latex')
ylabel('$\sum\Im\left(H\right)$','interpreter', 'latex')

for ii=1:3
    subplot(3,1,ii)
    grid on
end

global cursor1 cursor2;

cursor1=uicontrol('style','text','position',[10 0 60 15],'BackgroundColor',[1,1,1]);
cursor2=uicontrol('style','text','position',[80 0 60 15],'BackgroundColor',[1,1,1]);
set(uicontrol,'style','text','position',[10 15 60 15],'string','f (Hz)','BackgroundColor',[1,1,1]);
set(uicontrol,'style','text','position',[80 15 300 15],'string','   Gain     (click axes to display its coordinates)','HorizontalAlignment','left','BackgroundColor',[1,1,1]);
com1=['val=get(gca,''CurrentPoint'');'...
		'frequency=val(1,1);'...
		'Gain=val(1,2);'...
		'set(cursor1,''string'',num2str(frequency));'...
		'set(cursor2,''string'',num2str(Gain));'];
set(fig,'WindowButtonMotionFcn',com1);
