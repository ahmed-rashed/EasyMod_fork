function [cursor1,cursor2] = plot_ind_mode(I,freq)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function displays the three mode indicators - sum of FRFs, sum of FRFs
%  real part and sum of FRFs imaginary part.
%
%  Synthax:
%  [cursor1,cursor2] = plot_ind_mode(I,freq)
%
%  Input data:
%  freq: frequency vector,
%  I: structure containing the different indicators
%        I.ISUM: sum of included FRFs,
%        I.ISRe: sum of included real part FRFs,
%        I.ISIm: sum of included imaginary part FRFs.
%
%  Output data:
%  cursor1 and cursor2: parameters related to the cursor (x,y) of plot.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


disp('(1): sum of included FRFs ') ;
disp('(2): sum of included real part FRFs ') ;
disp('(3): sum of included imaginary part FRFs ') ;
di = input('Select?','s') ;

switch di
    case '1'
        gca = figure ;
        plot(freq,20*log10(I.ISUM)) ;
        grid on ;
        xlabel('Frequency  [Hz]') ;
        ylabel('sum of included FRFs  [dB]') ;
        cursor1 = uicontrol('style','text','position',[10 10 60 15]) ;
        cursor2 = uicontrol('style','text','position',[80 10 60 15]) ;
        set(uicontrol,'style','text','position',[10 30 60 15],'string','Frequency [Hz]:') ;
        set(uicontrol,'style','text','position',[80 30 60 15],'string','Gain [dB]:') ;
        com1 = [...
                'val = get(gca,''CurrentPoint'') ;'...
                'frequency = val(1,1) ;'...
                'Gain = val(1,2) ;'...
                'set(cursor1,''string'',num2str(frequency)) ;'...
                'set(cursor2,''string'',num2str(Gain)) ;'...
            ] ;
        set(gca,'WindowButtonMotionFcn',com1) ;
    case '2'
        gca = figure ;
        plot(freq,20*log10(abs(I.ISRe)))
        grid on
        xlabel('Frequency  [Hz]')
        ylabel('sum of included real part FRFs  [dB]')
        cursor1 = uicontrol('style','text','position',[10 10 60 15]) ;
        cursor2 = uicontrol('style','text','position',[80 10 60 15]) ;
        set(uicontrol,'style','text','position',[10 30 60 15],'string','Frequency [Hz]:') ;
        set(uicontrol,'style','text','position',[80 30 60 15],'string','Gain [dB]:') ;
        com1 = [...
                'val = get(gca,''CurrentPoint'') ;'...
                'frequency = val(1,1) ;'...
                'Gain = val(1,2) ;'...
                'set(cursor1,''string'',num2str(frequency)) ;'...
                'set(cursor2,''string'',num2str(Gain)) ;'...
            ] ;
        set(gca,'WindowButtonMotionFcn',com1) ;
    case '3'
        gca = figure ;
        plot(freq,20*log10(abs(I.ISIm)))
        grid on
        xlabel('Frequency  [Hz]')
        ylabel('sum of included imaginary part FRFs  [dB]')
        cursor1 = uicontrol('style','text','position',[10 10 60 15]) ;
        cursor2 = uicontrol('style','text','position',[80 10 60 15]) ;
        set(uicontrol,'style','text','position',[10 30 60 15],'string','Frequency [Hz]:') ;
        set(uicontrol,'style','text','position',[80 30 60 15],'string','Gain [dB]:') ;
        com1 = [...
                'val = get(gca,''CurrentPoint'') ;'...
                'frequency = val(1,1) ;'...
                'Gain = val(1,2) ;'...
                'set(cursor1,''string'',num2str(frequency)) ;'...
                'set(cursor2,''string'',num2str(Gain)) ;'...
            ] ;
        set(gca,'WindowButtonMotionFcn',com1) ;
    otherwise
        disp('Invalid request')
end
        