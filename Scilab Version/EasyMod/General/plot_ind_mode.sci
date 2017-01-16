function [cursor1,cursor2] = plot_ind_mode(I,frq)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function displays the three mode indicators - sum of FRFs, sum of FRFs
//  real part and sum of FRFs imaginary part.
//
//  Synthax:
//  [cursor1,cursor2] = plot_ind_mode(I,frq)
//
//  Input data:
//  frq: frequency vector,
//  I: structure containing the different indicators
//        I.ISUM: sum of included FRFs,
//        I.ISRe: sum of included real part FRFs,
//        I.ISIm: sum of included imaginary part FRFs.
//
//  Output data:
//  cursor1 and cursor2: unused in SciLab version.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


disp('(1): sum of included FRFs ') ;
disp('(2): sum of included real part FRFs ') ;
disp('(3): sum of included imaginary part FRFs ') ;
di = input('Select?','s') ;

select di
    case '1'
        figure ;
        for i =1:length(frq)
            if abs(I.ISUM(i)) == 0 then
                I.ISUM(i) = 1D-9 ;
            end
        end
        plot(frq,20*log10(I.ISUM)) ;
        set(gca(),"grid",[1 1]) ;
        xlabel('Frequency  [Hz]') ;
        ylabel('sum of included FRFs  [dB]') ;
    case '2'
        figure ;
        for i =1:length(frq)
            if abs(I.ISRe(i)) == 0 then
                I.ISRe(i) = 1D-9 ;
            end
        end
        plot(frq,20*log10(abs(I.ISRe)))
        set(gca(),"grid",[1 1]) ;
        xlabel('Frequency  [Hz]')
        ylabel('sum of included real part FRFs  [dB]')
    case '3'
        figure ;
        for i =1:length(frq)
            if abs(I.ISIm(i)) == 0 then
                I.ISIm(i) = 1D-9 ;
            end
        end
        plot(frq,20*log10(abs(I.ISIm)))
        set(gca(),"grid",[1 1]) ;
        xlabel('Frequency  [Hz]')
        ylabel('sum of included imaginary part FRFs  [dB]')
    else
        disp('Invalid request')
end
cursor1 = [] ;
cursor2 = [] ;


endfunction     