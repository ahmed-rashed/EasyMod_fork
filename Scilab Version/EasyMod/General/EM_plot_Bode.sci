function EM_plot_Bode(frq,H,fmin,fmax)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function plots the Bode's diagram of a FRF in a specific frequency range.
//
//  Synthax:
//  EM_plot_Bode(freq,H,fmin,fmax)
//
//  Input data:
//  frq: frequency vector,
//  H: FRF in vector or matrix form (ordered in column),
//  fmin: low frequency,
//  fmin: high frequency.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


N = length(frq) ;
if (argn(2)<4) 
    fmax = N ;
end
if (argn(2)<3) 
    fmin = frq(1) ;
end
imin = 1 ;
imax = 1 ;
for ifr = 1:N
  if (frq(ifr)<fmin) 
      imin = ifr ;
  end
  if (frq(ifr)<fmax) 
      imax = ifr ;
  end
end
subplot(211)
plot(frq(imin:imax),20*log(abs(H(imin:imax,:))))
title('Bode plot')
ylabel('Magnitude  [dB]')
set(gca(),"grid",[1 1])
subplot(212)
plot(frq(imin:imax),atan(imag(H(imin:imax,:)),real(H(imin:imax,:))))
ylabel('Phase  [rad]')
xlabel('Frequency  [Hz]')
set(gca(),"grid",[1 1])

endfunction    