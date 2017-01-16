function EM_plot_Nyquist(frq,H,fmin,fmax)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function plots the Nyquist's diagram of a FRF in a specific frequency range.
//
//  Synthax:
//  M_plot_Nyquits(frq,H,fmin,fmax)
//
//  Input data:
//  freq: frequency vector,
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
plot(real(H(imin:imax,:)),imag(H(imin:imax,:)),'-o')
title('Nyquist plot')
xlabel('Real part')
ylabel('Imaginary part')
set(gca(),"grid",[1 1])

endfunction