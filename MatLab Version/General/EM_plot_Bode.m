function EM_plot_Bode(freq,H,fmin,fmax)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function plots the Bode's diagram of a FRF in a specific frequency range.
%
%  Synthax:
%  EM_plot_Bode(freq,H,fmin,fmax)
%
%  Input data:
%  freq: frequency vector,
%  H: FRF in vector or matrix form (ordered in column),
%  fmin: low frequency,
%  fmin: high frequency.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


N=length(freq);
if (nargin<4) 
    fmax=N;
end
if (nargin<3) 
    fmin=freq(1);
end
imin=1;
imax=1;
for ifr=1:N
  if (freq(ifr)<fmin) 
      imin=ifr;
  end
  if (freq(ifr)<fmax) 
      imax=ifr;
  end
end
subplot(211)
plot(freq(imin:imax),20*log(abs(H(imin:imax,:))))
title('Bode plot')
ylabel('Magnitude  [dB]')
grid on
subplot(212)
plot(freq(imin:imax),angle(H(imin:imax,:)))
ylabel('Phase  [rad]')
xlabel('Frequency  [Hz]')
grid on
        