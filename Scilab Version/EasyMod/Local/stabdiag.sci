function stabdiag(ftemp,xitemp,testxi,FMAX,MaxMod,H,f)

// ------------------   This file is part of EasyMod   ----------------------------
//  Internal function
//
//  This function displays the stabilization chart.
//
//  Input data:
//  ftemp: updated eigenvalue matrix,
//  xitemp: updated damping matrix,
//  testxi: updated test matrix,
//  FMAX: maximum freqneucy covered by the data,
//  MaxMod: maximum number of modes to calculate in the iterations,
//  H: FRF matrix,
//  f: frequency vector.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


N = 1 ;
y = N ;
fmax = max(f) ;
figure ;
subplot(2,2,1) ;
plot(-1,-1,'o b',-1,-1,'g+',-1,-1,'r*') ;
a = get("current_axes") ;
a.data_bounds = [0,0;fmax,MaxMod] ;
set(gca(),"grid",[1 1]) ;
set(gca(),"auto_clear","off") ;

plot(ftemp(:,N),y,'o b') ;
title('Stabilization diagram') ;
xlabel('Frequency  [Hz]') ;
ylabel('Number of modes') ;
I = size(ftemp,1) ;

// Displaying one FRF on the chart
HdB = abs(H(:,1)) ;
facteur = (max(HdB)) ;
HdB1 = HdB/facteur*(MaxMod - 1) ;
set(gca(),"auto_clear","off") ;
plot(f,HdB1,'m') ;
n = 0 ;
f = 0 ;
d = 0 ;
for N = 2:MaxMod
   // Each frequency i of line N
   for i = 1:I
      if ftemp(i,N) == 0
      else
         if ftemp(i,N-1) == 0
            // First time frequency
            y = N ;
            n = n+1 ;
            set(gca(),"auto_clear","off") ;
            plot_n = plot(ftemp(i,N),y,'bo') ;
        else
            // At least second time frequency
            // --> stabilized frequency
            if testxi(i,N) == 0
               // First time damping
               y = N  ;
               f = f+1 ;
               set(gca(),"auto_clear","off") ;
               plot_f = plot(ftemp(i,N),y,'g+') ;
           else
               // At least second time damping
               // --> stabilized frequency and damping
               y  =  N ;
               d  =  d+1 ;
               set(gca(),"auto_clear","off") ;
               plot_d = plot(ftemp(i,N),y,'r*') ;
           end
         end
      end
   end
end
legend(' new mode',' frequency stabilization',' frequency - damping stabilization',-4) ;

endfunction