function [ftemp,xitemp,testxi,FNMOD,XIMOD] = rec(fn,xi,N, FMAX,ftemp,xitemp,testxi,FNMOD,XIMOD,prec1,prec2)

// ------------------   This file is part of EasyMod   ----------------------------
//  Internal function
//
//  This function applies the comparison between the frequency and damping
//  values obtained in the current step and the ones of the previous step
//
//  Input data:
//  fn: natural frequency vector,
//  xi: damping vector,
//  N: number of mode in the current step,
//  FMAX: maximum freqneucy covered by the data (for directly eliminating
//  the frequency out of range),
//  ftemp: eigenvalue matrix obtained in the previous step,
//  xitemp: damping matrix obtained in the previous step,
//  testxi: test matrix giving information about the damping stabilisation,
//  FNMOD: matrix where frequency values are saved before comparison,
//  XIMOD: matrix where damping values are saved before comparison,
//  prec1: tolerance in frequency,
//  prec2: tolerance in damping.
//
//  Output data:
//  ftemp: updated eigenvalue matrix,
//  xitemp: updated damping matrix,
//  testxi: updated test matrix,
//  FNMOD et XIMOD: updated matrices.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


//  Necessary functions:
//  -----------------------------------------------------------
//  dedoubl.m


// elimination of double (conjugate poles)
[fnmod,ximod] = dedoubl(fn,xi) ;
A = length(fnmod) ;
FNMOD(1:A,N) = fnmod(:) ;
XIMOD(1:A,N) = ximod(:) ;

if N == 1
   // first step
   A = length(fnmod) ;
   ftemp(1:A,1) = fnmod(:) ;
   xitemp(1:A,1) = ximod(:) ;
   // Next steps are identical
else
   // starting data
   A = length(fnmod) ;
   f = find(ftemp(:,N-1)) ;
   B = length(f) ;
   derf = max(f) ;   
   // Previous step is defined as start
   for a = 1:A
      // Vector test is created, whith values equal to:
      // * 1 at line b, if fnmod(a) is equal to ftemp(b)
      // * 0 otherwise.
      clear test
      test = zeros(1,B) ;
      if fnmod(a) < FMAX
         for i = 1:B
            b = f(i) ;
            if (abs(fnmod(a)-ftemp(b,N-1))/ftemp(b,N-1)) < prec1
               test(i) = 1 ;
            end
         end
         // Vector test is created
         // If vector is null, fnmod(a) does not correspond to any frequency already obtained.
         // --> fnmod(a) is added after ftemp
         // --> ximod(a) is added after xitemp
         // --> testxi(a) = 0 is imposed
         if norm(test) == 0
            derf = derf+1 ;
            ftemp(derf,N) = fnmod(a) ;
            xitemp(derf,N) = ximod(a) ;
            testxi(derf,N) = 0 ;
         end
         // If the vector contains only one 1, fnmod(a) correspond to one of the frequencies already obtained.
         // --> a means is performed
         // --> value is stored in ftemp  at the line of the corresponding frequency
         if norm(test) == 1
            clear x j
            x = find(test) ;
            j = f(x) ;
            ftemp(j,N) = (fnmod(a)+ftemp(j,N-1))/2 ;
            // The same test is performed for the damping
            if (abs(ximod(a)-xitemp(j,N-1))/xitemp(j,N-1))<prec2
               xitemp(j,N) = (ximod(a)+xitemp(j,N-1))/2 ;
               testxi(j,N) = 1 ;
            else
               xitemp(j,N) = ximod(a) ;
            end
         end
         // If the vector contains various 1, fnmod(a) correspond to the frequencies already obtained.
         // --> we look for whicg fnmod(a) is closer
         // --> a means is performed
         // --> value is stored in ftemp  at the line of the corresponding frequency
         if norm(test) > 1
            clear y k m compo pgd Y n
            y = find(test) ;
            for k = 1:length(y)
               m = f(y(k)) ;
               compo(k) = abs(fnmod(a)-ftemp(m,N-1)) ;
            end
            [pgd,Y] = min(compo) ;
            n = f(y(Y)) ;
            ftemp(n,N) = (fnmod(a)+ftemp(n,N-1))/2 ;
            // The same test is performed for the damping
            if (abs(ximod(a)-xitemp(n,N-1))/xitemp(n,N-1)) < prec2
               xitemp(n,N) = (ximod(a)+xitemp(n,N-1))/2 ;
               testxi(n,N) = 1 ;
            else
               xitemp(n,N) = ximod(a) ;
               testxi(n,N) = 0 ;
            end
         end
      else
         
      end
   end
   ftemp ;
   xitemp ;
end

endfunction