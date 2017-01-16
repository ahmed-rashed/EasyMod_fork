function [H,Ni,nddl,Np,dt] = AnMatIR(mat)

// ------------------   This file is part of EasyMod   ----------------------------
//  Internal function
//
//  Analysis of the impulse response matrix. It is composed of impulse
//  responses placed on the matrix at (i,j) location corresponding to their
//  DOF. The function calculates:
//  - a reduced matrix only with the useful columns,
//  - the number of inputs,
//  - the number of DOF (supposed to be the number of outputs),
//  - the number of measurement locations,
//  - the time resolution .
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


// Definition of the number of DOF
clear nddl
nddl = size(mat,1)-1 ; 
a = sprintf('The studied structure has %1.0f degrees of freedom',nddl);
disp(' ')
disp('-----------------------------------------------------------------------------');
disp(' ')
disp(a)

// nombre de points de mesure
clear dl Np I temps dt
dl = find(mat(nddl+1,:)) ;
Np = max(dl) ;
b = sprintf('%1.0f measurement samples are enumerated',Np) ;
disp(' ') ;
disp(b) ;
I = mat(1:nddl,:) ;
temps = mat(nddl+1,1:Np) ;
dt = temps(2)-temps(1) ;

// test on the impulse responses
clear ncol test i j 
ncol = size(mat,2)/Np ;
for i = 1:nddl
   for j = 1:ncol
      if norm(mat(i,(j-1)*Np+1:j*Np)) ~= 0
         test(i,j) =1 ;
      else
         test(i,j) =0 ;
      end
   end
end
disp('h matrix = ') ;
disp(test) ;
disp('(availability of the matrix h)') ;
disp('if h(i,j) = 1 then impulse response is available for measurement at point i and excitation at point j') ;

// Test for obtaining the number of input Ni
clear Ni s H c  
Ni = 0 ;
disp(' ') ;
for s = 1:ncol
   if norm(test(:,s)) == 0
   else
      Ni = Ni+1;
      H(1:nddl,(Ni-1)*Np+1:Ni*Np) = mat(1:nddl,(s-1)*Np+1:s*Np) ;
      c = sprintf('Input at DOF %1.0f',s) ;
      disp(c) ;
   end
end

endfunction