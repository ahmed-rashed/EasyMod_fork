function [result,Y]=releve(FTEMP,XITEMP,TESTXI,col)
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function takes in the modal parameters
%
%  Input data:
%  FTEMP: the complete matrix containing the frequency values,
%  XITEMP: the complete matrix containing the damping values,
%  TESTXI: the matrix giving information about the complete damping stabilization, 
%  col: the number of modes for which parameters are evaluated.
%
%  Output data:
%  result: results data in tabular forms
%      result(:,1)=natural frequency
%      result(:,2)=damping ratio
%      result(:,3)=stabilitisation state variable (1:
%      stabilization, 2: non stabilization),
%  Y: index matrix.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


clear y long f syn B Y g result
y=find(FTEMP(:,col));
long=length(y);
syn=zeros(long,3);
for i=1:long
   f=y(i);
   if FTEMP(f,col-1) == 0
   else
      syn(i,1:3)=[FTEMP(f,col) XITEMP(f,col) TESTXI(f,col)];
   end
end
[B,Y]=sort(syn(:,1));
for j=1:length(Y)
   g=Y(j);
   result(j,1:3)=[B(j) syn(g,2)*100 syn(g,3)];
end
disp(' ')
disp('-----------------------------------------------------------------------------');
n=sprintf('The follwing matrix (result) displays the final results for N=%1.0f',col);
disp('First column corresponds to the natural frequency (Hz)');
disp('Second column corresponds to the damping ratio (%)');
disp('last column corresponds to the stabilization state:');
disp('1 : if damping stabilization');
disp('0 : if frequency only stabilization');
disp(' ');