function [wmod,ximod]=dedoubl(w,xi)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function supresses the redundant values of frequency and damping
%  vectors while conserving the links between these two vectors.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


d=1;
n=length(w);
wmod(1)=w(1);
ximod(1)=xi(1);
for i=2:n
   clear testw
   testw(1)=0;
   for j=1:i-1
      if abs(w(i)-w(i-j)) < 1e-7
         testw(1+j)=1;
      else
         testw(1+j)=0;
      end
   end
   if norm(testw) == 0
      d=d+1;
      wmod(d)=w(i);
      ximod(d)=xi(i); 
   end
end
wmod=wmod';
ximod=ximod';