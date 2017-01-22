function [MatIRs]=gen_resp_impul(H,ff,str_array)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Calculation of impulse responses from FRF.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


om=2*pi*ff;
[M,N]=size(H);
for ind=1:N     
      % Calculation of the corresponding impulse response
      i=str_array(ind).response;
      j=str_array(ind).excitation;
      [x,temps]=InvFft(H(:,ind),om);
      % Memory allocation in the matrix MatIRs
      v=length(x);
      X(i,(j-1)*v+1:j*v)=x;
      % This matrix serves to visualize which impulse responses are available
      visu(i,j)=1;
end   

% The matrix is completed by zeros
m=size(visu,1);
p=size(visu,2);
diff=m-p;
if diff ~= 0
   for k=1:diff
      X(:,(p+k-1)*v+1:(p+k)*v)=zeros(m,v);
   end
end

% Time vector is added and it is identitical for all the responses
X(m+1,1:v)=temps;
MatIRs=X;