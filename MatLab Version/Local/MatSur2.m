function [GG] = MatSur2(H,No,Ni,Nt,p)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function builds the overdetermined matrix from
%    - the reduced matrix resulting from the analyse of the impulse
%    responses matrix,
%    - the number of outputs and inputs (No and Ni),
%    - the last measurement point taken into account (Npmax),
%    - the number of equations by output (Neq),
%    - the number of measurement points (Nt),
%    - the order (p).
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


Neq=floor(Nt/(p+1)); 
Npmax=Neq*(p+1); 
% Creation of g (or [Hn]i)
clear i n j g; 
for i=1:No
   for n=1:Npmax
      for j=1:Ni
            g(:,j,n,i)=H(i,((j-1)*Nt+n));
      end
   end
end

% Sub-matrices G assembly
clear i l m n G; 
for i=1:No
   for l=1:Neq 
      n=(p+1)*l;  
      for m=1:(p+1) 
         G(l,(m-1)*Ni+1:m*Ni,i)=g(:,:,n-(m-1),i); 
      end 
   end
end

% Overdetermined matrix GG assembly
clear i GG
for i=1:No
   GG((i-1)*Neq+1:i*Neq,:)=G(:,:,i); 
end