function [h_rows,N_inp,N_out,N_t,D_t]=AnMatIR(h_mat)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Analysis of the impulse response matrix. It is composed of impulse
%  responses placed on the matrix at (i,j) location corresponding to their
%  DOF. The function calculates:
%  - a reduced matrix only with the useful columns,
%  - the number of inputs,
%  - the number of DOF (supposed to be the number of outputs),
%  - the number of measurement locations,
%  - the time resolution .
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


% Definition of the number of DOF
N_out=size(h_mat,1)-1; 

% nombre de points de mesure
N_t=max(find(h_mat(end,:)));
D_t=h_mat(end,2)-h_mat(end,1);

% test on the impulse responses
N_col=size(h_mat,2)/N_t;
for i=1:N_out
   for j=1:N_col
      if norm(h_mat(i,(j-1)*N_t+1:j*N_t)) ~= 0
         test(i,j) =1;
      else
         test(i,j) =0;
      end
   end
end
disp('h matrix=');
disp(test);
disp('(availability of the matrix h)');
disp('if h(i,j)=1 then impulse response is available for measurement at point i and excitation at point j');

% Test for obtaining the number of input Ni
N_inp=0;
disp(' ');
for s=1:N_col
   if norm(test(:,s))
      N_inp=N_inp+1;
      h_rows(:,(N_inp-1)*N_t+1:N_inp*N_t)=h_mat(1:N_out,(s-1)*N_t+1:s*N_t);
      c=sprintf('Input at DOF %1.0f',s);
      disp(c);
   end
end