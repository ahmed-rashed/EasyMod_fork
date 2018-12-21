function [f_n_temp,zeta_temp,testzeta,FNMOD,ZETAMOD]=rec(f_n_r,zeta_r,N, f_max,f_n_temp,zeta_temp,testzeta,FNMOD,ZETAMOD,prec1,prec2)
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function applies the comparison between the frequency and damping
%  values obtained in the current step and the ones of the previous step
%
%  Input data:
%  f_n_r: natural frequency vector,
%  zeta_r: damping vector,
%  N: number of mode in the current step,
%  f_max: maximum freqneucy covered by the data (for directly eliminating the frequency out of range),
%  f_n_temp: eigenvalue matrix obtained in the previous step,
%  zeta_temp: damping matrix obtained in the previous step,
%  testzeta: test matrix giving information about the damping stabilisation,
%  FNMOD: matrix where frequency values are saved before comparison,
%  XIMOD: matrix where damping values are saved before comparison,
%  prec1: tolerance in frequency,
%  prec2: tolerance in damping.
%
%  Output data:
%  f_n_temp: updated eigenvalue matrix,
%  zeta_temp: updated damping matrix,
%  testzeta: updated test matrix,
%  FNMOD et XIMOD: updated matrices.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


%  Necessary functions:
%  -----------------------------------------------------------
%  dedoubl.m


% elimination of double (conjugate poles)
[f_n_mod,i_f_n_r]=uniquetol(f_n_r,1e-7/max(abs(f_n_r)));
zeta_mod=zeta_r(i_f_n_r);
A=length(f_n_mod);
FNMOD(1:A,N)=f_n_mod;
ZETAMOD(1:A,N)=zeta_mod;

if N == 1
   % first step
   f_n_temp(1:A,1)=f_n_mod;
   zeta_temp(1:A,1)=zeta_mod;
   % Next steps are identical
else
   % starting data
   f=find(f_n_temp(:,N-1));
   B=length(f);
   derf=max(f);   
   % Previous step is defined as start
   for a=1:A
      % Vector test is created, whith values equal to:
      % * 1 at line b, if fnmod(a) is equal to f_n_temp(b)
      % * 0 otherwise.
      test=false(1,B);
      if f_n_mod(a) < f_max
         for i=1:B
            b=f(i);
            if (abs(f_n_mod(a)-f_n_temp(b,N-1))/f_n_temp(b,N-1)) < prec1
               test(i)=true;
            end
         end
         % Vector test is created
         % If vector is null, fnmod(a) does not correspond to any frequency already obtained.
         % --> fnmod(a) is added after f_n_temp
         % --> ximod(a) is added after zeta_temp
         % --> testzeta(a)=0 is imposed
         if all(~test)
            derf=derf+1;
            f_n_temp(derf,N)=f_n_mod(a);
            zeta_temp(derf,N)=zeta_mod(a);
            testzeta(derf,N)=0;
         end
         nn=length(find(test));
         if nn==1
            % If the vector contains only one 1, fnmod(a) correspond to one of the frequencies already obtained.
            % --> a means is performed
            % --> value is stored in f_n_temp  at the line of the corresponding frequency

            j=f(test);
            f_n_temp(j,N)=(f_n_mod(a)+f_n_temp(j,N-1))/2;
            % The same test is performed for the damping
            if (abs(zeta_mod(a)-zeta_temp(j,N-1))/zeta_temp(j,N-1))<prec2
               zeta_temp(j,N)=(zeta_mod(a)+zeta_temp(j,N-1))/2;
               testzeta(j,N)=1;
            else
               zeta_temp(j,N)=zeta_mod(a);
            end
         elseif nn > 1
            % If the vector contains various 1, fnmod(a) correspond to the frequencies already obtained.
            % --> we look for whicg fnmod(a) is closer
            % --> a means is performed
            % --> value is stored in f_n_temp  at the line of the corresponding frequency
            y=find(test);
            comp=zeros(1,length(y));
            for k=1:length(y)
               m=f(y(k));
               comp(k)=abs(f_n_mod(a)-f_n_temp(m,N-1));
            end
            [~,Y]=min(comp);
            n=f(y(Y));
            f_n_temp(n,N)=(f_n_mod(a)+f_n_temp(n,N-1))/2;
            % The same test is performed for the damping
            if (abs(zeta_mod(a)-zeta_temp(n,N-1))/zeta_temp(n,N-1)) < prec2
               zeta_temp(n,N)=(zeta_mod(a)+zeta_temp(n,N-1))/2;
               testzeta(n,N)=1;
            else
               zeta_temp(n,N)=zeta_mod(a);
               testzeta(n,N)=0;
            end
         end
      else
         
      end
   end
end

