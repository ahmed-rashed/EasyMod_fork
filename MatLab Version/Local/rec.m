function [f_r_temp,zeta_r_temp,test_zeta,FNMOD,ZETAMOD]=rec(f_r_col,zeta_r_col,n_mode,f_max,f_r_temp,zeta_r_temp,test_zeta,FNMOD,ZETAMOD,prec1,prec2)
% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function applies the comparison between the frequency and damping
%  values obtained in the current step and the ones of the previous step
%
%  Input data:
%  f_r_col: natural frequency vector,
%  zeta_r_col: damping vector,
%  n_mode: number of mode in the current step,
%  f_max: maximum freqneucy covered by the data (for directly eliminating the frequency out of range),
%  f_r_temp: eigenvalue matrix obtained in the previous step,
%  zeta_r_temp: damping matrix obtained in the previous step,
%  test_zeta: test_vec matrix giving information about the damping stabilisation,
%  FNMOD: matrix where frequency values are saved before comparison,
%  XIMOD: matrix where damping values are saved before comparison,
%  prec1: tolerance in frequency,
%  prec2: tolerance in damping.
%
%  Output data:
%  f_r_temp: updated eigenvalue matrix,
%  zeta_r_temp: updated damping matrix,
%  test_zeta: updated test_vec matrix,
%  FNMOD and XIMOD: updated matrices.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


%  Necessary functions:
%  -----------------------------------------------------------
%  dedoubl.m


% elimination of double (conjugate poles)
[f_r_unique,i_f_r_unique]=uniquetol(f_r_col,1e-7/max(abs(f_r_col)));
zeta_r_unique=zeta_r_col(i_f_r_unique);
N=length(f_r_unique);
FNMOD(1:N,n_mode)=f_r_unique;
ZETAMOD(1:N,n_mode)=zeta_r_unique;

if n_mode == 1
   % first step
   f_r_temp(1:N,1)=f_r_unique;
   zeta_r_temp(1:N,1)=zeta_r_unique;
else
   % starting data
   ind=find(f_r_temp(:,n_mode-1)~=0);
   N_f_r_last=length(ind);
   ind_max=max(ind);   
   % Previous step is defined as start
   for n=1:N
      % Vector test_vec is created, with values equal to:
      % * 1 at line ind(i), if fnmod(n) is equal to f_r_temp(ind(i))
      % * 0 otherwise.
      test_vec=false(1,N_f_r_last);
      if f_r_unique(n) < f_max
         for i=1:N_f_r_last
            if (abs(f_r_unique(n)-f_r_temp(ind(i),n_mode-1))/f_r_temp(ind(i),n_mode-1)) < prec1
               test_vec(i)=true;
            end
         end
         % Vector test_vec is created
         % If vector is null, fnmod(n) does not correspond to any frequency already obtained.
         % --> fnmod(n) is added after f_r_temp
         % --> ximod(n) is added after zeta_r_temp
         % --> test_zeta(n)=0 is imposed
         if all(~test_vec)
            ind_max=ind_max+1;
            f_r_temp   (ind_max,n_mode)=f_r_unique(n);
            zeta_r_temp(ind_max,n_mode)=zeta_r_unique(n);
            test_zeta  (ind_max,n_mode)=0;
         end
         nn=length(find(test_vec));
         if nn==1
            % If the vector contains only one 1, fnmod(n) correspond to one of the frequencies already obtained.
            % --> n means is performed
            % --> value is stored in f_r_temp  at the line of the corresponding frequency

            j=ind(test_vec);
            f_r_temp(j,n_mode)=(f_r_unique(n)+f_r_temp(j,n_mode-1))/2;
            % The same test_vec is performed for the damping
            if (abs(zeta_r_unique(n)-zeta_r_temp(j,n_mode-1))/zeta_r_temp(j,n_mode-1))<prec2
               zeta_r_temp(j,n_mode)=(zeta_r_unique(n)+zeta_r_temp(j,n_mode-1))/2;
               test_zeta(j,n_mode)=1;
            else
               zeta_r_temp(j,n_mode)=zeta_r_unique(n);
            end
         elseif nn > 1
            % If the vector contains various 1, fnmod(n) correspond to the frequencies already obtained.
            % --> we look for which fnmod(n) is closer
            % --> n means is performed
            % --> value is stored in f_r_temp  at the line of the corresponding frequency
            y=find(test_vec);
            comp=zeros(1,length(y));
            for k=1:length(y)
               m=ind(y(k));
               comp(k)=abs(f_r_unique(n)-f_r_temp(m,n_mode-1));
            end
            [~,Y]=min(comp);
            nn=ind(y(Y));
            f_r_temp(nn,n_mode)=(f_r_unique(n)+f_r_temp(nn,n_mode-1))/2;
            % The same test_vec is performed for the damping
            if (abs(zeta_r_unique(n)-zeta_r_temp(nn,n_mode-1))/zeta_r_temp(nn,n_mode-1)) < prec2
               zeta_r_temp(nn,n_mode)=(zeta_r_unique(n)+zeta_r_temp(nn,n_mode-1))/2;
               test_zeta(nn,n_mode)=1;
            else
               zeta_r_temp(nn,n_mode)=zeta_r_unique(n);
               test_zeta(nn,n_mode)=0;
            end
         end
      else
         
      end
   end
end

