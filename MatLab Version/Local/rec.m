function [f_r_temp,zeta_r_temp,bDampingFreq_Stabilized]=rec(f_r_col,zeta_r_col,n_mode,f_max,f_r_temp,zeta_r_temp,bDampingFreq_Stabilized,prec1,prec2)
%  This function applies the comparison between the frequency and damping
%  values obtained in the current step and the ones of the previous step
%
%  Input data:
%  f_r_col: natural frequency vector,
%  zeta_r_col: damping vector,
%  n_mode: number of mode in the current step,
%  f_max: maximum freqneucy covered by the data (for directly eliminating the frequency out of range),
%  prec1: tolerance in frequency,
%  prec2: tolerance in damping.
%
%  Output data:
%  f_r_temp: updated eigenvalue matrix,
%  zeta_r_temp: updated damping matrix,
%  bDampingFreq_Stabilized: test matrix,
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT

% elimination of double (conjugate poles)
[f_r_unique,i_f_r_unique]=uniquetol(f_r_col,1e-7/max(abs(f_r_col)));
zeta_r_unique=zeta_r_col(i_f_r_unique);
N_f_r_unique=length(f_r_unique);

if n_mode==1    % first step
   f_r_temp(1:N_f_r_unique,1)=f_r_unique;
   zeta_r_temp(1:N_f_r_unique,1)=zeta_r_unique;
else    % Next steps are identical
   ind_vec=find(f_r_temp(:,n_mode-1)~=0);
   ind_max=max(ind_vec);

   % Previous step is defined as start
   for n_f_r_unique=1:N_f_r_unique
      if f_r_unique(n_f_r_unique) < f_max
         test_vec=(abs(f_r_unique(n_f_r_unique)-f_r_temp(ind_vec,n_mode-1))./f_r_temp(ind_vec,n_mode-1)) < prec1;
         ind1_vec=find(test_vec);
         nn=length(ind1_vec);
         if nn==0  % Otherwise, f_r_unique(n_f_r_unique) does not correspond to any frequency already obtained.
            ind_max=ind_max+1;
            f_r_temp   (ind_max,n_mode)=f_r_unique(n_f_r_unique);
            zeta_r_temp(ind_max,n_mode)=zeta_r_unique(n_f_r_unique);
            bDampingFreq_Stabilized(ind_max,n_mode)=false;
         elseif nn==1   % f_r_unique(n_f_r_unique) correspond to one of the frequencies already obtained.
            % --> n_f_r_unique means is performed
            % --> value is stored in f_r_temp  at the line of the corresponding frequency
            j_vec=ind_vec(ind1_vec);
            f_r_temp(j_vec,n_mode)=(f_r_unique(n_f_r_unique)+f_r_temp(j_vec,n_mode-1))/2;
            % The same test_vec is performed for the damping
            if (abs(zeta_r_unique(n_f_r_unique)-zeta_r_temp(j_vec,n_mode-1))/zeta_r_temp(j_vec,n_mode-1)) < prec2
               zeta_r_temp(j_vec,n_mode)=(zeta_r_unique(n_f_r_unique)+zeta_r_temp(j_vec,n_mode-1))/2;
               bDampingFreq_Stabilized(j_vec,n_mode)=true;
            else
               zeta_r_temp(j_vec,n_mode)=zeta_r_unique(n_f_r_unique);
            end
         elseif nn > 1  % f_r_unique(n_f_r_unique) correspond to the frequencies already obtained.
            % --> we look for which f_r_unique(n_f_r_unique) is closer
            % --> n_f_r_unique means is performed
            % --> value is stored in f_r_temp  at the line of the corresponding frequency
            comp=abs(f_r_unique(n_f_r_unique)-f_r_temp(ind_vec(ind1_vec),n_mode-1));
            [~,Y]=min(comp);
            nnn=ind_vec(ind1_vec(Y));
            f_r_temp(nnn,n_mode)=(f_r_unique(n_f_r_unique)+f_r_temp(nnn,n_mode-1))/2;
            % The same test_vec is performed for the damping
            if (abs(zeta_r_unique(n_f_r_unique)-zeta_r_temp(nnn,n_mode-1))/zeta_r_temp(nnn,n_mode-1)) < prec2
               zeta_r_temp(nnn,n_mode)=(zeta_r_unique(n_f_r_unique)+zeta_r_temp(nnn,n_mode-1))/2;
               bDampingFreq_Stabilized(nnn,n_mode)=true;
            else
               zeta_r_temp(nnn,n_mode)=zeta_r_unique(n_f_r_unique);
               bDampingFreq_Stabilized(nnn,n_mode)=false;
            end
         end
      end
   end
end

