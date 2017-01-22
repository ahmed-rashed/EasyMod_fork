function [freq_local,H_local,H_gen_local,infoMODE,circ_prop]=circle_fit(H,freq,fmin,fmax);

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Identification based on the circle-fit method (SDOF method) giving the
%  natural frequency, the loss factor and the modal constant.
%
%  Synthax:
%  [freq_local,H_local,H_gen_local,infoMODE,circ_prop]=circle_fit(H,freq,fmin,fmax)
%
%  Input data:
%  H: FRF vector that we want to analyse,
%  freq: frequency vector,
%  fmin and fmax: frequency range embracing the mode that we want to
%                 analyse.
%
%  Output data:
%  freq_local: studied frequency vector extrated from the overall frequency,
%  H_local: studied FRF vector extrated from the FRF,
%  H_gen_local: synthetized FRF from modal parameters,
%  infoMODE: structure containing the different identified parameters
%                infoMODE.frequencyk=natural frequency
%                infoMODE.etak=loss factor
%                infoMODE.Bijk=modal constant.
%  circ_prop: structure containing information about Nyquist's circle
%          circ_prop.x0=abscissa of circle center
%          circ_prop.y0=ordinate of circle center
%          circ_prop.R0=circle radius
%          circ_prop.D =damping evolution.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


%  Necessary functions: 
%  -----------------------------------------------------------
%  err_fit_circle.m
%  modify_ref_angle.m
%  interpol.m
%  estimate_natural_frequency.m


problem=0;
temp=find(freq>fmin);
index_low=temp(1,1)-1;
temp=find(freq>fmax);
index_high=temp(1,1)-1;

% Data extraction
H_local=H(index_low:index_high);
freq_local=freq(index_low:index_high);
ww_mode=2*pi*freq_local;
x=real(H_local);
y=imag(H_local);

% In the case of insufficient number of samples
N_pts=length(freq_local); 
if N_pts < 5
    disp('!!!')
    disp('Insufficient number of samples in the studied frequency range')
    disp(' ')
    problem=1;
else
  % Best circle finding
  [x0,y0,R0]=err_fit_circle(x,y);
  theta_brut=atan2(y-y0,x-x0);
  theta_brut=mod(theta_brut,2*pi);
  [theta,theta_ref]=modify_ref_angle(theta_brut);
  for jj=1:(N_pts-1)
     delta(jj)=theta(jj+1)-theta(jj);
     [f_interpol,theta_interpol]=interpol(freq_local(jj),theta(jj),freq_local(jj+1),theta(jj+1));
     f_delta(jj)=f_interpol;
     theta_delta(jj)=theta_interpol; 
  end
  for jj=1:N_pts-2
     delta2(jj)=delta(jj+1)-delta(jj);
     [f_interpol,theta_interpol]=interpol(f_delta(jj),theta_delta(jj),f_delta(jj+1),theta_delta(jj+1));
     f_delta2(jj)=f_interpol;
     theta_delta2(jj)=theta_interpol;
  end

  % Local and overall maximum finding
  sign_ref=sign(delta2(1,1));
  kk=0; % maximum counter
  
  ind_zero=[];
  for index=1:length(delta2)
      if sign(delta2(1,index)) == -sign_ref;
          kk=kk+1;
          sign_ref=-sign_ref;
          ind_zero(kk)=index;
      end
  end
   
  % Case where there is no resonance
  if isempty(ind_zero) == 1
      disp('!!!')
      disp('No resonance in the studied frequency range')
      disp(' ')
      problem=1;
  else
     [fr, theta_r]=estimate_natural_frequency(ind_zero,delta2,f_delta2,theta_delta2);
     ss=find(delta==max(delta));
     kkk=find(ind_zero==ss);
  
     % Natural frequency caractéristics
     eig_theta=theta_r(kkk);
     eig_frequency=fr(kkk);
     wr=2*pi*eig_frequency;
 
     if isempty(eig_frequency) == 1
        disp('!!!') 
        disp('No resonance in the studied frequency range')
        disp(' ')
        problem=1;
     else
         
        % Damping (loss factor) calculation
        temp=find(freq_local>eig_frequency);
        k1=temp(1,1)-1;
        k2=temp(1,1);
        wb=freq_local(k1:-1:1)*2*pi;  % half-power low frequency 
        wa=freq_local(k2:N_pts)*2*pi;  % half-power high frequency
        theta_b=(eig_theta-theta(k1:-1:1));
        theta_a=(theta(k2:N_pts)-eig_theta);
        for jj=1:length(wb)
           for kk=1:length(wa)
              if  theta_b(jj) > pi
                  continue
              elseif theta_a(kk) > pi
                  continue
              end
              N=(wa(kk)^2-wb(jj)^2);
              D=wr^2*(tan(theta_a(kk)/2)+ tan(theta_b(jj)/2));
              xi(jj,kk)=N/D;
           end
        end
        loss_factor=mean(mean(xi));

        % Mode shape calculation
        ref=2*pi-(theta_r(kkk)+theta_ref);
        xa=R0*cos(ref)+x0;
        ya=R0*sin(ref)+y0;
        phi=-atan2(xa-x0,ya-y0);
        Bmod=2*R0*loss_factor*wr*wr;
        B=complex(Bmod*cos(phi),Bmod*sin(phi));
             
        % Residue calculation
        xb=R0*cos(ref+pi)+x0;
        yb=R0*cos(ref+pi)+x0;
        Rijk=complex(xb,yb);
        % FRF synthetization from mosal parameters  
        DEN=complex(wr^2-ww_mode.^2,loss_factor*wr^2);
        if problem == 0
            H_gen_local=(B./DEN);
        end
     end
  end
end

% Data saving
if problem == 0
    circ_prop=struct('x0',x0,'y0',y0,'R0',R0,'D',xi,'theta',theta,'theta_ref',theta_ref,'eig_theta',eig_theta);
    infoMODE=struct('frequencyk',eig_frequency,'etak',loss_factor,'Bijk',B);
elseif problem == 1
    circ_prop=struct('x0',NaN,'y0',NaN,'R0',NaN,'D',NaN,'theta',NaN,'theta_ref',NaN,'eig_theta',NaN);
    infoMODE=struct('frequencyk',NaN,'etak',NaN,'Bijk',0);
    H_gen_local=0;
end



