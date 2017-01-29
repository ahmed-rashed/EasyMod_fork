function [lsce_res,infoMODE]=lsce(H_oneSided_cols,f_col,infoFRF)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Identification based on the least-square complex exponential (LSCE)
%  method (MDOF method) giving the natural frequency, the loss factor 
%  and the modal constant in all the frequency range.
%
%  Synthax:
%  [lsce_res,infoFRF,infoMODE]=lsce(H,freq,infoFRF) 
%
%  Input data:
%  H: FRF matrix containing all the FRFs (column number=FRF number),
%  freq: frequency vector,
%  infoFRF: structure containing information on FRFs 
%            infoFRF(jj).response=jjth FRF response node
%            infoFRF(jj).dir_response=jjth FRF response direction (1=X, 2=Y,
%             3=Z, 4=RotX, 5=RotY, 6=RotZ)
%            infoFRF(jj).excitation=jjth FRF excitation node
%            infoFRF(jj).dir_excitation=jjth FRF excitation direction
%            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ). 
%
%  Output data:
%  lsce_res: results data in tabular forms
%      lsce_res(:,1)=natural frequency
%      lsce_res(:,2)=damping ratio
%      lsce_res(:,3)=stabilitisation state variable (1: stabilization, 2: non stabilization),
%  infoMODE: structure containing the different identified parameters
%                infoMODE.frequencyk=natural frequency
%                infoMODE.etak=loss factor
%                infoMODE.Bijk=modal constant,
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


%  Necessary functions:
%  -----------------------------------------------------------
%  gen_resp_impul.m
%  AnMatIR.m
%  Modord.m
%  MatSur2.m
%  rec.m
%  releve.m
%  stabdiag.m
%  mode_lsce.m
%  mode_stab.m


% Calculation of the impulse response matrix from frequency data
D_f=f_col(2)-f_col(1);
N_f_max=size(H_oneSided_cols,1)-1;
f_max=N_f_max*D_f;
[h_cols,N_inp,N_out,N_t,D_t]=calculate_h(H_oneSided_cols,D_f,infoFRF);

% Calculation parameters capture
%		- maximum iteration
disp(' ');
disp('-----------------------------------------------------------------------------');
disp(' ');
MaxMod=input('Model order - maximum iteration to analyse (defaults: 30):\n');    %Number of modes
if isempty(MaxMod), MaxMod=30; end;
%		- the tolerance in frequency
disp(' ');
prec1=input('Tolerance (%) in frequency (defaults: 1%):\n');
if isempty(prec1), prec1=1; end;
prec1=prec1/100;
%		- the tolerance in damping
disp(' ');
prec2=input('Tolerance (%) in damping (defaults: 1%):\n');
if isempty(prec2), prec2=1; end; 
prec2=prec2/100;

% Definition of matrices allowed from the modal parameters at each step
WD=zeros(2*MaxMod,MaxMod-1);
WN=zeros(2*MaxMod,MaxMod-1);
ZETA=zeros(2*MaxMod,MaxMod-1);
LL=zeros(N_inp*(MaxMod-1),MaxMod*2);
Z=zeros(2*MaxMod,MaxMod-1);

% Definition of matrices allowed from the comparison results at each step
FTEMP=zeros(2*MaxMod,MaxMod-1);
FNMOD=zeros(2*MaxMod,MaxMod-1);
ZETAMOD=zeros(2*MaxMod,MaxMod-1);
XITEMP=zeros(2*MaxMod,MaxMod-1);
TESTXI=zeros(2*MaxMod,MaxMod-1);

InvCond=zeros(1,MaxMod);
% err=zeros(1,MaxMod);
epsilon=zeros(1,MaxMod);
% Check a while loop here
for N=1:MaxMod
   % Data
   p=modord(N,N_inp); % p represents the order of the linear differential equations 
   % Overdetermined matrix definition
   G_cols=MatSur2(h_cols,N_inp,p);

   % Overdetermined equation A*beta=B formulation
   B=-G_cols(:,1:N_inp);
   A=G_cols(:,N_inp+1:end);
   % Resolving (using the pseudoinverse function)
   beta_cols=pinv(A)*B;     %[Maia, equation (4.12)]
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Error calculation 
   %		- Inverse of conditionnning number
   InvCond(N)=1/cond(A);
%    %		- Singular normalized  values
%    SingVals=svd(A); 
%    err(N)=1/(max(SingVals)/min(SingVals));
   %		- of least squares
   epsilon(N)=norm(B-(A*beta_cols));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % Eigenvalue problem resolving
   [L,z_col]=PbValPp(beta_cols,N_inp,p);
   % Modal parameters extraction
   lambda_r=log(z_col)./D_t;
   w_d_r=imag(lambda_r);
   delta_r=real(lambda_r);
   w_n_r=sqrt(w_d_r.^2+delta_r.^2);
   zeta_r=-delta_r./w_n_r;
   % Stabilization validation (frequency and damping)
   [FTEMP,XITEMP,TESTXI,FNMOD,ZETAMOD]=rec(w_n_r/(2*pi),zeta_r,N,f_max,FTEMP,XITEMP,TESTXI,FNMOD,ZETAMOD,prec1,prec2);
%    % Data saving
%    WD(1:length(w_d_r),N)=w_d_r;
%    WN(1:length(w_n_r),N)=w_n_r;
%    ZETA(1:length(zeta_r),N)=zeta_r;
%    LL((N-1)*N_inp+1:N*N_inp,1:size(L,2))=L;
   Z(1:length(z_col),N)=z_col;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stabilization chart
stabdiag(FTEMP,XITEMP,TESTXI,f_max,MaxMod,H_oneSided_cols,f_col);

% Least squares error chart visualization
subplot(2,1,2);
yyaxis left
semilogy(epsilon,'-.');
ylabel('Least squares error');
grid on;
zoom on;

% Error chart visualization
subplot(2,2,4);
semilogy(InvCond,'-.');
xlabel('Number of modes');
ylabel('Conditioning error');
grid on;
zoom on;

% Selection of the model size
disp(' ');
disp('-----------------------------------------------------------------------------');
disp(' ');
N_modes=input('Based on the figures, select the number of modes (defult is the number of iterations):\n');
if isempty(N_modes), N_modes=MaxMod; end;
%warning off

% Results statement
lsce_res=releve(FTEMP,XITEMP,TESTXI,N_modes);
kk=find(lsce_res(:,1)>0); % Only non-null frequencies are considered
M=size(lsce_res,1);
lsce_res=lsce_res(kk:M,:);
disp('     f_n        zeta   total stabilization');
disp(lsce_res);

% Mode shape extraction
[A_r_phys,w_n_r_phys,zeta_r_phys]=mode_lsce(h_cols,D_t,Z(:,N_modes),N_modes);
infoMODE=mode_stab(infoFRF,A_r_phys,w_n_r_phys/2/pi,zeta_r_phys,lsce_res,prec1,prec2);