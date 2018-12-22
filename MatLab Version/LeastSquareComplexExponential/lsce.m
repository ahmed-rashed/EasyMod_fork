function [result,infoMODE]=lsce(H_oneSided_cols,f_col,infoFRF)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Identification based on the least-square complex exponential (LSCE)
%  method (MDOF method) given the natural frequency, the loss factor 
%  and the modal constant in all the frequency range.
%
%  Synthax:
%  [result,infoFRF,infoMODE]=lsce(H,freq,infoFRF) 
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
%  result: results data in tabular forms
%      result(:,1)=natural frequency
%      result(:,2)=damping ratio
%      result(:,3)=stabilitisation state variable (1: stabilization, 2: non stabilization),
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

if ~iscolumn(f_col),error('f_col must be a column vector.'),end
N_f_max=length(f_col);
if size(H_oneSided_cols,1)~=N_f_max,error('Number of rows of H_oneSided_cols must be the same as f_col.'),end

N_inputs=max([infoFRF.excitation]);

% Sampling parameters
f_max=f_col(end);
N_t=2*N_f_max-2;  %This function assumes N_t=even. This is a valid assumption since N_t is usually power of 2.
D_f=f_col(2)-f_col(1);
[D_t,~,~]=samplingParameters_T_N(1/D_f,N_t);

h_cols=ifft_one_sided(H_oneSided_cols)*N_t; % impulse response functions

% Maximum iteration
%%%%%%%%%%%%%%Modifications by Ahmed Rashed
% disp(' ');
% disp('-----------------------------------------------------------------------------');
% disp(' ');
% N_modes_expected=input('Model order (number of expected modes) (defaults: 30):\n');    %Number of modes
% if isempty(N_modes_expected), N_modes_expected=30; end
N_modes_expected=10;

% The tolerance in frequency
%%%%%%%%%%%%%%Modifications by Ahmed Rashed
% disp(' ');
% prec1=input('Tolerance (%) in frequency (defaults: 1%):\n');
% if isempty(prec1), prec1=1; end
prec1=1;

prec1_percent=prec1/100;


% The tolerance in damping
%%%%%%%%%%%%%%Modifications by Ahmed Rashed
% disp(' ');
% prec2=input('Tolerance (%) in damping (defaults: 1%):\n');
% if isempty(prec2), prec2=1; end
prec2=1;

prec2_percent=prec2/100;

% Definition of matrices allowed from the modal parameters at each step
Z=zeros(2*N_modes_expected,N_modes_expected-1);

% Definition of matrices allowed from the comparison results at each step
f_r_temp=zeros(2*N_modes_expected,N_modes_expected-1);
FNMOD=zeros(2*N_modes_expected,N_modes_expected-1);
ZETAMOD=zeros(2*N_modes_expected,N_modes_expected-1);
zeta_r_temp=zeros(2*N_modes_expected,N_modes_expected-1);
test_zeta=zeros(2*N_modes_expected,N_modes_expected-1);

InvConditionNumber=zeros(1,N_modes_expected);
% err=zeros(1,N_modes_expected);
leastSquareError=zeros(1,N_modes_expected);
for n_mode=1:N_modes_expected
    p=modord(n_mode,N_inputs); % p represents the order of the linear differential equations
    % Overdetermined matrix definition
    h_mat_rows=MatSur2(h_cols,N_inputs,p);

    % Overdetermined equation h_T_transpose*beta=h_T_dash_transpose formulation
    h_T_dash_transpose=-h_mat_rows(:,1:N_inputs);
    h_T_transpose=h_mat_rows(:,N_inputs+1:end);

    Beta_T_transpose=pinv(h_T_transpose)*h_T_dash_transpose;     %solving using pseudoinverse [Maia, equation (4.36)]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Error calculation 
    InvConditionNumber(n_mode)=1/cond(h_T_transpose);    %Inverse of conditionnning number

    %    SingVals=svd(h_T_transpose); %Singular normalized  values
    %    err(n_mode)=1/(max(SingVals)/min(SingVals));
    %		- of least squares

    leastSquareError(n_mode)=norm(h_T_dash_transpose-(h_T_transpose*Beta_T_transpose));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Eigenvalue problem solving
    [eigVecs,V_r_col]=PbValPp(Beta_T_transpose,N_inputs,p);

    % Modal parameters extraction
    lambda_r_vec=log(V_r_col)/D_t;
    %[w_r_col,zeta_r_col]=MDOF_Modal_Param_Visc(lambda_r_vec);
    w_d_r_col=imag(lambda_r_vec);
    delta=real(lambda_r_vec);
    w_r_col=sqrt(w_d_r_col.^2+delta.^2);
    zeta_r_col=-(delta./w_r_col);
    

    % Stabilization validation (frequency and damping)
    [f_r_temp,zeta_r_temp,test_zeta,FNMOD,ZETAMOD]=rec(w_r_col/(2*pi),zeta_r_col,n_mode,f_max,f_r_temp,zeta_r_temp,test_zeta,FNMOD,ZETAMOD,prec1_percent,prec2_percent);
    % Data saving
    Z(1:N_inputs*p,n_mode)=V_r_col;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stabilization chart
figure;
subplot(2,1,1);
stabdiag(f_r_temp,test_zeta,N_modes_expected,H_oneSided_cols,f_col);

% Least squares error chart visualization
subplot(2,1,2);
yyaxis left
semilogy(leastSquareError,'.-');
ylabel('Least squares error');
grid on;

% Error chart visualization
yyaxis right
semilogy(InvConditionNumber,'.-');
ylabel('Conditioning error');
xlabel('Number of expected modes');
grid on;

% Selection of the model size
%%%%%%%%%%%%%%Modifications by Ahmed Rashed
% disp(' ');
% disp('-----------------------------------------------------------------------------');
% disp(' ');
% N_modes=input('Based on the figures, select the number of modes (defult is the number of iterations):\n');
% if isempty(N_modes), N_modes=N_modes_expected; end
N_modes=3;

%warning off

% Results statement
ind=find((f_r_temp(:,N_modes)~=0) & (f_r_temp(:,N_modes-1)~=0));
result=[f_r_temp(ind,N_modes) zeta_r_temp(ind,N_modes) test_zeta(ind,N_modes)];
disp(' ')
disp('-----------------------------------------------------------------------------');
disp(['The follwing matrix (result) displays the final results for N=',int2str(N_modes)]);
disp('First column corresponds to the natural frequency (Hz)');
disp('Second column corresponds to the damping ratio');
disp('last column corresponds to the stabilization state:');
disp('1 : if damping stabilization');
disp('0 : if frequency only stabilization');
disp(' ');
disp('     f_n        zeta   total stabilization');
disp(result);

% Mode shape extraction
[A_r_phys,w_n_r_phys,zeta_r_phys]=mode_lsce(h_cols,D_t,Z(:,N_modes),N_modes);
infoMODE=mode_stab(infoFRF,A_r_phys,w_n_r_phys/2/pi,zeta_r_phys,result,prec1_percent,prec2_percent);