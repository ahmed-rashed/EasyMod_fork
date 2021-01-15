function [lsce_result,infoMODE]=lsce(Receptance_cols,f_col,infoFRF)
%Modifications by Ahmed Rashed
%
% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  Identification based on the least-square complex exponential (LSCE)
%  method (MDOF method) given the natural frequency, the loss factor 
%  and the modal constant in all the frequency range.
%
%  Synthax:
%  [lsce_result,infoFRF,infoMODE]=lsce(H,freq,infoFRF) 
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
%  lsce_result: results data in tabular forms
%      lsce_result(:,1)=natural frequency
%      lsce_result(:,2)=damping ratio
%      lsce_result(:,3)=stabilitisation state variable (1: stabilization, 2: non stabilization),
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
if size(Receptance_cols,1)~=N_f_max,error('Number of rows of Receptance_cols must be the same as f_col.'),end

N_inputs=max([infoFRF.excitation]);

% Sampling parameters
f_max=f_col(end);
N_t=2*N_f_max-2;  %This function assumes N_t=even. This is a valid assumption since N_t is usually power of 2.
D_f=f_col(2)-f_col(1);
[D_t,f_s,~]=samplingParameters_T_N(1/D_f,N_t);

h_cols=ifft(Receptance_cols,N_t,1,'symmetric')*N_t; % impulse response functions

% Maximum iteration
N_modes_expected=input('Model order (number of expected modes): (default: 10)');
if isempty(N_modes_expected)
    N_modes_expected=10;
end

% The tolerance in frequency
prec_f_r=input('Percentage frequency Tolerance (%): (default: 1%)')/100;
if isempty(prec_f_r)
    prec_f_r=1/100;
end

% The tolerance in damping
prec_zeta_r=input('Percentage damping Tolerance (%): (default: 5%)')/100;
if isempty(prec_zeta_r)
    prec_zeta_r=5/100;
end

% Matrices for comparing results at each step
f_r_mat=nan(2*N_modes_expected,N_modes_expected);
zeta_r_mat=nan(2*N_modes_expected,N_modes_expected);
f_r_stabilized_mat=false(2*N_modes_expected,N_modes_expected);
f_zeta_r_stabilized_mat=f_r_stabilized_mat;

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
    [~,V_r_col]=PbValPp(Beta_T_transpose,N_inputs,p);

    % Modal parameters extraction
    lambda_r_col=log(V_r_col)/D_t;
%     [w_r_col,zeta_r_col]=MDOF_Modal_Param_Visc(lambda_r_col);
    w_d_r_col=imag(lambda_r_col); 
    delta_r_col=real(lambda_r_col); 
    w_r_col=sqrt(w_d_r_col.^2+delta_r_col.^2); 
    zeta_r_col=-(delta_r_col./w_r_col); 
    ind1_col=find(w_r_col<2*pi*f_max);
    [f_r_unique_col,ind2_col]=uniquetol(w_r_col(ind1_col)/2/pi,1e-7/max(abs(w_r_col/2/pi)));    % Eliminate repeated natural frequencies (resulting from conjugate poles)
    zeta_r_unique_col=zeta_r_col(ind1_col(ind2_col));

    if n_mode>1 % Stabilization validation (frequency and damping)
        f_r_old_col=f_r_new_col;
        zeta_r_old_col=zeta_r_new_col;
        [f_r_new_col,zeta_r_new_col,f_r_stabilized_col,f_zeta_r_stabilized_col]=rec(f_r_unique_col,zeta_r_unique_col,f_r_old_col,zeta_r_old_col,prec_f_r,prec_zeta_r);
    else    %n_mode=1
        f_r_new_col=f_r_unique_col;
        zeta_r_new_col=zeta_r_unique_col;
        f_r_stabilized_col=false(size(zeta_r_unique_col));
        f_zeta_r_stabilized_col=f_r_stabilized_col;
    end
    N_r_unique=length(f_r_unique_col);
    f_r_mat(1:N_r_unique,n_mode)=f_r_new_col;
    zeta_r_mat(1:N_r_unique,n_mode)=zeta_r_new_col;
    f_r_stabilized_mat(1:N_r_unique,n_mode)=f_r_stabilized_col;
    f_zeta_r_stabilized_mat(1:N_r_unique,n_mode)=f_zeta_r_stabilized_col;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stabilization chart
stabdiag(f_r_mat,f_r_stabilized_mat,f_zeta_r_stabilized_mat,Receptance_cols,f_col,leastSquareError,InvConditionNumber);

% Selection of the model size
N_modes=input('Based on the figures, select the number of modes: (defult is the number of iterations)\n');
if isempty(N_modes)
    N_modes=N_modes_expected;
else
    if N_modes>N_modes_expected
        error('Number of modes cannot exceed the number of iterations.')
    end
end

warning off

% Results statement
lsce_result=[f_r_new_col(f_r_stabilized_col),zeta_r_new_col(f_r_stabilized_col),f_zeta_r_stabilized_col(f_r_stabilized_col)];

disp('        f_n        zeta         f_n & zeta stabilization state');
disp(lsce_result);

% Mode shape extraction
[A_r_phys,w_n_r_phys,zeta_r_phys]=mode_lsce(h_cols,D_t,V_r_col,N_modes);

infoMODE=mode_stab(infoFRF,A_r_phys,w_n_r_phys/2/pi,zeta_r_phys,lsce_result,prec_f_r,prec_zeta_r);