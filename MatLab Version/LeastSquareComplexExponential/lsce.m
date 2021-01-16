function [f_r_col,zeta_r_col,f_zeta_r_stabilized_col,A_r]=lsce(Receptance_cols,f_col,N_inputs, ...
                                                                    f_r_tol,zeta_r_tol)   %Optional arguments
%Created by Ahmed Rashed based on EasyMod

if ~iscolumn(f_col),error('f_col must be a column vector.'),end
N_f_max=length(f_col);
if size(Receptance_cols,1)~=N_f_max,error('Number of rows of Receptance_cols must be the same as f_col.'),end

if nargin<4 % Natural frequency tolerance
    f_r_tol=1/100;
else
    if isempty(f_r_tol)
        f_r_tol=1/100;
    end
end

if nargin<5 % zeta tolerance
    zeta_r_tol=5/100;
else
    if isempty(zeta_r_tol)
        zeta_r_tol=5/100;
    end
end

% Sampling parameters
f_max=f_col(end);
N_t=2*N_f_max-2;  %This function assumes N_t=even. This is a valid assumption since N_t is usually power of 2.
D_f=f_col(2)-f_col(1);
[D_t,~,~]=samplingParameters_T_N(1/D_f,N_t);

h_cols=ifft(Receptance_cols,N_t,1,'symmetric')*N_t; % impulse response functions

% Maximum iteration
N_modes_expected=input('Model order (number of expected modes): (default: 10)');
if isempty(N_modes_expected)
    N_modes_expected=10;
end

% Matrices for comparing results at each step
f_r_cell=cell(N_modes_expected,1);
zeta_r_cell=cell(N_modes_expected,1);
stabilized_f_r_cell=cell(N_modes_expected,1);
stabilized_f_zeta_r_cell=stabilized_f_r_cell;
InvConditionNumber_col=zeros(N_modes_expected,1);
% err=zeros(1,N_modes_expected);
leastSquareError_col=zeros(N_modes_expected,1);
V_r_cell=cell(1,N_modes_expected,1);
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
    InvConditionNumber_col(n_mode)=1/cond(h_T_transpose);    %Inverse of conditionnning number

    %    SingVals=svd(h_T_transpose); %Singular normalized  values
    %    err(n_mode)=1/(max(SingVals)/min(SingVals));
    %		- of least squares

    leastSquareError_col(n_mode)=norm(h_T_dash_transpose-(h_T_transpose*Beta_T_transpose));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Eigenvalue problem solving
    [~,V_r_cell{n_mode}]=PbValPp(Beta_T_transpose,N_inputs,p);

    % Modal parameters extraction
    lambda_r_col=log(V_r_cell{n_mode})/D_t;
%     [w_r_temp1_col,zeta_r_temp1_col]=MDOF_Modal_Param_Visc(lambda_r_col);
    w_d_r_temp1_col=imag(lambda_r_col); 
    delta_r_temp1_col=real(lambda_r_col); 
    w_r_temp1_col=sqrt(w_d_r_temp1_col.^2+delta_r_temp1_col.^2); 
    zeta_r_temp1_col=-(delta_r_temp1_col./w_r_temp1_col); 
    ind1_col=find(w_r_temp1_col<2*pi*f_max);
    [f_r_temp2_col,ind2_col]=uniquetol(w_r_temp1_col(ind1_col)/2/pi,1e-7/max(abs(w_r_temp1_col/2/pi)));    % Eliminate repeated natural frequencies (resulting from conjugate poles) and sorts them
    zeta_r_temp2_col=zeta_r_temp1_col(ind1_col(ind2_col));

    if n_mode>1 % Stabilization validation (frequency and damping)
        f_r_old_col=f_r_new_col;
        zeta_r_old_col=zeta_r_new_col;
        [f_r_new_col,zeta_r_new_col,f_r_stabilized_new_col,f_zeta_r_stabilized_new_col]=rec(f_r_temp2_col,zeta_r_temp2_col,f_r_old_col,zeta_r_old_col,f_r_tol,zeta_r_tol);
    else    %n_mode=1
        f_r_new_col=f_r_temp2_col;
        zeta_r_new_col=zeta_r_temp2_col;
        f_r_stabilized_new_col=false(size(zeta_r_temp2_col));
        f_zeta_r_stabilized_new_col=f_r_stabilized_new_col;
    end
    f_r_cell{n_mode}=f_r_new_col;
    zeta_r_cell{n_mode}=zeta_r_new_col;
    stabilized_f_r_cell{n_mode}=f_r_stabilized_new_col;
    stabilized_f_zeta_r_cell{n_mode}=f_zeta_r_stabilized_new_col;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stabilization chart
stabilizationDiagram(f_r_cell,stabilized_f_r_cell,stabilized_f_zeta_r_cell,Receptance_cols,f_col,leastSquareError_col,InvConditionNumber_col);

% Selection of the model size
N_modes=input('Based on the figures, select the number of modes: (default is the number of iterations)');
if isempty(N_modes)
    N_modes=N_modes_expected;
else
    if N_modes>N_modes_expected
        error('Number of modes cannot exceed the number of iterations.')
    elseif 2*N_modes>N_t
        error('Number of samples is insufficient for the requested number of modes!')
    end
end

warning off

% Results
f_r_col=f_r_cell{N_modes}(stabilized_f_r_cell{N_modes});
zeta_r_col=zeta_r_cell{N_modes}(stabilized_f_r_cell{N_modes});
f_zeta_r_stabilized_col=stabilized_f_zeta_r_cell{N_modes}(stabilized_f_r_cell{N_modes});

% Mode shape extraction
A_r=mode_lsce(h_cols,V_r_cell{n_mode},N_modes,ind2_col);

warning on