clearvars
close all
clc

% Model parameters:
M=[1 0 0; 0 1 0; 0 0 1];
C=[40 0 0; 0 40 0; 0 0 40];
K=[237315 -161000 0; -161000 398315 -161000; 0 -161000 398315];

% FRF storage
f_max=200;
N_f_max=400;
D_f=f_max/N_f_max;

f_col=(0:N_f_max-1).'*D_f;

ii_row=[1,1,1];
jj_row=[1,2,3];
n_FRF=length(ii_row);

[EigVectors_Normalized,EigValues_vec]=MDOF_Eig_Visc(M, C, K);
Receptance_cols=MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized, 2*pi*f_col, ii_row, jj_row);
[w_r_col,zeta_r_col]=MDOF_Modal_Param_Visc(EigValues_vec);
f_r_col=w_r_col/2/pi
zeta_r_col

for ii=1:n_FRF
    infoFRF(ii).response=ii;
    infoFRF(ii).dir_response=3;
    infoFRF(ii).excitation=1;
    infoFRF(ii).dir_excitation=3;
    
    unv58write(Receptance_cols(:,ii),1,3,ii,3,0,D_f,['3DL_H',num2str(jj_row(ii)),num2str(ii_row(ii)),'.unv']);
end

% % FRF Visualization
% f0=figure;
% f1=figure;
% for ii=1:n_FRF
%     for fig=[f0,f1]
%         H_label=['\alpha_{',int2str(ii_row(ii)),int2str(jj_row(ii)),'}'];
%         figure(fig)
%         ax_mag_h=subplot(2,n_FRF,ii);
%         plot_FRF_mag_phase(f_col,Receptance_cols(:,ii),0,ax_mag_h,gobjects,'',H_label);
%         hold on
%         
%         subplot(2,n_FRF,n_FRF+ii)
%         plot_FRF_Nyq(Receptance_cols(:,ii),[],H_label);
%         hold on
%     end
% end

% %% SDOF Methods
% ShowInternalDetails_circleFit=false;
% ShowInternalDetails_Dobson=false;
% f_mode_min=[40 80 110];
% f_mode_max=[60 100 130];
% n_modes=length(f_mode_min);
% f2=figure;
% infoFRF1=[];
% infoFRF2=[];
% for ii=1:n_FRF
%     label_str=['\alpha_{',int2str(infoFRF(ii).response),',',int2str(infoFRF(ii).dir_excitation),'}'];
%     
%     Receptance_Calculated1=zeros(N_f_max,1);
%     Receptance_Calculated2=zeros(N_f_max,1);
%     for jj=1:n_modes
%         LocalZone_ind=find((f_col>=f_mode_min(jj)) & (f_col<=f_mode_max(jj)));
%         Receptance_local_col=Receptance_cols(LocalZone_ind,ii);
%         f_local_col=f_col(LocalZone_ind);
%         
% %         %Circle Fit
% %         [f_r,eta_r,A_r,circ_prop]=FRF_CircleFit(f_local_col,Receptance_local_col,ShowInternalDetails_circleFit,label_str);
% %         Receptance_Calculated1=Receptance_Calculated1+A_r./complex((2*pi*f_r)^2-(2*pi*f_col).^2,eta_r*(2*pi*f_r)^2);
% %         infoFRF1=add_data(ii,f_r,eta_r,A_r,infoFRF1,jj);
%         
%         %Dobson method
%         [f_r,eta_r,A_r]=DobsonMethod(f_local_col,Receptance_local_col,ShowInternalDetails_Dobson,label_str);
%         Receptance_Calculated2=Receptance_Calculated2+A_r./complex((2*pi*f_r)^2-(2*pi*f_col).^2,eta_r*(2*pi*f_r)^2);
%         infoFRF2=add_data(ii,f_r,eta_r,A_r,infoFRF2,jj);
%     end
%     
%     figure(f1)
%     subplot(2,n_FRF,ii)
%     semilogy(f_col,abs(Receptance_Calculated1),f_col,abs(Receptance_Calculated2))
%     legend('Measured','Circle fit','Dobson method');
%    
%     subplot(2,n_FRF,n_FRF+ii)
%     plot(real(Receptance_Calculated1),imag(Receptance_Calculated1),real(Receptance_Calculated2),imag(Receptance_Calculated2))
%     legend('Measured','Circle fit','Dobson method');
% end
% % Circle Fit Results saving
% % infoMODE1=save_result_modal(infoFRF1);
% % unv55write(infoMODE1,'3DL_circle_fit.unv',1);
% 
% % Dobson method Results saving
% % infoMODE2=save_result_modal(infoFRF2);
% % unv55write(infoMODE2,'3DL_line_fit.unv',1);

%% Least-square complex exponential
[RES,infoMODE3]=lsce(Receptance_cols,f_col,infoFRF);

% % Results saving
% unv55write(infoMODE3,'3DL_LSCE.unv',1)

% Compare with Matlab implementation
N_t=2*N_f_max-2;
[D_t,f_s,~]=samplingParameters_T_N(1/D_f,N_t);
figure
modalsd(Receptance_cols,f_col,f_s,'MaxModes',10)