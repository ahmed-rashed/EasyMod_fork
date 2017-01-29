clearvars
close all
clc

% Model parameters:
M=[1 0 0; 0 1 0; 0 0 1];
C=[40 0 0; 0 40 0; 0 0 40];
K=[237315 -161000 0; -161000 398315 -161000; 0 -161000 398315];

% FRF storage
f_max=200;
N=400;
D_f=f_max/N;

f_col=(0:N-1).'*D_f;

ii_row=[1,1,1];
jj_row=[1,2,3];
n_FRF=length(ii_row);

[EigVectors_Normalized, EigValues_mat]=MDOF_Eig_Visc(M, C, K);
Receptance_cols=MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, 2*pi*f_col, ii_row, jj_row);
for ii=1:n_FRF
    unv58write(Receptance_cols(:,ii),1,3,ii,3,0,D_f,['3DL_H',num2str(jj_row(ii)),num2str(ii_row(ii)),'.unv']);
end

% FRF Visualization
f0=figure;
f1=figure;
for ii=1:n_FRF
    for fig=[f0,f1]
        H_label=['\alpha_{',int2str(ii_row(ii)),int2str(jj_row(ii)),'}'];
        figure(fig)
        ax_mag_h=subplot(2,n_FRF,ii);
        plot_FRF_mag_phase(f_col,Receptance_cols(:,ii),0,ax_mag_h,gobjects,'',H_label);
        hold on
        
        subplot(2,n_FRF,n_FRF+ii)
        plot_FRF_Nyq(Receptance_cols(:,ii),H_label);
        %coloured_line_3d(real(Receptance(:,ii)),imag(Receptance(:,ii)),zeros(size(Receptance(:,ii))),f_col);view(2)
        hold on
    end
end

%% Circle-fit
ShowInternalDetails=true;
f_mode_min=[40 80 110];
f_mode_max=[60 100 130];
n_modes=length(f_mode_min);
f2=figure;
for ii=1:n_FRF
    Receptance_Calculated=zeros(N,1);
    for jj=1:n_modes
        Receptance_temp=Receptance_cols(:,ii);
        LocalZone_flag=(f_col>=f_mode_min(jj)) & (f_col<=f_mode_max(jj));
        Receptance_local=Receptance_temp(LocalZone_flag);
        freq_local=f_col(LocalZone_flag);
        
        %Circle Fit
        [Mode_Inf,circ_prop]=FRF_CircleFit(freq_local,Receptance_local,ShowInternalDetails)
        
        Receptance_Calculated=Receptance_Calculated+Mode_Inf.A_r./(complex((2*pi*Mode_Inf.f_r)^2-(2*pi*f_col).^2,Mode_Inf.eta_r*(2*pi*Mode_Inf.f_r)^2));
        
        %Receptance_local visualization
        figure(f2)
        ax=subplot(n_FRF,n_modes,(ii-1)*n_modes+jj);
        visualizeLocalReceptance(freq_local,Receptance_local,circ_prop,ii_row(ii),jj_row(ii),ax);
        
        %infoFRF=add_data(i,infoMODE.frequencyk,infoMODE.etak,infoMODE.Bijk,infoFRF,j);
    end
    
    figure(f1)
    subplot(2,n_FRF,ii)
    semilogy(f_col,abs(Receptance_Calculated))
    legend('Measured','Calculated');
   
    subplot(2,n_FRF,n_FRF+ii)
    %coloured_line_3d(real(Receptance(:,ii)),imag(Receptance(:,ii)),zeros(size(Receptance(:,ii))),f_col)
    plot(real(Receptance_Calculated),imag(Receptance_Calculated))
    legend('Measured','Calculated');
end
% % Results saving
% infoMODE1=save_result_modal(infoFRF);
% unv55write(infoMODE1,'3DL_circle_fit.unv',1);

%% Line-fit
N=size(Receptance_cols,2);
bornes_min[40 75 110];
bornes_max[60 100 130];
Nbr_modelength(bornes_min);
for i=1:N
    for j=1:Nbr_mode
        [freq_local,H_local,H_gen_local,infoMODE,circ_prop]=line_fit(Receptance_cols(:,i),f_col,bornes_min(j),bornes_max(j));
        plot_line_fit(freq_local,H_local,H_gen_local,infoMODE,circ_prop);
        infoFRF2=add_data(i,infoMODE.frequencyk,infoMODE.etak,infoMODE.Bijk,infoFRF,j);
    end
end
% Results saving
infoMODE2=save_result_modal(infoFRF2);
unv55write(infoMODE2,'3DL_line_fit.unv',1);

%% Least-square complex exponential
[RES,infoFRF3,infoMODE3]=lsce(Receptance_cols,f_col,infoFRF);

% Results saving
unv55write(infoMODE3,'3DL_LSCE.unv',1)