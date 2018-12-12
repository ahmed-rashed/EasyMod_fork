%% these are the examples of EasyMod 2.1.0 User's guide available at http://hosting.umons.ac.be/html/mecara/EasyMod/easymod.pdf
clearvars
global cursor1 cursor2;

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

% Check here
f_col=((0:N-1)+1).'*D_f;

ii_row=[1,1,1];
jj_row=[1,2,3];
n_FRF=length(ii_row);

receptance_cols=MDOF_FRF_slow(@(w)MDOF_FRF_Point_Visc(M, C, K, w), 2*pi*f_col, size(M,1), ii_row, jj_row);
for ii=1:n_FRF
    unv58write(-(2*pi*f_col).^2.*receptance_cols(:,ii),1,3,ii,3,0,D_f,['3DL_H',num2str(jj_row(ii)),num2str(ii_row(ii)),'.unv']);
end

% FRF loading
[H11,f_col,infoFRF(1)]=unv58read('3DL_H11.unv');
[H21,f_col,infoFRF(2)]=unv58read('3DL_H21.unv');
[H31,f_col,infoFRF(3)]=unv58read('3DL_H31.unv');
Accelerance_cols=[H11,H21,H31];

% FRF Visualization
figure;
for ii=1:n_FRF
    H_label=['\alpha_{',int2str(ii_row(ii)),int2str(jj_row(ii)),'}'];
    ax_mag_h=subplot(4,n_FRF,ii);
    ax_r=subplot(4,n_FRF,n_FRF+ii);
    ax_i=subplot(4,n_FRF,2*n_FRF+ii);
    plot_FRF_mag_phase(f_col,Accelerance_cols(:,ii),0,ax_mag_h,ax_r,'',H_label);
    plot_FRF_r_i(f_col,Accelerance_cols(:,ii),ax_r,ax_i,'',H_label);
    
    subplot(4,n_FRF,3*n_FRF+ii);
    plot_FRF_Nyq(Accelerance_cols(:,ii),[],H_label);
end

plot_ind_mode(Accelerance_cols,f_col);

%% Circle-fit
f_mode_min=[40 75 110];
f_mode_max=[60 100 130];
n_modes=length(f_mode_min);
infoFRF1=infoFRF;
for ii=1:n_FRF
    for jj=1:n_modes
        [freq_local,H_local,H_gen_local,infoMODE,circ_prop]=circle_fit(Accelerance_cols(:,ii),f_col,f_mode_min(jj),f_mode_max(jj));
        plot_circle_fit(freq_local,H_local,H_gen_local,infoMODE,circ_prop);
        infoFRF1(ii).infoMODE(jj)= infoMODE;
    end
end
% Results saving
infoMODE1=save_result_modal(infoFRF1);
unv55write(infoMODE1,'3DL_circle_fit.unv',1);

%% Line-fit
ShowInternalDetails=true;
n_FRF=size(Accelerance_cols,2);
f_mode_min=[40 75 110];
f_mode_max=[60 100 130];
n_modes=length(f_mode_min);
infoFRF2=infoFRF;
for ii=1:n_FRF
    for jj=1:n_modes
        LocalZone_ind=find((f_col>=f_mode_min(jj)) & (f_col<=f_mode_max(jj)));
        
        H_local=Accelerance_cols(LocalZone_ind,ii);
        [infoMODE,line_prop]=DobsonMethod(f_col(LocalZone_ind),Accelerance_cols(LocalZone_ind,ii),ShowInternalDetails);
        H_gen_local=A_r./complex((2*pi*f_r)^2-(2*pi*f_col).^2,eta_r*(2*pi*f_r)^2);
        %plot_line_fit(H_local,H_gen_local,infoMODE,line_prop);
        infoFRF2(ii).infoMODE(jj)=infoMODE;
    end
end
% Results saving
infoMODE2=save_result_modal(infoFRF2);
unv55write(infoMODE2,'3DL_line_fit.unv',1);

%% Least-square complex exponential
clc
% check here; I think receptance should be used instead of accelerance
Accelerance_oneSided_cols=Accelerance_cols;
Accelerance_oneSided_cols(2:end,:)=2*Accelerance_cols(2:end,:);
[RES,infoMODE3]=lsce(Accelerance_oneSided_cols,f_col,infoFRF);

% Results saving
unv55write(infoMODE3,'3DL_LSCE.unv',1)

%% Geometry definition and saving
Nodes = [1 0 0 1 ; 2 0 0 2 ; 3 0 0 3];
Connexions = [1 2 3];
unv15and82write(Nodes,Connexions,'3DL_geometry.unv',8,8);
pause

% Mode shapes visualization
figure
Methode1 = [real(infoMODE1.Bijk(3,:));real(infoMODE1.Bijk(6,:));real(infoMODE1.Bijk(9,:))];
Methode2 = [real(infoMODE2.Bijk(3,:));real(infoMODE2.Bijk(6,:));real(infoMODE2.Bijk(9,:))];
Methode3 = [real(infoMODE3.Bijk(3,:));real(infoMODE3.Bijk(6,:));real(infoMODE3.Bijk(9,:))];
subplot(3,1,1)
h1 = stem(Nodes(:,2),Methode1(:,1),'r');
hold on
h2 = stem(Nodes(:,2),Methode2(:,1),'b');
hold on
h3 = stem(Nodes(:,2),Methode3(:,1),'g');
legend([h1(1);h2(1);h3(1)],'Circle-fit','Line-fit','LSCE')
xlabel('Node positon')
ylabel('Mode shape - first mode')
axis([1 3 -4 4])

subplot(3,1,2)
h1 = stem(Nodes(:,2),Methode1(:,2),'r');
hold on
h2 = stem(Nodes(:,2),Methode2(:,2),'b');
hold on
h3 = stem(Nodes(:,2),Methode3(:,2),'g');
legend([h1(1);h2(1);h3(1)],'Circle-fit','Line-fit','LSCE')
xlabel('Node positon')
ylabel('Mode shape - second mode')
axis([1 3 -4 4])

subplot(3,1,3)
h1 = stem(Nodes(:,2),Methode1(:,3),'r');
hold on
h2 = stem(Nodes(:,2),Methode2(:,3),'b');
hold on
h3 = stem(Nodes(:,2),Methode3(:,3),'g');
legend([h1(1);h2(1);h3(1)],'Circle-fit','Line-fit','LSCE')
xlabel('Node positon')
ylabel('Mode shape - third mode')
axis([1 3 -4 4])

%% Modal assurance criterion between the identification methods
[DELTA_FREQ,MAC,coupled] = mac(infoMODE1,infoMODE2,0.8,0.1,1) ;
[DELTA_FREQ,MAC,coupled] = mac(infoMODE1,infoMODE3,0.8,0.1,1) ;
[DELTA_FREQ,MAC,coupled] = mac(infoMODE2,infoMODE3,0.8,0.1,1) ;

%% Mode collinearity
modegauss(infoMODE1) ;
modegauss(infoMODE2) ;
modegauss(infoMODE3) ;
