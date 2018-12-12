function DobsonInternalDetails(f_local_vec,alpha_local_vec,f_r,eta_r,A_r,line_prop,label_str)
w_local_vec=2*pi*f_local_vec;
w_r=2*pi*f_r;

B_r=0;  %Correct this

alpha_gen_local_vec=A_r./complex(w_r^2-w_local_vec.^2,eta_r*w_r^2)+B_r;

% Lines properties
Delta=line_prop.Delta;
tr=line_prop.tr;
ti=line_prop.ti;
ur=line_prop.ur;
dr=line_prop.dr;
ui=line_prop.ui;
di=line_prop.di;

figure

% Bode plot
ax_mag_h=subplot(6,2,[2,4,6,8]-1);
ax_phase_h=subplot(6,2,[10,12]-1);
[ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(f_local_vec,[alpha_local_vec,alpha_gen_local_vec],false,ax_mag_h,ax_phase_h,'',label_str);
temp1=sort([get(ax_phase_h,'XTick'),f_r]);
ind_r=find(temp1==f_r);
set(ax_phase_h,'XTick',temp1);
set(ax_mag_h,'XTick',temp1);
temp2=get(ax_phase_h,'XTickLabel');
temp2{ind_r}=['{\itf_{r}}=',num2str(f_r,'%.2f' )];
set(ax_phase_h,'XTickLabel',temp2);
ax_phase_h.XAxis.FontName='Times';

legend(ax_mag_h,'Measured','Dobson method','Location','southeast');

% Various lines visualization
subplot(4,4,3+[0,4]);
plot(f_local_vec.^2,real(Delta));
ylabel('$\Re\left(\Delta\right)$', 'interpreter', 'latex')
axis tight

subplot(4,4,4+[0,4]);
plot(f_local_vec.^2,imag(Delta));
ylabel('$\Im\left(\Delta\right)$', 'interpreter', 'latex')
axis tight

subplot(4,4,11+[0,4]);
plot(f_local_vec.^2,tr);
hold on;
plot(f_local_vec.^2,ur*w_local_vec.^2+dr);
ylabel('$t_{\mathrm{R}}$', 'interpreter', 'latex');
xlabel('$f^{2}$  (Hz\textsuperscript{2})', 'interpreter', 'latex');
axis tight

subplot(4,4,12+[0,4]);
plot(f_local_vec.^2,ti);
hold on;
plot(f_local_vec.^2,ui*w_local_vec.^2+di);
ylabel('$t_{\mathrm{I}}$', 'interpreter', 'latex');
xlabel('$f^{2}$  (Hz\textsuperscript{2})', 'interpreter', 'latex');
axis tight