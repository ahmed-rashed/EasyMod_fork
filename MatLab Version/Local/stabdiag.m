function stabdiag(f_r_mat,f_r_stabilized_mat,f_zeta_r_stabilized_mat,Receptance_cols,f_col,leastSquareError,InvConditionNumber)
% Created by Ahmed Rashed

if any(size(f_r_mat)~=size(f_r_stabilized_mat)),error('f_r_mat and f_r_stabilized_mat must have the same size'),end
if any(size(f_r_mat)~=size(f_zeta_r_stabilized_mat)),error('f_r_mat and f_r_stabilized_mat must have the same size'),end

N_modes=size(f_r_mat,2);
n_mode_row=(1:N_modes);
n_mode_mat=ones(size(f_r_mat,1),1)*n_mode_row;

figure;
axesColorRoder=colororder;

subplot(1,5,2:5);
yyaxis left
ind_f_zeta_r_stabilized_col=find(f_zeta_r_stabilized_mat);
plot(f_r_mat(ind_f_zeta_r_stabilized_col),n_mode_mat(ind_f_zeta_r_stabilized_col),'+','Color',axesColorRoder(1,:),'DisplayName','Stable in frequency \& damping');
hold on

ind_f_r_stabilized_col=find(f_r_stabilized_mat);
ind_f_r_stabilized_only_col=setdiff(ind_f_r_stabilized_col,ind_f_zeta_r_stabilized_col);
plot(f_r_mat(ind_f_r_stabilized_only_col),n_mode_mat(ind_f_r_stabilized_only_col),'o','Color',axesColorRoder(1,:),'DisplayName','Stable in frequency');

ind_f_r_not_stabilized_col=setdiff(find(~isnan(f_r_mat)),ind_f_r_stabilized_col);
plot(f_r_mat(ind_f_r_not_stabilized_col),n_mode_mat(ind_f_r_not_stabilized_col),'.','Color',axesColorRoder(1,:),'DisplayName','Not stable in frequency');
hold off

N_modes_lim=[0,1.05*N_modes];
set(gca,'YTickLabel',[],'YLim',N_modes_lim);

yyaxis right
semilogy(f_col,sum(abs(Receptance_cols),2),'DisplayName','$\sum\left|\alpha\left(f\right)\right|$');
ylabel('$\sum\left|\alpha\left(f\right)\right|$','interpreter','latex');

title('Stabilization diagram');
xlabel('$f$ (Hz)','interpreter','latex');
grid on;
legend('show','interpreter','latex')


ax1=subplot(1,5,1,'XColor',axesColorRoder(3,:),'YColor',axesColorRoder(1,:),'NextPlot','add','YGrid','on','XScale','log');
plot(leastSquareError,n_mode_row,'.-','color',axesColorRoder(3,:));    % Least squares error visualization
xlabel('Least square error');
ylabel('Number of expected modes');
ylim(N_modes_lim)

ax2=axes('Position',ax1.Position,'Color','none','XAxisLocation','top','YAxisLocation','right','XColor',axesColorRoder(4,:),'YColor',axesColorRoder(1,:),'NextPlot','add','YTickLabel',[],'XScale','log');
plot(ax2,InvConditionNumber,n_mode_row,'.-','Color',axesColorRoder(4,:));  % Condition number visualization
xlabel(ax2,'Conditioning error');
ylim(N_modes_lim)