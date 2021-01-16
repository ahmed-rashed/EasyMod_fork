function stabilizationDiagram(f_r_cell,stabilized_f_r_cell,stabilized_f_zeta_r_cell,Receptance_cols,f_col,leastSquareError_col,InvConditionNumber_col)
% Created by Ahmed Rashed

N_modes=length(f_r_cell);
if length(stabilized_f_r_cell)~=N_modes,error('f_r_cell and stabilized_f_r_cell must have the same size'),end
if size(stabilized_f_zeta_r_cell)~=N_modes,error('f_r_cell and stabilized_f_r_cell must have the same size'),end

n_mode_col=(1:N_modes).';

f_r_combined_col=cell2mat(f_r_cell);
n_mode_combined_col=cell2mat(cellfun(@(x,y) y*ones(length(x),1),f_r_cell,num2cell(n_mode_col),'UniformOutput',false));
stabilized_f_r_combined_col=cell2mat(stabilized_f_r_cell);
stabilized_f_zeta_r_combined_col=cell2mat(stabilized_f_zeta_r_cell);

figure;
axesColorRoder=colororder;

subplot(1,5,2:5);
yyaxis left
ind_f_zeta_r_stabilized_col=find(stabilized_f_zeta_r_combined_col);
plot(f_r_combined_col(ind_f_zeta_r_stabilized_col),n_mode_combined_col(ind_f_zeta_r_stabilized_col),'+','Color',axesColorRoder(1,:),'DisplayName','Stable in frequency \& damping');
hold on

ind_f_r_stabilized_col=find(stabilized_f_r_combined_col);
ind_f_r_stabilized_only_col=setdiff(ind_f_r_stabilized_col,ind_f_zeta_r_stabilized_col);
plot(f_r_combined_col(ind_f_r_stabilized_only_col),n_mode_combined_col(ind_f_r_stabilized_only_col),'o','Color',axesColorRoder(1,:),'DisplayName','Stable in frequency');

ind_f_r_not_stabilized_col=setdiff(find(~isnan(f_r_combined_col)),ind_f_r_stabilized_col);
plot(f_r_combined_col(ind_f_r_not_stabilized_col),n_mode_combined_col(ind_f_r_not_stabilized_col),'.','Color',axesColorRoder(1,:),'DisplayName','Not stable in frequency');
hold off

N_modes_lim=[0,1.05*N_modes];
set(gca,'YTickLabel',[],'YLim',N_modes_lim);

yyaxis right
semilogy(f_col,sum(abs(Receptance_cols),2),'DisplayName','$\sum\left|\alpha\left(f\right)\right|$');
ylabel('$\sum\left|\alpha\left(f\right)\right|$','interpreter','latex');

title('Stabilization diagram');
xlabel('$f$ (Hz)','interpreter','latex');
grid on;
legend('show','interpreter','latex','location','southeast')

ax1=subplot(1,5,1,'XColor',axesColorRoder(3,:),'YColor',axesColorRoder(1,:),'NextPlot','add','YGrid','on','XScale','log');
plot(leastSquareError_col,n_mode_col,'.-','color',axesColorRoder(3,:));    % Least squares error visualization
xlabel('Least square error');
ylabel('Number of expected modes');
ylim(N_modes_lim)

ax2=axes('Position',ax1.Position,'Color','none','XAxisLocation','top','YAxisLocation','right','XColor',axesColorRoder(4,:),'YColor',axesColorRoder(1,:),'NextPlot','add','YTickLabel',[],'XScale','log');
plot(ax2,InvConditionNumber_col,n_mode_col,'.-','Color',axesColorRoder(4,:));  % Condition number visualization
xlabel(ax2,'Conditioning error');
ylim(N_modes_lim)