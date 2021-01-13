function stabdiag(f_r_temp,bDampingFreq_Stabilized,N_modes,H_oneSided_cols,f_col)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function displays the stabilization chart.
%
%  Input data:
%  f_r_temp: updated eigenvalue matrix,
%  bDampingFreq_Stabilized: updated test matrix,
%  N_modes: maximum number of modes to calculate in the iterations,
%  H_oneSided_cols: FRF matrix,
%  f_col: frequency vector.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT

% Displaying one FRF on the chart
yyaxis right
handle=plot(f_col,abs(H_oneSided_cols(:,1)));
set(get(get(handle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
ylabel('$|\alpha(f)|$','interpreter','latex');

yyaxis left
n_mode=1;
plot_n=plot(f_r_temp(1,n_mode),n_mode,'b o','DisplayName','New mode');
hold on;
for n_mode=2:N_modes
    for i_f_r_temp=1:size(f_r_temp,1)    % Each frequency i of line n_mode
        if f_r_temp(i_f_r_temp,n_mode)==0
            continue
        elseif f_r_temp(i_f_r_temp,n_mode-1)==0  % First time frequency --> New mode
            plot_n.XData=[plot_n.XData,f_r_temp(i_f_r_temp,n_mode)];
            plot_n.YData=[plot_n.YData,n_mode];
        else    % Successive to first time frequency --> frequency already stabilized
            if bDampingFreq_Stabilized(i_f_r_temp,n_mode)==0    % First damping
                if ~exist('plot_f','var')
                    plot_f=plot(f_r_temp(i_f_r_temp,n_mode),n_mode,'g +','DisplayName','frequency stabilized');
                else
                    plot_f.XData=[plot_f.XData,f_r_temp(i_f_r_temp,n_mode)];
                    plot_f.YData=[plot_f.YData,n_mode];
                end
            else % Successive to first damping --> stabilized frequency and damping
                if ~exist('plot_f_zeta','var')
                    plot_f_zeta=plot(f_r_temp(i_f_r_temp,n_mode),n_mode,'r *','DisplayName','frequency & damping stabilized');
                else
                    plot_f_zeta.XData=[plot_f_zeta.XData,f_r_temp(i_f_r_temp,n_mode)];
                    plot_f_zeta.YData=[plot_f_zeta.YData,n_mode];
                end
            end
        end
    end
end

title('Stabilization diagram');
xlabel('$f$ (Hz)','interpreter','latex');
ylabel('Number of expected modes');
grid on;
legend('show')