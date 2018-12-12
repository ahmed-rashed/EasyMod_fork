function stabdiag(ftemp,xitemp,testxi,FMAX,MaxMod,H_cols,f)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  This function displays the stabilization chart.
%
%  Input data:
%  ftemp: updated eigenvalue matrix,
%  xitemp: updated damping matrix,
%  testxi: updated test matrix,
%  FMAX: maximum freqneucy covered by the data,
%  MaxMod: maximum number of modes to calculate in the iterations,
%  H_cols: FRF matrix,
%  f: frequency vector.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


N=1;
y=N;

figure;
subplot(2,1,1);
plot(ftemp(:,N),y,'o b');
title('Stabilization diagram');
xlabel('Frequency  [Hz]');
ylabel('Number of modes');
% axis([0 max(f) 0 MaxMod]);
grid on;
hold on;
I=size(ftemp,1);

% Displaying one FRF on the chart
HdB=abs(H_cols(:,1));
facteur=max(HdB);
HdB1=HdB/facteur*(MaxMod-1);
% hold on;
plot(f,HdB1,'m');
n=0;
f=0;
d=0;
for N=2:MaxMod
   % Each frequency i of line N
   for i=1:I
      if ftemp(i,N) == 0
      else
         if ftemp(i,N-1) == 0
            % First time frequency
            y=N;
            n=n+1;
            hold on;
            plot_n=plot(ftemp(i,N),y,'b o');
        else
            % At least second time frequency
            % --> stabilized frequency
            if testxi(i,N) == 0
               % First time damping
               y=N ;
               f=f+1;
               hold on;
               plot_f=plot(ftemp(i,N),y,'g +');
           else
               % At least second time damping
               % --> stabilized frequency and damping
               y = N;
               d = d+1;
               hold on;
               plot_d=plot(ftemp(i,N),y,'r *');
           end
         end
      end
   end
end

if exist('plot_d','var') 
    legend_plot=[plot_n(1);plot_f(1);plot_d(1)];
    legend(legend_plot,' new mode',' frequency stabilization',' frequency - damping stabilization','Location','SouthEast');
else
    legend_plot=[plot_n(1);plot_f(1)];
    legend(legend_plot,' new mode',' frequency stabilization','Location','SouthEast');
end
