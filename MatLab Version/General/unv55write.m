function unv55write(infoMODE,filename,ind_complex)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function writes the information about the modal parameters in a 55
%  type UFF file. 
%
%  Synthax:
%  unv55write(infoMODE,filename,ind_complex) 
%
%  Input data:
%  filename: file name where the information is written,
%  ind_complex: boolean variable sprcifying if the saved mode shapes are
%  real (=1) or complex (=0),
%  infoMODE: structure containing the different identified parameters
%                infoMODE.frequencyk=natural frequency
%                infoMODE.etak=loss factor
%                infoMODE.Bijk=modal constant.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


%  Necessary functions:
%  -----------------------------------------------------------
%  cplxtoreal.m


if isempty(ind_complex), ind_complex=1; end;
[M,N]=size(infoMODE.frequencyk);
Nbr_mode=N;
Nbr_point=size(infoMODE.Bijk,1)/3;
name=[filename];
fid=fopen(name,'wt');
time=clock;
for ind=1:Nbr_mode;
    fprintf(fid,'    %i\n',-1);
    fprintf(fid,'    %i\n',55);
    fprintf(fid,'%s\n','Processing');
    fprintf(fid,'Analysis data universal file, EasyMod project "%s"\n',name);
    fprintf(fid,'%i-%i-%i %i:%i:%i \n',time(3),time(2),time(1),time(4),time(5),round(time(6)));
    fprintf(fid,'MODE NO. %i, FREQUENCY %8.6f(Hz), DAMPING %8.5f\n',ind,infoMODE.frequencyk(ind),infoMODE.etak(ind)/2);
    % RECORD 5
    fprintf(fid,'         %i         %i        % i        %i\n',4,3,8,13);
    % RECORD 6
    if ind_complex == 1
        fprintf(fid,'         %i         %i          %i        %i         %i         %i\n',1,2,2,8,2,3);
    else
        fprintf(fid,'         %i         %i         %i         %i         %i         %i\n',1,3,2,8,5,3);
    end
    % RECORD 7
    if ind_complex == 1
        fprintf(fid,'         %i         %i         %i         %i\n',2,4,1,ind);
    else
        fprintf(fid,'         %i         %i         %i         %i\n',2,6,1,ind);
    end
    % RECORD 8
    if ind_complex == 1
        fprintf(fid,'  %7.4e  %7.4e  %7.4e  %7.4e\n',infoMODE.frequencyk(ind),1,infoMODE.etak(ind)/2,infoMODE.etak(ind));
    else
        freqk=infoMODE.frequencyk(ind);
        xik=infoMODE.etak(ind)/2;
        pk=complex(-xik*2*pi*freqk,2*pi*freqk*sqrt(1-xik^2));
        ModeA=complex(0,1);
        ModeB=pk*ModeA;
        if (real(ModeB)>=0) && (imag(ModeB)>=0)
            fprintf(fid,' %7.4e  %7.4e  %7.4e  %7.4e  %7.4e  %7.4e\n',real(pk),imag(pk),real(ModeA),imag(ModeA),real(ModeB),imag(ModeB));
        elseif (real(ModeB)>=0) && (imag(ModeB)<0)
            fprintf(fid,' %7.4e  %7.4e  %7.4e  %7.4e  %7.4e %7.4e\n',real(pk),imag(pk),real(ModeA),imag(ModeA),real(ModeB),imag(ModeB));
        elseif (real(ModeB)<0) && (imag(ModeB)>=0)
            fprintf(fid,' %7.4e  %7.4e  %7.4e  %7.4e %7.4e  %7.4e\n',real(pk),imag(pk),real(ModeA),imag(ModeA),real(ModeB),imag(ModeB));
        elseif (real(ModeB)<0) && (imag(ModeB)<0)
            fprintf(fid,' %7.4e  %7.4e  %7.4e  %7.4e %7.4e %7.4e\n',real(pk),imag(pk),real(ModeA),imag(ModeA),real(ModeB),imag(ModeB));
        end
    end
    % RECORD 9 et RECORD 10
    if ind_complex == 1  
        for jj=1:Nbr_point
            if jj < 10
                fprintf(fid,'         %i\n',jj);    
            elseif jj < 100
                fprintf(fid,'        %i\n',jj);
            else
                fprintf(fid,'       %i\n',jj);
            end
            MODE_REAL=cplxtoreal(infoMODE);
            fprintf(fid,' %+7.4e %+7.4e %+7.4e\n',MODE_REAL.Bijk(3*(jj-1)+1,ind),MODE_REAL.Bijk(3*(jj-1)+2,ind),MODE_REAL.Bijk(3*(jj-1)+3,ind));
        end  
    else
        for jj=1:Nbr_point
            if jj < 10
                fprintf(fid,'         %i\n',jj);
            elseif jj < 100
                fprintf(fid,'        %i\n',jj);
            else
                fprintf(fid,'       %i\n',jj);
            end
            fprintf(fid,' %+7.4e %+7.4e %+7.4e %+7.4e %+7.4e %+7.4e\n',real(infoMODE.Bijk(3*(jj-1)+1,ind)),imag(infoMODE.Bijk(3*(jj-1)+1,ind)),real(infoMODE.Bijk(3*(jj-1)+2,ind)),imag(infoMODE.Bijk(3*(jj-1)+2,ind)),real(infoMODE.Bijk(3*(jj-1)+3,ind)),imag(infoMODE.Bijk(3*(jj-1)+3,ind)));
        end
    end
    fprintf(fid,'    %i\n',-1);
end
fclose(fid); 