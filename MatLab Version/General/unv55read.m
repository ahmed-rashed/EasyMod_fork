function infoMODE=unv55read(filename,No)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function reads the 55 type UFF file containing the information
%  about the modal parameters.
%
%  Synthax:
%  infoMODE=unv55read(filename,No) 
%
%  Input data:
%  filename: file name where the information is readed,
%  No: number of experimental nodes.
%
%  Output data:
%  infoMODE: structure containing the different identified parameters
%                infoMODE.frequencyk=natural frequency
%                infoMODE.etak=loss factor
%                infoMODE.Bijk=modal constant.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


nom_unv=[filename];
frequency='FREQUENCY';
freq=[];
amort=[];
temp=[];
tempo=[];
modal=[];
fichier=fopen(nom_unv,'r');
num=0;
while feof(fichier) ~= 1
    valec=fscanf(fichier,'%s',[1]);
    if strcmp('Residuals',valec) == 1 
        break         
    end
    if strcmp(frequency,valec) == 1 
        num=num+1; % counter
        frq=fscanf(fichier,'%g',[1]);
        freq(num)=frq; % natural frequency storing
        valec=fscanf(fichier,'%s',[2]);
        temp=fscanf(fichier,'%g',[1]);
        xi(num)=temp; % damping ratio storing
        valec=fscanf(fichier,'%i',[8]);
        ind_cplx=fscanf(fichier,'%i',[1]);
        NDV=fscanf(fichier,'%i',[1]); % number of data per node
        if ind_cplx == 2 
            valec=fscanf(fichier,'%g',[8]);
            Nbr=No*(NDV+1);
            temp=fscanf(fichier,'%g',[Nbr]);
            temp=temp';
            [M,N]=size(temp);
            for index=1:NDV+1:N-NDV
                tempo=[tempo temp(index+1:index+NDV)];
            end
            psi(num,:)=tempo;
            tempo=[];
        elseif ind_cplx == 5
            valec=fscanf(fichier,'%g',[10]);
            Nbr=No*(2*NDV+1);
            temp=fscanf(fichier,'%g',[Nbr]);
            temp=temp';
            [M,N]=size(temp);
            for index=1:2*NDV+1:N-2*NDV
                x=complex(temp(index+1),temp(index+2));
                y=complex(temp(index+3),temp(index+4));
                z=complex(temp(index+5),temp(index+6));
                xyz=[x y z];
                tempo=[tempo xyz];
            end
            psi(num,:)=tempo;
            tempo=[];
        end
    end
end
fclose(fichier);
psi=psi';
infoMODE=struct('frequencyk',freq,'etak',xi*2,'Bijk',psi);