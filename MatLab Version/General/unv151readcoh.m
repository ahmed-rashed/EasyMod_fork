function [coh,fref,infoFRF]=unv151readcoh(filename)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function reads the 151 type UFF file in order to extract only the
%  coherences.
%
%  Synthax:
%  [coh,freq,infoFRF]=unv151readcoh(filename) 
%
%  Input data:
%  filename: file name where the information is readed,
%
%  Output data:
%  H: coherence matrix,
%  freq: frequency vector,
%  infoFRF: structure containing information on the associated FRFs 
%            infoFRF(jj).response=jjth FRF response node
%            infoFRF(jj).dir_response=jjth FRF response direction (1=X, 2=Y,
%             3=Z, 4=RotX, 5=RotY, 6=RotZ)
%            infoFRF(jj).excitation=jjth FRF excitation node
%            infoFRF(jj).dir_excitation=jjth FRF excitation direction
%            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ). 
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


name_file=[filename];
fid=fopen(name_file,'r');
ss=['opening file ',name_file];
disp(ss);
key='Coherence';
NONE='NONE';
nbr_frf=0;
while feof(fid) ~= 1
    valec=fscanf(fid,'%s',[1]);
    if strcmp(key,valec) == 1
        nbr_frf=nbr_frf+1; % FRF counter
        for vv=1:5
            temp=fgetl(fid);
        end
        % RECORD 6
        valec=fscanf(fid,'%i',[4]);
        valec=fscanf(fid,'%s',[1]);
        resp_node=fscanf(fid,'%i',[1]);
        resp_dir=fscanf(fid,'%i',[1]);
        valec=fscanf(fid,'%s',[1]);
        reference_node=fscanf(fid,'%i',[1]);
        reference_dir=fscanf(fid,'%i',[1]);
        % RECORD 7
        precision=fscanf(fid,'%i',[1]);
        Nbr_data=fscanf(fid,'%i',[1]);
        abscissa_spacing=fscanf(fid,'%i',[1]);
        abscissa_minimum=fscanf(fid,'%g',[1]);
        abscissa_increment=fscanf(fid,'%g',[1]);
        valec=fscanf(fid,'%g',[1]);
        % RECORD 8
        valec=fscanf(fid,'%i',[4]);
        valec=fscanf(fid,'%s',[2]);
        % RECORD 9
        valec=fscanf(fid,'%i',[4]);
        valec=fscanf(fid,'%s',[2]);
        % RECORD 10
        valec=fscanf(fid,'%i',[4]);
        valec=fscanf(fid,'%s',[2]);
        % RECORD 11
        valec=fscanf(fid,'%i',[4]);
        valec=fscanf(fid,'%s',[2]);
        % RECORD 12
        if precision<=4; % real data
            even=1;
            temp=fscanf(fid,'%g',[Nbr_data*even]);
            N=length(temp);
            coh(:,nbr_frf)=temp(1:N);
        else % complex data
            even=2;
            temp=fscanf(fid,'%g',[Nbr_data*even]);
            N=length(temp);
            R=temp(1:2:N-1);
            I=temp(2:2:N);
            coh(:,nbr_frf)=complex(R,I);
        end
        abscissa_maximum=(Nbr_data-1)*abscissa_increment+abscissa_minimum;
        ff=[abscissa_minimum:abscissa_increment:abscissa_maximum];
        freq=ff';
        infoFRF(nbr_frf)=struct('response',resp_node,'dir_response',resp_dir,'excitation',reference_node,'dir_excitation',reference_dir);
    end
end