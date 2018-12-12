function [H_cols,f_col,infoFRF]=unv58read(filename)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function reads the 58 type UFF file containing the information of a
%  FRF.
%
%  Synthax:
%  [H,freq,infoFRF]=unv58read(filename) 
%
%  Input data:
%  filename: file name where the information is readed,
%
%  Output data:
%  H: FRF vector containing the FRF,
%  freq: frequency vector,
%  infoFRF: structure containing information on the FRF 
%            infoFRF(jj).response=jjth FRF response node
%            infoFRF(jj).dir_response=jjth FRF response direction (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ)
%            infoFRF(jj).excitation=jjth FRF excitation node
%            infoFRF(jj).dir_excitation=jjth FRF excitation direction
%            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ). 
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS

f_id=fopen(filename,'r');
ss=['opening file ',filename];
disp(ss);
H_cols=[];
ii=1;
while feof(f_id) ~= 1
    valec=fscanf(f_id,'%s',[1]);
    if strcmp('FRF',valec)
        for vv=1:5
            fgetl(f_id);
        end
        % RECORD 6
        valec=fscanf(f_id,'%i',[4]);
        valec=fscanf(f_id,'%s',[1]);
        resp_node=fscanf(f_id,'%i',[1]);
        resp_dir=fscanf(f_id,'%i',[1]);
        valec=fscanf(f_id,'%s',[1]);
        reference_node=fscanf(f_id,'%i',[1]);
        reference_dir=fscanf(f_id,'%i',[1]);
        infoFRF(ii)=struct('response',resp_node,'dir_response',resp_dir,'excitation',reference_node,'dir_excitation',reference_dir);
        % RECORD 7
        precision=fscanf(f_id,'%i',[1]);
        N=fscanf(f_id,'%i',[1]);
        abscissa_spacing=fscanf(f_id,'%i',[1]);
        abscissa_minimum=fscanf(f_id,'%g',[1]);
        abscissa_increment=fscanf(f_id,'%g',[1]);
        valec=fscanf(f_id,'%g',[1]);
        % RECORD 8
        valec=fscanf(f_id,'%i',[4]);
        valec=fscanf(f_id,'%s',[2]);
        % RECORD 9
        valec=fscanf(f_id,'%i',[4]);
        valec=fscanf(f_id,'%s',[2]);
        % RECORD 10
        valec=fscanf(f_id,'%i',[4]);
        valec=fscanf(f_id,'%s',[2]);
        % RECORD 11
        valec=fscanf(f_id,'%i',[4]);
        valec=fscanf(f_id,'%s',[2]);
        % RECORD 12
        if precision <= 4
            even=1;
        else
            even=2;
        end
        temp=fscanf(f_id,'%g',N*even);
        
        if even == 1
            R=temp(1:N);
            I=0;
        else
            R=temp(1:even:even*N-1);
            I=temp(2:even:even*N);
        end

        H_col_temp=complex(R,I);
        % Changing sign if excitation or response direction are negative
        if infoFRF(ii).dir_response < 0
            infoFRF(ii).dir_response=-infoFRF(ii).dir_response;
            H_col_temp=-H_col_temp;
        end
        if infoFRF(ii).dir_excitation < 0
            infoFRF(ii).dir_excitation=-infoFRF(ii).dir_excitation;
            H_col_temp=-H_col_temp;
        end

        H_cols=[H_cols,H_col_temp];
        ii=ii+1;
    end
end
fclose(f_id);

f_col=(abscissa_minimum:abscissa_increment:(N-1)*abscissa_increment+abscissa_minimum)';
