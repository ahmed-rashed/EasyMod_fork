function [H,freq,infoFRF] = unv58read(filename)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function reads the 58 type UFF file containing the information of a
%  FRF.
%
%  Synthax:
%  [H,freq,infoFRF] = unv58read(filename) 
%
%  Input data:
%  filename: file name where the information is readed,
%
%  Output data:
%  H: FRF vector containing the FRF,
%  freq: frequency vector,
%  infoFRF: structure containing information on the FRF 
%            infoFRF(jj).response = jjth FRF response node
%            infoFRF(jj).dir_response = jjth FRF response direction (1=X, 2=Y,
%             3=Z, 4=RotX, 5=RotY, 6=RotZ)
%            infoFRF(jj).excitation = jjth FRF excitation node
%            infoFRF(jj).dir_excitation = jjth FRF excitation direction
%            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ). 
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


name_file = [filename] ;
fid = fopen(name_file,'r') ;
ss = ['opening file ',name_file] ;
disp(ss) ;
FRF = 'FRF' ;
NONE = 'NONE' ;
while feof(fid) ~= 1
    valec = fscanf(fid,'%s',[1]) ;
    if strcmp(FRF,valec) == 1 ;
        for vv = 1:5
            fgetl(fid) ;
        end
        % RECORD 6
        valec = fscanf(fid,'%i',[4]) ;
        valec = fscanf(fid,'%s',[1]) ;
        resp_node = fscanf(fid,'%i',[1]) ;
        resp_dir = fscanf(fid,'%i',[1]) ;
        valec = fscanf(fid,'%s',[1]) ;
        reference_node = fscanf(fid,'%i',[1]) ;
        reference_dir = fscanf(fid,'%i',[1]) ;
        % RECORD 7
        precision = fscanf(fid,'%i',[1]) ;
        Nbr_data = fscanf(fid,'%i',[1]) ;
        abscissa_spacing = fscanf(fid,'%i',[1]) ;
        abscissa_minimum = fscanf(fid,'%g',[1]) ;
        abscissa_increment = fscanf(fid,'%g',[1]) ;
        valec = fscanf(fid,'%g',[1]) ;
        % RECORD 8
        valec = fscanf(fid,'%i',[4]) ;
        valec = fscanf(fid,'%s',[2]) ;
        % RECORD 9
        valec = fscanf(fid,'%i',[4]) ;
        valec = fscanf(fid,'%s',[2]) ;
        % RECORD 10
        valec = fscanf(fid,'%i',[4]) ;
        valec = fscanf(fid,'%s',[2]) ;
        % RECORD 11
        valec = fscanf(fid,'%i',[4]) ;
        valec = fscanf(fid,'%s',[2]) ;
        % RECORD 12
        if precision <= 4
            even = 1 ;
        else
            even = 2 ;
        end
        temp = fscanf(fid,'%g',[Nbr_data*even]) ;  
    end       
end
fclose(fid) ;
N = length(temp) ;
if even == 1
    R = temp(1:N) ;
    I = 0 ;
else
    R = temp(1:even:N-1) ;
    I = temp(2:even:N) ;
end
H = complex(R,I) ;
ff = [abscissa_minimum:abscissa_increment:(Nbr_data-1)*abscissa_increment+abscissa_minimum] ;
freq = ff' ;
infoFRF = struct('response',resp_node,'dir_response',resp_dir,'excitation',reference_node,'dir_excitation',reference_dir) ;

% Changing sign if excitation or response direction are negative
if infoFRF.dir_response < 0
    infoFRF.dir_response = - infoFRF.dir_response ;
    H(:,ii) = - H(:,ii) ;
end
if infoFRF.dir_excitation < 0
    infoFRF.dir_excitation = - infoFRF.dir_excitation ;
    H(:,ii) = - H(:,ii) ;
end