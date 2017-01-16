function [H,frq,infoFRF] = unv58read(filename)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function reads the 58 type UFF file containing the information of a
//  FRF.
//
//  Synthax:
//  [H,frq,infoFRF] = unv58read(filename) 
//
//  Input data:
//  filename: file name where the information is readed,
//
//  Output data:
//  H: FRF vector containing the FRF,
//  frq: frequency vector,
//  infoFRF: structure containing information on the FRF 
//            infoFRF(jj).response = jjth FRF response node
//            infoFRF(jj).dir_response = jjth FRF response direction (1=X, 2=Y,
//             3=Z, 4=RotX, 5=RotY, 6=RotZ)
//            infoFRF(jj).excitation = jjth FRF excitation node
//            infoFRF(jj).dir_excitation = jjth FRF excitation direction
//            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ). 
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


name_file = [filename] ;
fid = mopen(name_file,'r') ;
ss = ['opening file ',name_file] ;
disp(ss) ;
FRF = 'FRF' ;
NONE = 'NONE' ;
while meof(fid) == 0 ,
    valec = mfscanf(fid,'%s') ;
    if strcmp(valec,FRF) == 0 
        for vv = 1:5
            mgetl(fid,1) ;
        end
        // RECORD 6
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%s') ;
        resp_node = mfscanf(fid,'%i') ;
        resp_dir = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%s') ;
        reference_node = mfscanf(fid,'%i') ;
        reference_dir = mfscanf(fid,'%i') ;
        // RECORD 7
        precision = mfscanf(fid,'%i') ; 
        Nbr_data = mfscanf(fid,'%i') ;
        abscissa_spacing = mfscanf(fid,'%i') ;
        abscissa_minimum = mfscanf(fid,'%g') ;
        abscissa_increment = mfscanf(fid,'%g') ;
        valec = mfscanf(fid,'%g') ;
        // RECORD 8
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%s') ;
        valec = mfscanf(fid,'%s') ;
        // RECORD 9
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%s') ;
        valec = mfscanf(fid,'%s') ;
        // RECORD 10
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%s') ;
        valec = mfscanf(fid,'%s') ;
        // RECORD 11
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%i') ;
        valec = mfscanf(fid,'%s') ;
        valec = mfscanf(fid,'%s') ;
        // RECORD 12
        if precision <= 4
            even = 1 ;
        else
            even = 2 ;
        end
        temp = [] ;
        for i = 1:(Nbr_data*even+1)
            temp = [temp mfscanf(fid,'%g')] ; 
            i = i+1 ;
        end
        valec = mfscanf(fid,'%s') ;
        valec = mfscanf(fid,'%s') ;
    end
end
mclose(fid) ;
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
frq = ff' ;
H = H.' ;
infoFRF = struct('response',resp_node,'dir_response',resp_dir,'excitation',reference_node,'dir_excitation',reference_dir) ;

// Changing sign if excitation or response direction are negative
if infoFRF.dir_response < 0
    infoFRF.dir_response = - infoFRF.dir_response ;
    H(:,ii) = - H(:,ii) ;
end
if infoFRF.dir_excitation < 0
    infoFRF.dir_excitation = - infoFRF.dir_excitation ;
    H(:,ii) = - H(:,ii) ;
end

endfunction