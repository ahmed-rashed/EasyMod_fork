function [H,frq,infoFRF] = unv151read(filename)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function reads the 151 type UFF file in order to extract only the
//  information of all FRFs.
//
//  Synthax:
//  [H,frq,infoFRF] = unv151read(filename) 
//
//  Input data:
//  filename: file name where the information is readed,
//
//  Output data:
//  H: FRF matrix containing all the FRFs,
//  frq: frequency vector,
//  infoFRF: structure containing information on the FRFs 
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
nbr_frf = 0 ;
while meof(fid) == 0
    valec = mfscanf(fid,'%s') ;
    if isempty(valec) then
        break ;
    end
    if strcmp(FRF,valec) == 0
        nbr_frf = nbr_frf+1 ; // FRF counter
        for vv = 1:5
            temp = mgetl(fid,1) ;
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
        if precision <= 4 ; // real data
            even = 1 ;
            temp = [] ;
            for i = 1:(Nbr_data*even+1)
                temp = [temp mfscanf(fid,'%g')] ;
                i = i+1 ;
            end
            N = length(temp) ;
            H(:,nbr_frf) = temp(1:N) ;
        else // complex data
            even = 2 ;
            temp = [] ;
            for i = 1:(Nbr_data*even+1)
                temp = [temp mfscanf(fid,'%g')] ;
                i = i+1 ;
            end
            N = length(temp) ;
            R = temp(1:2:N-1)' ;
            I = temp(2:2:N)' ;
            H(:,nbr_frf) = complex(R,I) ;
        end
        abscissa_maximum = (Nbr_data-1)*abscissa_increment+abscissa_minimum ;
        ff = [abscissa_minimum:abscissa_increment:abscissa_maximum] ;
        frq = ff' ;
        infoFRF(nbr_frf) = struct('response',resp_node,'dir_response',resp_dir,'excitation',reference_node,'dir_excitation',reference_dir) ;
        valec = mfscanf(fid,'%s') ;
        valec = mfscanf(fid,'%s') ;
    end
end
mclose(fid) ;
infoFRF = infoFRF' ;

// Changing sign if excitation or response direction are negative
for ii = 1:size(infoFRF,1)
    if infoFRF(ii).dir_response < 0
        infoFRF(ii).dir_response = - infoFRF(ii).dir_response ;
        H(:,ii) = - H(:,ii) ;
    end
    if infoFRF(ii).dir_excitation < 0
        infoFRF(ii).dir_excitation = - infoFRF(ii).dir_excitation ;
        H(:,ii) = - H(:,ii) ;
    end
end

endfunction