function infoMODE = unv55read(filename,No)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function reads the 55 type UFF file containing the information
//  about the modal parameters.
//
//  Synthax:
//  infoMODE = unv55read(filename,No) 
//
//  Input data:
//  filename: file name where the information is readed,
//  No: number of experimental nodes.
//
//  Output data:
//  infoMODE: structure containing the different identified parameters
//                infoMODE.frequencyk = natural frequency
//                infoMODE.etak = loss factor
//                infoMODE.Bijk = modal constant.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


nom_unv = [filename] ;
frequency = 'FREQUENCY' ;
frq = [] ;
amort = [] ;
temp = [] ;
tempo = [] ;
modal = [] ;
fichier = mopen(nom_unv,'r') ;
num = 0 ;
valec = mfscanf(fichier,'%s') ;
while meof(fichier) == 0
    if strcmp(valec,'Residuals') == 0 
        break         
    end
    if strcmp(frequency,valec) == 0
        num = num+1 ; // counter
        frq0 = mfscanf(fichier,'%g') ;
        frq(num) = frq0 ; // natural frequency storing
        valec = mfscanf(fichier,'%s') ;
        valec = mfscanf(fichier,'%s') ;
        temp = mfscanf(fichier,'%g') ;
        xi(num) = temp ; // damping ratio storing
        valec = mfscanf(fichier,'%i') ;
        valec = mfscanf(fichier,'%i') ;
        valec = mfscanf(fichier,'%i') ;
        valec = mfscanf(fichier,'%i') ;
        valec = mfscanf(fichier,'%i') ;
        valec = mfscanf(fichier,'%i') ;
        valec = mfscanf(fichier,'%i') ;
        valec = mfscanf(fichier,'%i') ;
        ind_cplx = mfscanf(fichier,'%i') ;
        NDV = mfscanf(fichier,'%i') ; // number of data per node
        if ind_cplx == 2
            for i=1:8
                valec = mfscanf(fichier,'%g') ;
            end
            Nbr = No*(NDV+1) ;
            temp = [] ;
            for i=1:Nbr
                temp = [temp mfscanf(fichier,'%g')] ;
            end            
            [M,N] = size(temp) ;
            for index = 1:NDV+1:N-NDV
                tempo = [tempo temp(index+1:index+NDV)] ;
            end
            psi(num,:) = tempo ;
            tempo = [] ;
        elseif ind_cplx == 5
            for i=1:10
                valec = mfscanf(fichier,'%g') ;
            end
            Nbr = No*(2*NDV+1) ;
            temp = [] ;
            for i=1:Nbr
                temp = [temp mfscanf(fichier,'%g')] ;
            end
            [M,N] = size(temp) ;
            for index = 1:2*NDV+1:N-2*NDV
                x = complex(temp(index+1),temp(index+2)) ;
                y = complex(temp(index+3),temp(index+4)) ;
                z = complex(temp(index+5),temp(index+6)) ;
                xyz = [x y z] ;
                tempo = [tempo xyz] ;
            end
            psi(num,:) = tempo ;
            tempo = [] ;
        end
    end
    valec = mfscanf(fichier,'%s') ;
end
mclose(fichier) ;
frq = frq' ;
xi = xi' ;
psi = psi.' ;
infoMODE = struct('frequencyk',frq,'etak',xi*2,'Bijk',psi) ;

endfunction