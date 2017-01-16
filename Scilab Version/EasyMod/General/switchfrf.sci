function [Hswitch,freqswitch,infoFRFswitch] = switchfrf(H,frq,infoFRF)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function switches the FRF Hij(w)=Hji(w), according to the
//  Betti-Maxwell theorem.
//
//  Synthax:
//  [Hswitch,freqswitch,infoFRFswitch] = switchfrf(H,frq,infoFRF) 
//
//  Input data:
//  H: FRF matrix containing all the FRFs (column number = FRF number),
//  frq: frequency vector,
//  infoFRF: structure containing information on FRFs 
//            infoFRF(jj).response = jjth FRF response node
//            infoFRF(jj).dir_response = jjth FRF response direction (1=X, 2=Y,
//             3=Z, 4=RotX, 5=RotY, 6=RotZ)
//            infoFRF(jj).excitation = jjth FRF excitation node
//            infoFRF(jj).dir_excitation = jjth FRF excitation direction
//            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ). 
//
//  Output data:
//  Hswitch: switched FRF matrix containing all the FRFs (column number =
//  FRF number), 
//  freqswitch: frequency vector,
//  infoFRFswitch: structure containing information on switchedFRFs 
//            infoFRFswitch(jj).response = jjth FRF response node
//            infoFRFswitch(jj).dir_response = jjth FRF response direction
//            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ)
//            infoFRFswitch(jj).excitation = jjth FRF excitation node
//            infoFRFswitch(jj).dir_excitation = jjth FRF excitation direction
//            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ). 
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


[MM,NN] = size(H) ;
for index = 1:NN
    infoFRFswitch(index).excitation = infoFRF(index).response ;
    infoFRFswitch(index).response = infoFRF(index).excitation ;
    infoFRFswitch(index).dir_response = infoFRF(index).dir_response ;
    infoFRFswitch(index).dir_excitation = infoFRF(index).dir_excitation ;
end
Hswitch = H ;
freqswitch = frq ;

endfunction