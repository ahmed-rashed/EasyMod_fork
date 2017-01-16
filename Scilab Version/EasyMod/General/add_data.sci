function infoFRF = add_data(index_FRF,frequencyk,etak,Bijk,infoFRF,index_mode)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  Storage of modal parameters for each FRF and for each DOF.
//
//  Synthax:
//  infoFRF = add_data(index_FRF,frequencyk,etak,Bijk,infoFRF,index_mode)
//
//  Input data:
//  index_FRF: index of studied FRF,
//  frequencyk: natural frequency,
//  etak: loss factor,
//  Bijk: modal constant,
//  infoFRF: structure containing information on FRFs 
//            infoFRF(jj).response = jjth FRF response node
//            infoFRF(jj).dir_response = jjth FRF response direction (1=X, 2=Y,
//             3=Z, 4=RotX, 5=RotY, 6=RotZ)
//            infoFRF(jj).excitation = jjth FRF excitation node
//            infoFRF(jj).dir_excitation = jjth FRF excitation direction
//            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ), 
//  index_mode: index of studied mode.
//
//  Output data:
//  infoFRF: initial structure where modal parameters are added
//   modaux identifiés
//      FRF(index_FRF).infoMODE(index_mode,1) = natural frequency,
//      FRF(index_FRF).infoMODE(index_mode,2) = loss factor,
//      FRF(index_FRF).infoMODE(index_mode,3) = modal constant.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


infoFRF(index_FRF).infoMODE(index_mode,1) = frequencyk ;
infoFRF(index_FRF).infoMODE(index_mode,2) = etak ;
infoFRF(index_FRF).infoMODE(index_mode,3) = Bijk ;

endfunction