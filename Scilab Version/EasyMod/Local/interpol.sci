function [f_interpol,theta_interpol] = interpol(f1,theta1,f2,theta2)

// ------------------   This file is part of EasyMod   ----------------------------
//  Internal function
//
//  Calculation of natural frequency by linear interpolation.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


f_interpol = (f2+f1)/2 ;
a = (theta2-theta1)/(f2-f1) ;
b = theta2-a*f2 ;
theta_interpol = a*f_interpol+b ;

endfunction