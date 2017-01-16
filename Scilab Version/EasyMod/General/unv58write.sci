function unv58write(H,numin,dir_excitation,numout,dir_response,fmin,finc,filename)

// ------------------   This file is part of EasyMod   ----------------------------
//  User function
//
//  This function writes the information of a FRF in a 58 type UFF file.
//
//  Synthax:
//  unv58write(H,numin,dir_excitation,numout,dir_response,fmin,finc,filename) 
//
//  Input data:
//  H: FRF vector,
//  numin: intput node DOF,
//  dir_excitation: intput node direction,
//  numout: output node DOF,
//  dir_response: output node direction,
//  fmin: start frequency,
//  finc: frequency resolution,
//  filename: file name where the information is written.
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


if size(H,1) == 1
    H = H.' ;
end
N = size(H,1) ;
Fname = mopen(filename,'w') ;
// Data set delimiters
fprintf(Fname,'%6.0f\n',-1) ;
fprintf(Fname,'%6.0f\n',58) ;
// Record 1 (ID Line:  Model Identification)
fprintf(Fname,'FRF\n') ;
// Record 2 (ID Line:  Run Indentification)
fprintf(Fname,'NONE\n') ;
// Record 3 (ID Line:  Run Date and Time)
da = date() ;
d = getdate() ;
fprintf(Fname,'%s %s:%s:%s\n',da,string(d(7)),string(d(8)),string(d(9))) ;
// Record 4 
fprintf(Fname,'NONE\n') ;
// Record 5 
fprintf(Fname,'NONE\n') ;
// Record 6 (DOF identification)
ft = 4 ;	// Function Type (4: for a frequency dependent function)
fi = 0 ;	// Function Identification Number
seqno = 0 ;	// Sequence Number
lcno = 0 ;	// Load Case Identification Number
resnam = 'NONE' ;	// Response Entity Name
resnod = numout ;	// Response node
resdof = dir_response ;	// Response Direction
refnam = 'NONE' ;	// Reference Entity Name
refnod = numin ;	// Reference node
refdof = dir_excitation ;	// Reference Direction
fprintf(Fname,'%5.0f%10.0f%5.0f%10.0f%5.6s%16.0f%4.0f%5.6s%16.0f%4.0f\n',...
ft,fi,seqno,lcno,resnam,resnod,resdof,refnam,refnod,refdof) ;
// Record 7 (Data Form) 
odt = 5 ;	// Ordinate Data Type
ndp = N ;	// Number of Data (-pairs, if uneven abscissa ; -values if even abscissa)
abspac = 1 ;	// Abscissa Spacing
absmin = fmin ;	// Abscissa Minimum
absinc = finc ;	// Abscissa Increment
zaxis = 0.0 ;	// Z-axis Value
fprintf(Fname,'%10.0f%10.0f%10.0f%13.4e%13.4e%13.4e\n',...
odt,ndp,abspac,absmin,absinc,zaxis) ;
// Record 8 (Abscissa Data Characteristics)
datatype = 18 ;	// Specific Data Type (18: for a frequency dependent function)
lue = 0 ;	// Length Units Exponent  
fue = 0 ;	// Force Units Exponent 
tue = 0 ;	// Temperature Units Exponent 
axislab = 'NONE' ;	// Axis Label
axisunlab = 'NONE' ;	// Axis Units Label
fprintf(Fname,'%10.0f%5.0f%5.0f%5.0f%5.6s%21.6s\n',datatype,lue,fue,tue,axislab,axisunlab) ;
// Record 9 (Ordinate Numerator Data Characteristics)
datatype = 12 ;	// Specific Data Type (13: for a frequency dependent function)
lue = 0 ;	// Length Units Exponent    (0: for force ; 1: for moment)
fue = 0 ;	// Force Units Exponent
tue = 0 ;	// Temperature Units Exponent
axislab = 'NONE' ;	// Axis Label
axisunlab = 'NONE' ;	// Axis Units Label

fprintf(Fname,'%10.0f%5.0f%5.0f%5.0f%5.6s%21.6s\n',datatype,lue,fue,tue,axislab,axisunlab) ;
// Record 10 (Ordinate Denominator Data Characteristics) 
datatype = 13 ;	// Specific Data Type (8: for a frequency dependent function ; 11 for  viscous damping)
lue = 0 ;	// Length Units Exponent      (0 for translational ; 1 for rotational displacement
fue = 0 ;	// Force Units Exponent
tue = 0 ;	// Temperature Units Exponent
axislab = 'NONE' ;	// Axis Label
axisunlab = 'NONE' ;	// Axis Units Label
fprintf(Fname,'%10.0f%5.0f%5.0f%5.0f%5.6s%21.6s\n',datatype,lue,fue,tue,axislab,axisunlab) ;
// Record 11 (Z-axis Data Characteristics)
datatype = 0 ;	// Specific Data Type (18: for a frequency dependent function)
lue = 0 ;	// Length Units Exponent
fue = 0 ;	// Force Units Exponent
tue = 0 ;	// Temperature Units Exponent
axislab = 'NONE' ;	// Axis Label
axisunlab = 'NONE' ;	// Axis Units Label
fprintf(Fname,'%10.0f%5.0f%5.0f%5.0f%5.6s%21.6s\n',datatype,lue,fue,tue,axislab,axisunlab) ;
// Record 12
R = real(H(:,1)) ;
I = imag(H(:,1)) ;
for i = 1:N
    fprintf(Fname,'%13.4e',R(i)) ;
    fprintf(Fname,'%13.4e',I(i)) ;
    if pmodulo(i,3) == 0
        fprintf(Fname,'\n') ; 
    end
end
if pmodulo(i,3) == 0
else
    fprintf(Fname,'\n') ;
end   
// Data Set Delimiter
fprintf(Fname,'%6.0f\n',-1) ; 
mclose(Fname) ;

endfunction