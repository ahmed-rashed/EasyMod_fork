function unv58write(H,numin,dir_excitation,numout,dir_response,fmin,finc,filename)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function writes the information of a FRF in a 58 type UFF file.
%
%  Synthax:
%  unv58write(H,numin,dir_excitation,numout,dir_response,fmin,finc,filename) 
%
%  Input data:
%  H: FRF vector,
%  numin: intput node DOF,
%  dir_excitation: intput node direction,
%  numout: output node DOF,
%  dir_response: output node direction,
%  fmin: start frequency,
%  finc: frequency resolution,
%  filename: file name where the information is written.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


if size(H,1) == 1
    H = H.' ;
end
N = size(H,1) ;
Fname = fopen(filename,'w') ;
% Data set delimiters
fprintf(Fname,'%6.0f\n',-1) ;
fprintf(Fname,'%6.0f\n',58) ;
% Record 1 (ID Line:  Model Identification)
fprintf(Fname,'FRF\n') ;
% Record 2 (ID Line:  Run Indentification)
fprintf(Fname,'NONE\n') ;
% Record 3 (ID Line:  Run Date and Time)
da = datestr(now) ;
fprintf(Fname,'%20.20s\n',da) ;
% Record 4 
fprintf(Fname,'NONE\n') ;
% Record 5 
fprintf(Fname,'NONE\n') ;
% Record 6 (DOF identification)
ft = 4 ;	% Function Type (4: for a frequency dependent function)
fi = 0 ;	% Function Identification Number
seqno = 0 ;	% Sequence Number
lcno = 0 ;	% Load Case Identification Number
resnam = 'NONE' ;	% Response Entity Name
resnod = numout ;	% Response node
resdof = dir_response ;	% Response Direction
refnam = 'NONE' ;	% Reference Entity Name
refnod = numin ;	% Reference node
refdof = dir_excitation ;	% Reference Direction
ss = blanks(82) ; 
ss(81:82) = '\n' ;
ss(01:5) = sprintf('%5.0f',ft) ;
ss(11:15) = sprintf('%5.0f',fi) ;
ss(16:20) = sprintf('%5.0f',seqno) ;
ss(21:30) = sprintf('%10.0f',lcno) ;
ss(31:35) = sprintf('%5.6s',resnam) ;
ss(42:51) = sprintf('%10.0f',resnod) ;
ss(52:55) = sprintf('%4.0f',resdof) ;
ss(56:60) = sprintf('%5.6s',refnam) ;
ss(67:76) = sprintf('%10.0f',refnod) ;
ss(77:80) = sprintf('%4.0f',refdof) ;
fprintf(Fname,ss) ;
% Record 7 (Data Form) 
odt = 5 ;	% Ordinate Data Type
ndp = N ;	% Number of Data (-pairs, if uneven abscissa ; -values if even abscissa)
abspac = 1 ;	% Abscissa Spacing
absmin = fmin ;	% Abscissa Minimum
absinc = finc ;	% Abscissa Increment
zaxis = 0.0 ;	% Z-axis Value
ss = blanks(80) ; 
ss(79:80) = '\n' ;
ss(01:10) = sprintf('%10.0f',odt) ;
ss(11:20) = sprintf('%10.0f',ndp) ;
ss(21:30) = sprintf('%10.0f',abspac) ;
ss(31:43) = sprintf('%13.4e',absmin) ;
ss(44:56) = sprintf('%13.4e',absinc) ;
ss(57:69) = sprintf('%13.4e',zaxis) ;
fprintf(Fname,ss) ;
% Record 8 (Abscissa Data Characteristics)
datatype = 18 ;	% Specific Data Type (18: for a frequency dependent function)
lue = 0 ;	% Length Units Exponent  
fue = 0 ;	% Force Units Exponent 
tue = 0 ;	% Temperature Units Exponent 
axislab = 'NONE' ;	% Axis Label
axisunlab = 'NONE' ;	% Axis Units Label
ss = blanks(80) ; 
ss(79:80) = '\n' ;
ss(01:10) = sprintf('%10.0f',datatype) ;
ss(11:15) = sprintf('%5.0f',lue) ;
ss(16:20) = sprintf('%5.0f',fue) ;
ss(21:25) = sprintf('%5.0f',tue) ;
ss(26:30) = sprintf('%5.6s',axislab) ;
ss(47:51) = sprintf('%5.6s',axisunlab) ;    
fprintf(Fname,ss) ;
% Record 9 (Ordinate Numerator Data Characteristics)
datatype = 12 ;	% Specific Data Type (13: for a frequency dependent function)
lue = 0 ;	% Length Units Exponent    (0: for force ; 1: for moment)
fue = 0 ;	% Force Units Exponent
tue = 0 ;	% Temperature Units Exponent
axislab = 'NONE' ;	% Axis Label
axisunlab = 'NONE' ;	% Axis Units Label
ss = blanks(80) ; 
ss(79:80) = '\n' ;
ss(01:10) = sprintf('%10.0f',datatype) ;
ss(11:15) = sprintf('%5.0f',lue) ;
ss(16:20) = sprintf('%5.0f',fue) ;
ss(21:25) = sprintf('%5.0f',tue) ;
ss(26:30) = sprintf('%5.6s',axislab) ;
ss(47:51) = sprintf('%5.6s',axisunlab) ;  
fprintf(Fname,ss) ;
% Record 10 (Ordinate Denominator Data Characteristics) 
datatype = 13 ;	% Specific Data Type (8: for a frequency dependent function ; 11 for  viscous damping)
lue = 0 ;	% Length Units Exponent      (0 for translational ; 1 for rotational displacement
fue = 0 ;	% Force Units Exponent
tue = 0 ;	% Temperature Units Exponent
axislab = 'NONE' ;	% Axis Label
axisunlab = 'NONE' ;	% Axis Units Label
ss = blanks(80) ; 
ss(79:80) = '\n' ;
ss(01:10) = sprintf('%10.0f',datatype) ;
ss(11:15) = sprintf('%5.0f',lue) ;
ss(16:20) = sprintf('%5.0f',fue) ;
ss(21:25) = sprintf('%5.0f',tue) ;
ss(26:30) = sprintf('%5.6s',axislab) ;
ss(47:51) = sprintf('%5.6s',axisunlab) ;  
fprintf(Fname,ss) ;
% Record 11 (Z-axis Data Characteristics)
datatype = 0 ;	% Specific Data Type (18: for a frequency dependent function)
lue = 0 ;	% Length Units Exponent
fue = 0 ;	% Force Units Exponent
tue = 0 ;	% Temperature Units Exponent
axislab = 'NONE' ;	% Axis Label
axisunlab = 'NONE' ;	% Axis Units Label
ss = blanks(80) ; 
ss(79:80) = '\n' ;
ss(01:10) = sprintf('%10.0f',datatype) ;
ss(11:15) = sprintf('%5.0f',lue) ;
ss(16:20) = sprintf('%5.0f',fue) ;
ss(21:25) = sprintf('%5.0f',tue) ;
ss(26:30) = sprintf('%5.6s',axislab) ;
ss(47:51) = sprintf('%5.6s',axisunlab) ; 
fprintf(Fname,ss) ;
% Record 12
R = real(H(:,1)) ;
I = imag(H(:,1)) ;
for i = 1:N
    fprintf(Fname,'%13.4e',R(i)) ;
    fprintf(Fname,'%13.4e',I(i)) ;
    if mod(i,3) == 0
        fprintf(Fname,'\n') ; 
    end
end
if mod(i,3) == 0
else
    fprintf(Fname,'\n') ;
end   
% Data Set Delimiter
fprintf(Fname,'%6.0f\n',-1) ; 
fclose(Fname) ;