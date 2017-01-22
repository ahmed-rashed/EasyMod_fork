function unv15and82write(Nodes,Connexions,filename,NodeColor,LineColor)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function records the geometry information (nodes and connexions) in
%  UFF file bringing together the file type 15 (for the nodes) and the file
%  type 82 (for the connexions).
%
%  Synthax:
%  unv15and82write(Nodes,Connexions,filename,NodeColor,LineColor) 
%
%  Input data:
%  Nodes: (Nx4) matrix defining the different nodes of the studied
%  structure (first column: node number, 2nd-4th columns: x,y,z coordinates
%  of the corresponding node)
%  Connexions: vector defining the connexions between the nodes
%  (=succession of nodes),
%  filename: file name where the information is saved,
%  NodeColor: node color,
%  LineColor: connexion color.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


N=size(Nodes,1);
if size(Nodes,2) ~= 4
    error('Nodes matrix uncorrect');
end
Connexions=Connexions(:);
C=length(Connexions);

Fname= fopen(filename,'w');

% Data Set Delimiters
fprintf(Fname,'%6.0f\n',-1); 
fprintf(Fname,'%6.0f\n',15);
   
% Record  - Fields 1 to 7
for i=1:N
    fprintf(Fname,'%10i',Nodes(i,1));
    fprintf(Fname,'%10i',0);
    fprintf(Fname,'%10i',0);
    fprintf(Fname,'%10i',NodeColor);
    fprintf(Fname,'%13.5e',Nodes(i,2));
    fprintf(Fname,'%13.5e',Nodes(i,3));
    fprintf(Fname,'%13.5e',Nodes(i,4));
    fprintf(Fname,'\n'); 
end
   
% Data Set Delimiter
fprintf(Fname,'%6.0f\n',-1); 

% Data Set Delimiters
fprintf(Fname,'%6.0f\n',-1); 
fprintf(Fname,'%6.0f\n',82);
 
% Record 1
fprintf(Fname,'%10i',1);
fprintf(Fname,'%10i',C);
fprintf(Fname,'%10i\n',LineColor);

% Record 2
fprintf(Fname,'%1.80s\n','NONE');
   
% Record  - Fields 1 to 7
for i=1:C
    fprintf(Fname,'%10i',Connexions(i));
    if mod(i,8) == 0
         fprintf(Fname,'\n'); 
    end
end
if mod(C,8) ~= 0
    fprintf(Fname,'\n');
end
   
% Data Set Delimiter
fprintf(Fname,'%6.0f\n',-1); 

fclose(Fname);