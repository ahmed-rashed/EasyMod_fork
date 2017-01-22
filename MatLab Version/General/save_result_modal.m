function infoMODE=save_result_modal(infoFRF)

% ------------------   This file is part of EasyMod   ----------------------------
%  User function
%
%  This function saves the information related to the modal parameters.
%
%  Synthax:
%  infoMODE=save_result_modal(infoFRF) 
%
%  Input data:
%  infoFRF: structure containing information on FRFs 
%            infoFRF(jj).response=jjth FRF response node
%            infoFRF(jj).dir_response=jjth FRF response direction (1=X, 2=Y,
%             3=Z, 4=RotX, 5=RotY, 6=RotZ)
%            infoFRF(jj).excitation=jjth FRF excitation node
%            infoFRF(jj).dir_excitation=jjth FRF excitation direction
%            (1=X, 2=Y, 3=Z, 4=RotX, 5=RotY, 6=RotZ) 
%            infoFRF(jj).infoMODE(index_mode).frequencyk=natural
%            frequency identified on the FRF(jj) for the mode 'index_mode'
%            infoFRF(jj).infoMODE(index_mode).etak=damping
%            identified on the FRF(jj) for the mode 'index_mode
%            infoFRF(jj).infoMODE(index_mode).Bijk=modal
%            constant identified on the FRF(jj) for the mode 'index_mode.
%
%  Output data:
%  infoMODE: structure containing the different identified parameters
%                infoMODE.frequencyk=natural frequency
%                infoMODE.etak=loss factor
%                infoMODE.Bijk=modal constant.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


[M,N]=size(infoFRF);
Nbr_frf=N;
[M,N]=size(infoFRF(1).infoMODE);
Nbr_mode=N;
Nbr_FRF_estime=ones(Nbr_mode)*Nbr_frf;
for i=1:Nbr_frf
    for j=1:Nbr_mode
        if isnan(infoFRF(i).infoMODE(j).frequencyk)
            infoFRF(i).infoMODE(j).frequencyk=0;
            infoFRF(i).infoMODE(j).etak=0;
            Nbr_FRF_estime(j)=Nbr_FRF_estime(j) - 1;
        end
    end
end
for index=1:Nbr_frf
     node_disp(index)=[infoFRF(index).response];
     disp_dir(index)=[infoFRF(index).dir_response];
end
node_fin=max(node_disp);
node=[1:1:node_fin];
Nbr_node=length(node);
natural_frequency=zeros(1,Nbr_mode);
visqueux_factor=zeros(1,Nbr_mode);
psi=zeros(3*Nbr_node,Nbr_mode);
for kk=1:Nbr_mode;
    for jj=1:Nbr_frf;
        temp_frequency=infoFRF(jj).infoMODE(kk).frequencyk;
        natural_frequency(kk)=natural_frequency(kk)+temp_frequency;
        temp_factor=infoFRF(jj).infoMODE(kk).etak;
        visqueux_factor(kk)=visqueux_factor(kk)+temp_factor;
        B=infoFRF(jj).infoMODE(kk).Bijk;
        psi((infoFRF(jj).response-1)*3+infoFRF(jj).dir_response,kk)=B;
    end
end
for i=1:Nbr_mode
   natural_frequency(i)=natural_frequency(i)/(Nbr_FRF_estime(i));
   visqueux_factor(i)=visqueux_factor(i)/(Nbr_FRF_estime(i));
end
% Normalization
for index=1:Nbr_mode
   temp=find(abs(psi(:,index))>0);
   psi(:,index)=psi(:,index)/psi(temp(1,1),index);
end
infoMODE=struct('frequencyk',natural_frequency,'etak',visqueux_factor,'Bijk',psi);
