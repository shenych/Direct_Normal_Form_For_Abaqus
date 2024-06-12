function [PHI,Omega]=Compute_PHI(master_modes,Nmodes,dof_of_disp,dof_of_rotation,Nnodes,path,meshfile)

file_name=append('inp_files/Compute_PHI.inp');

copyfile(meshfile,file_name)

fileID = fopen(file_name,'a');

text='';

text=append(text,'** STEP: Step-1',newline,... 
'*Step, name=linearmode, nlgeom=NO, perturbation',newline,...
'*Frequency, eigensolver=Lanczos, sim, acoustic coupling=on, normalization=mass',newline,...
'50, , , , ,',newline,...
',*Restart, write, frequency=0',newline,...
'*Output, field, variable=PRESELECT',newline,...
'*Output, history, variable=PRESELECT',newline,...
'*End Step');

fprintf(fileID,text);
fclose(fileID);




textpath=append('"',path);

file_name='inp_files/Compute_PHI_test.py';
copyfile("../SRC_DNF/Compute_PHI.py",file_name)
fileID = fopen(file_name,'r+');
fgetl(fileID);
textmodes='master_modes=[';
for ii=1:Nmodes
textmodes=append(textmodes,num2str(master_modes(ii)),',');
end
textmodes=append(textmodes,'];');
textpy=['import os',newline,'os.chdir(r',textpath,'")',newline,textmodes,'Nmodes=',num2str(Nmodes),newline,'dof_RF=',num2str(dof_of_disp),newline,'dof_RM=',num2str(dof_of_rotation),newline,'Current_path=',textpath,'/"',newline];
fprintf(fileID,textpy);
fclose(fileID);



system('abaqus cae noGUI=inp_files/Compute_PHI_test');
dof_per_node=dof_of_disp+dof_of_rotation;
for i=1:Nmodes
    mat_name=append('Compute_PHI',num2str(i),'.rpt');
    [Omega(i),PHI(:,i)]=read_phi(mat_name,dof_per_node,Nnodes);
end

end



function [Omega,PHI]=read_phi(file_name,dof_per_node,Nnodes)


            fileID = fopen(file_name);
            A=textscan(fileID,'%s');
            fclose(fileID);

            N_account=Nnodes;
            for i=1:length(A{1,1})
                t(i)=str2double(A{1,1}{i});
            end


            for i=1:length(t)
                if t(i)==1&&t(i+dof_per_node+1)==2&&t(i+2*(dof_per_node+1))==3
                    break
                end
            end
            %i: where meaningful data starts

            for j=1:length(t)
                if t(j)==N_account&&t(j-(dof_per_node+1))==N_account-1
                    break
                end
            end
            %j+dof_per_node: where meaningful data ends

            data_raw=t(i:(j+dof_per_node));


            kicount=1;
            for k=1:length(data_raw)
                if mod(k,(dof_per_node+1))~=1 %dof+1
                    PHI(kicount,1)=data_raw(k);
                    kicount=kicount+1;
                end
            end


Index = find(contains(A{1,1},'Freq'));

Omega=(str2num(cell2mat(A{1,1}(Index+2))))*2*pi;

end
