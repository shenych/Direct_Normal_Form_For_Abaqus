function generate_STEP_inp(PHI_m,q,dof_of_disp,dof_of_rotation,dof_per_node,Nnodes,Nmodes,path,meshfile)



fileID0 = fopen(meshfile);
OA=textscan(fileID0,'%s');
fclose(fileID0);

Index = find(contains(OA{1,1},'*Instance'));
Index2 = find(contains(OA{1,1},'name='));

for i=1:length(Index2)
if Index2(i)>Index
instance_name=OA{1,1}(Index2(i));
break
end
end
Istname=num2str(cell2mat(instance_name));
Istname=Istname(6:end-1);





file_name=append('inp_files/raw.inp');

copyfile(meshfile,file_name)

fileID = fopen(file_name,'a');

Nset='';







for T=1:Nnodes %number of nodes
    chr=append('*Nset, nset=STEPNode',num2str(T),',instance=',Istname);%%%%%%%%%%%%%
    chr=[chr newline num2str(T) ',' newline ];
    Nset=append(Nset,chr);
end

fprintf(fileID,Nset);

fclose(fileID);



for iloop=1:Nmodes  %master_modes_selected

    for ne_or_po=[1,2]

        Dsp=PHI_m(:,iloop)*(-1)^ne_or_po*q(iloop);

        file_name=append('inp_files/Mode_',num2str(iloop),'_',num2str(ne_or_po),'.inp');

        copyfile("inp_files/raw.inp",file_name)

        fileID = fopen(file_name,'a');


        text='';
        text=append(newline,text,'** STEP: Step-2',newline,'*Step, name=Step-STEP, nlgeom=YES',newline,'*Static',newline,'0.1, 1., 1e-05, 0.1',newline,...
            '*Solution Technique, type=QUASI-NEWTON, reform kernel=20',...
            newline,'** BOUNDARY CONDITIONS',newline);

        for i=1:Nnodes
            Node=i;

                Boundary='';
                F(Node,:)=Dsp((dof_per_node*(i-1)+1):(dof_per_node*i),1);
                for dof=1:dof_per_node
                    chr=append('STEPNode',num2str(Node),', ',num2str(dof),', ',num2str(dof),', ',num2str(F(i,dof),6),newline);
                    Boundary=append(Boundary,chr);
                end
                text=append(text,'** Name: STEPBC-',num2str(Node),' Type: Displacement/Rotation',newline,'*Boundary',newline,Boundary);

        end


        text=append(text,newline,'**OUTPUT REQUESTS',newline,'*Restart, write, frequency=0',newline,'*Output, field, variable=PRESELECT',...
            newline,'*Output, history, variable=PRESELECT',newline,'*End Step');


        fprintf(fileID,text);

        fclose(fileID);

    end

end



for iloop=1:Nmodes  %master_modes_selected
    for jloop=iloop+1:Nmodes


        for Dsp_case=1:3

            if Dsp_case==1
                Dsp=PHI_m(:,iloop)*q(iloop)+PHI_m(:,jloop)*q(jloop);
            else
                if Dsp_case==2
                    Dsp=PHI_m(:,iloop)*q(iloop)-PHI_m(:,jloop)*q(jloop);
                else
                    Dsp=-PHI_m(:,iloop)*q(iloop)-PHI_m(:,jloop)*q(jloop);
                end
            end


            file_name=append('inp_files/Mode_',num2str(iloop),'_Mode_',num2str(jloop),'_',num2str(Dsp_case),'.inp');

            copyfile("inp_files/raw.inp",file_name)

            fileID = fopen(file_name,'a');


            text='';
            text=append(newline,text,'** STEP: Step-2',newline,'*Step, name=Step-STEP, nlgeom=YES',newline,'*Static',newline,'0.1, 1., 1e-05, 0.1',newline,'** BOUNDARY CONDITIONS',newline);

            for i=1:Nnodes
                Node=i;
                    Boundary='';
                    F(Node,:)=Dsp((dof_per_node*(i-1)+1):(dof_per_node*i),1);
                    for dof=1:dof_per_node
                        chr=append('STEPNode',num2str(Node),', ',num2str(dof),', ',num2str(dof),', ',num2str(F(i,dof),6),newline);
                        Boundary=append(Boundary,chr);
                    end
                    text=append(text,'** Name: STEPBC-',num2str(Node),' Type: Displacement/Rotation',newline,'*Boundary',newline,Boundary);
            end

            text=append(text,newline,'**OUTPUT REQUESTS',newline,'*Restart, write, frequency=0',newline,'*Output, field, variable=PRESELECT',...
                newline,'*Output, history, variable=PRESELECT',newline,'*End Step');

            fprintf(fileID,text);
            fclose(fileID);


        end
    end
end







for iloop=1:Nmodes  %master_modes_selected
    for jloop=iloop+1:Nmodes
        for kloop=jloop+1:Nmodes

            Dsp=PHI_m(:,iloop)*q(iloop)+PHI_m(:,jloop)*q(jloop)+PHI_m(:,kloop)*q(kloop);



            file_name=append('inp_files/Mode_',num2str(iloop),'_Mode_',num2str(jloop),'_Mode_',num2str(kloop),'.inp');

            copyfile("inp_files/raw.inp",file_name)

            fileID = fopen(file_name,'a');


            text='';
            text=append(newline,text,'** STEP: Step-2',newline,'*Step, name=Step-STEP, nlgeom=YES',newline,'*Static',newline,'0.1, 1., 1e-05, 0.1',newline,...
                          '*Solution Technique, type=QUASI-NEWTON, reform kernel=20',...
                newline,'** BOUNDARY CONDITIONS',newline);

            for i=1:Nnodes
                Node=i;
                    Boundary='';
                    F(Node,:)=Dsp((dof_per_node*(i-1)+1):(dof_per_node*i),1);
                    for dof=1:dof_per_node
                        chr=append('STEPNode',num2str(Node),', ',num2str(dof),', ',num2str(dof),', ',num2str(F(i,dof),6),newline);
                        Boundary=append(Boundary,chr);
                    end
                    text=append(text,'** Name: STEPBC-',num2str(Node),' Type: Displacement/Rotation',newline,'*Boundary',newline,Boundary);
            end

            text=append(text,newline,'**OUTPUT REQUESTS',newline,'*Restart, write, frequency=0',newline,'*Output, field, variable=PRESELECT',...
                newline,'*Output, history, variable=PRESELECT',newline,'*End Step');

            fprintf(fileID,text);
            fclose(fileID);


        end
    end
end






textpath=append('"',path);

file_name='inp_files/STEP_test.py';
copyfile("../SRC_DNF/STEP_full.py",file_name)
fileID = fopen(file_name,'r+');
fgetl(fileID);
textpy=['import os',newline,'os.chdir(r',textpath,'")',newline,'Nmodes=',num2str(Nmodes),newline,'dof_RF=',num2str(dof_of_disp),newline,'dof_RM=',num2str(dof_of_rotation),newline,'Current_path=',textpath,'/"',newline];
fprintf(fileID,textpy);
fclose(fileID);


end