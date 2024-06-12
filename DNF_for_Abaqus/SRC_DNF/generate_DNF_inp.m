function [qa,qb]=generate_DNF_inp(PHI_m,q,a_ten,b_ten,dof_of_disp,dof_of_rotation,dof_per_node,Nnodes,Nmodes,BCdofs,path,disp_applied,meshfile)




for i=1:length(BCdofs)
            a_ten=[a_ten(1:BCdofs(i)-1,:,:);zeros(1,Nmodes,Nmodes);a_ten(BCdofs(i):end,:,:)];
            b_ten=[b_ten(1:BCdofs(i)-1,:,:);zeros(1,Nmodes,Nmodes);b_ten(BCdofs(i):end,:,:)];
end


a_ten_q=a_ten;
b_ten_q=b_ten;

if dof_of_rotation==3
a_ten_q([(dof_of_disp+1):dof_per_node:end, (dof_of_disp+2):dof_per_node:end,(dof_of_disp+3):dof_per_node:end,],:,:)=[];
else if dof_of_rotation==2
a_ten_q([(dof_of_disp+1):dof_per_node:end, (dof_of_disp+2):dof_per_node:end],:,:)=[];
end
end

if dof_of_rotation==1
a_ten_q([(dof_of_disp+1):dof_per_node:end],:,:)=[];
end


if dof_of_rotation==3
b_ten_q([(dof_of_disp+1):dof_per_node:end, (dof_of_disp+2):dof_per_node:end,(dof_of_disp+3):dof_per_node:end,],:,:)=[];
else if dof_of_rotation==2
b_ten_q([(dof_of_disp+1):dof_per_node:end, (dof_of_disp+2):dof_per_node:end],:,:)=[];
end
end

if dof_of_rotation==1
b_ten_q([(dof_of_disp+1):dof_per_node:end],:,:)=[];
end

for icount=1:Nmodes
    for jcount=1:Nmodes
    maxa=max(abs(a_ten_q(:,icount,jcount)));
    qa(icount,jcount)=sqrt(disp_applied/maxa)*0.1;
    maxb=max(abs(b_ten_q(:,icount,jcount)));
    qb(icount,jcount)=sqrt(disp_applied/maxb)*0.1;
    end
end



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
chr=append('*Nset, nset=DNFNode',num2str(T),',instance=',Istname);
chr=[chr newline num2str(T) ',' newline ];
Nset=append(Nset,chr);
end
fprintf(fileID,Nset);

fclose(fileID);

%%


Tensortype={'a','b','s','d'};

for aorb=1:2

for jloop=1:Nmodes  %master_modes_selected
    for kloop=jloop:Nmodes  %master_modes_selected   
        for ne_or_po=[1,2]


if aorb==1
Dsp=a_ten(:,jloop,kloop)*(-1)^ne_or_po*qa(jloop,kloop)^2;
else
Dsp=b_ten(:,jloop,kloop)*(-1)^ne_or_po*qb(jloop,kloop)^2;
end

file_name=append('inp_files/',Tensortype{aorb},'_ten_',num2str(jloop),num2str(kloop),'_',num2str(ne_or_po),'.inp');

copyfile("inp_files/raw.inp",file_name)

fileID = fopen(file_name,'a');

text='';
text=append(newline,text,'** STEP: Step-2',newline,'*Step, name=Step-DNF, nlgeom=YES',newline,'*Static',newline,'0.1, 1., 1e-10, 0.1',newline,...
            '*Solution Technique, type=QUASI-NEWTON, reform kernel=10',...
            newline,'** BOUNDARY CONDITIONS',newline);

for i=1:Nnodes
    Node=i;

Boundary='';
F(Node,:)=Dsp((dof_per_node*(i-1)+1):(dof_per_node*i),1);
for dof=1:dof_per_node
    chr=append('DNFNode',num2str(Node),', ',num2str(dof),', ',num2str(dof),', ',num2str(F(i,dof),9),newline);%ABAQUS only takes six digitals.
    Boundary=append(Boundary,chr);
end
text=append(text,'** Name: DNFBC-',num2str(Node),' Type: Displacement/Rotation',newline,'*Boundary',newline,Boundary);


end


text=append(text,newline,'**OUTPUT REQUESTS',newline,'*Restart, write, frequency=0',newline,'*Output, field, variable=PRESELECT',...
    newline,'*Output, history, variable=PRESELECT',newline,'*End Step');

fprintf(fileID,text);

fclose(fileID);

        end

    end

end
end





for aorb=1:2

for iloop=1:Nmodes


for jloop=1:Nmodes  %master_modes_selected
    
for kloop=jloop:Nmodes

for Dsp_case=1:2

if aorb==1
    Dsp=(-1)^Dsp_case*(PHI_m(:,iloop)*q(iloop)+a_ten(:,jloop,kloop)*qa(jloop,kloop)^2);
else
    Dsp=(-1)^Dsp_case*(PHI_m(:,iloop)*q(iloop)+b_ten(:,jloop,kloop)*qb(jloop,kloop)^2);
end

file_name2=append('inp_files/Phi_',num2str(iloop),'_',Tensortype{aorb},'_ten_',num2str(jloop),num2str(kloop),'_',num2str(Dsp_case),'.inp');

copyfile("inp_files/raw.inp",file_name2)
fileID2 = fopen(file_name2,'a');



text='';
text=append(newline,text,'** STEP: Step-2',newline,'*Step, name=Step-DNF, nlgeom=YES',newline,'*Static',newline,'0.1, 1., 1e-10, 0.1',newline,...
    '*Solution Technique, type=QUASI-NEWTON, reform kernel=10',...
    newline,'** BOUNDARY CONDITIONS',newline);

for i=1:Nnodes
    Node=i;

Boundary='';
F(Node,:)=Dsp((dof_per_node*(i-1)+1):(dof_per_node*i),1);
for dof=1:dof_per_node
    chr=append('DNFNode',num2str(Node),', ',num2str(dof),', ',num2str(dof),', ',num2str(F(i,dof),9),newline); %ABAQUS only takes six digitals.
    Boundary=append(Boundary,chr);
end
text=append(text,'** Name: DNFBC-',num2str(Node),' Type: Displacement/Rotation',newline,'*Boundary',newline,Boundary);

end



text=append(text,newline,'**OUTPUT REQUESTS',newline,'*Restart, write, frequency=0',newline,'*Output, field, variable=PRESELECT',...
    newline,'*Output, history, variable=PRESELECT',newline,'*End Step');

fprintf(fileID2,text);
fclose(fileID2);

end

end
end
end
end

textpath=append('"',path);

file_name='inp_files/DNF_test.py';
copyfile("../SRC_DNF/DNF_full.py",file_name)
fileID = fopen(file_name,'r+');
fgetl(fileID);
textpy=['import os',newline,'os.chdir(r',textpath,'")',newline,'Nmodes=',num2str(Nmodes),newline,'dof_RF=',num2str(dof_of_disp),newline,'dof_RM=',num2str(dof_of_rotation),newline,'Current_path=',textpath,'/"',newline];
fprintf(fileID,textpy);
fclose(fileID);


end
