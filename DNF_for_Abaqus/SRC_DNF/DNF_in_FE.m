function [AH,BH,G,H,a_ten,b_ten,r_ten,Omega,PHI,qa,qb,q]=DNF_in_FE(master_modes,dof_of_disp,dof_of_rotation,disp_applied,path,meshfile)


cd(path)

ifolder=1;
foldername=append('temp_files_',num2str(ifolder));
while exist(foldername,'dir')==7
    ifolder=ifolder+1;
    foldername=append('temp_files_',num2str(ifolder));    
end
mkdir(foldername)
mkdir(append(foldername,'/inp_files'))
mkdir(append(foldername,'/rpt_files'))


Current_path=append(path,foldername);
addpath(genpath(Current_path))  %path

cd(foldername)
%% compute M, K, PHI, Omega and exclude the clamped dof.
dof_per_node=dof_of_disp+dof_of_rotation;
Nmodes=length(master_modes);

[M,K,BCdofs,Nnodes]=Compute_MK(dof_per_node,meshfile);



[PHI,Omega]=Compute_PHI(master_modes,Nmodes,dof_of_disp,dof_of_rotation,Nnodes,Current_path,meshfile);




%% generating inp files for running the STEP

PHI_m=PHI;
PHI(BCdofs,:)=[];
PHI_q=PHI;

if dof_of_rotation==3
PHI_q(:,[(dof_of_disp+1):dof_per_node:end, (dof_of_disp+2):dof_per_node:end,(dof_of_disp+3):dof_per_node:end,])=[];
else if dof_of_rotation==2
PHI_q(:,[(dof_of_disp+1):dof_per_node:end, (dof_of_disp+2):dof_per_node:end])=[];
end
end

if dof_of_rotation==1
PHI_q(:,[(dof_of_disp+1):dof_per_node:end])=[];
end


for icount=1:Nmodes
    [maxphi(icount),maxidx(icount)]=max(abs(PHI_q(:,icount)));
    q(icount)=disp_applied/maxphi(icount);
end


generate_STEP_inp(PHI_m,q,dof_of_disp,dof_of_rotation,dof_per_node,Nnodes,Nmodes,Current_path,meshfile);

system('abaqus cae noGUI=inp_files/STEP_test');

%% Compute nonlinear coeffients



for iloop=1:Nmodes  %master_modes_selected

FL_p=K*PHI(:,iloop)*q(iloop);
mat_name=append("Mode_",num2str(iloop),"_2",".rpt");
FNL_p=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
TFL_p=(FNL_p-FL_p);


FL_n=-K*PHI(:,iloop)*q(iloop);
mat_name=append("Mode_",num2str(iloop),"_1",".rpt");
FNL_n=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
TFL_n=(FNL_n-FL_n);


Gten(:,iloop,iloop)=(TFL_p+TFL_n)/q(iloop)^2/2;
Hten(:,iloop,iloop,iloop)=(TFL_p-TFL_n)/q(iloop)^3/2;

       for vloop=1:Nmodes
G(vloop,iloop,iloop)=PHI(:,vloop)'*Gten(:,iloop,iloop);
H(vloop,iloop,iloop,iloop)=PHI(:,vloop)'*Hten(:,iloop,iloop,iloop);
       end 

end


for iloop=1:Nmodes
    for jloop=iloop+1:Nmodes



FL1=K*(PHI(:,iloop)*q(iloop)+PHI(:,jloop)*q(jloop));
mat_name=append("Mode_",num2str(iloop),"_Mode_",num2str(jloop),"_1.rpt");
FNL1=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
TFL1=(FNL1-FL1);



FL2=K*(PHI(:,iloop)*q(iloop)-PHI(:,jloop)*q(jloop));
mat_name=append("Mode_",num2str(iloop),"_Mode_",num2str(jloop),"_2.rpt");
FNL2=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
TFL2=(FNL2-FL2);



FL3=K*(-PHI(:,iloop)*q(iloop)-PHI(:,jloop)*q(jloop));
mat_name=append("Mode_",num2str(iloop),"_Mode_",num2str(jloop),"_3.rpt");
FNL3=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
TFL3=(FNL3-FL3);



Gten(:,iloop,jloop)=(TFL1+TFL3-2*Gten(:,iloop,iloop)*q(iloop)^2-2*Gten(:,jloop,jloop)*q(jloop)^2)/(4*q(iloop)*q(jloop));
Gten(:,jloop,iloop)=Gten(:,iloop,jloop);


for vloop=1:Nmodes
G(vloop,iloop,jloop)=PHI(:,vloop)'*Gten(:,iloop,jloop);
G(vloop,jloop,iloop)=G(vloop,iloop,jloop);

H(vloop,iloop,iloop,jloop)=((PHI(:,vloop)'*-(TFL2+TFL3))+2*G(vloop,iloop,iloop)*q(iloop)^2+2*G(vloop,jloop,jloop)*q(jloop)^2-2*H(vloop,jloop,jloop,jloop)*q(jloop)^3)/(6*q(iloop)^2*q(jloop));%%%%%%%%%%
H(vloop,iloop,jloop,iloop)=H(vloop,iloop,iloop,jloop);
H(vloop,jloop,iloop,iloop)=H(vloop,iloop,iloop,jloop);


H(vloop,iloop,jloop,jloop)=((PHI(:,vloop)'*(TFL1+TFL2))-2*G(vloop,iloop,iloop)*q(iloop)^2-2*G(vloop,jloop,jloop)*q(jloop)^2-2*H(vloop,iloop,iloop,iloop)*q(iloop)^3)/(6*q(iloop)*q(jloop)^2);%%%%%%%%%
H(vloop,jloop,iloop,jloop)=H(vloop,iloop,jloop,jloop);
H(vloop,jloop,jloop,iloop)=H(vloop,iloop,jloop,jloop);

end


    end
end



for vloop=1:1:Nmodes
for iloop=1:1:Nmodes
    for jloop=iloop+1:1:Nmodes
        for kloop=jloop+1:1:Nmodes

FL=K*(PHI(:,iloop)*q(iloop)+PHI(:,jloop)*q(jloop)+PHI(:,kloop)*q(kloop));
mat_name=append("Mode_",num2str(iloop),"_Mode_",num2str(jloop),"_Mode_",num2str(kloop),".rpt");
FNL=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
TFL=(FNL-FL);


H(vloop,iloop,jloop,kloop)=(PHI(:,vloop)'*TFL-G(vloop,iloop,iloop)*q(iloop)^2-G(vloop,jloop,jloop)*q(jloop)^2-G(vloop,kloop,kloop)*q(kloop)*q(kloop)...
-2*G(vloop,iloop,jloop)*q(iloop)*q(jloop)-2*G(vloop,iloop,kloop)*q(iloop)*q(kloop)-2*G(vloop,jloop,kloop)*q(jloop)*q(kloop)...
-H(vloop,iloop,iloop,iloop)*q(iloop)^3-H(vloop,jloop,jloop,jloop)*q(jloop)^3-H(vloop,kloop,kloop,kloop)*q(kloop)^3-3*H(vloop,iloop,iloop,jloop)*q(iloop)^2*q(jloop)...
-3*H(vloop,iloop,jloop,jloop)*q(iloop)*q(jloop)^2-3*H(vloop,iloop,kloop,kloop)*q(iloop)*q(kloop)^2-3*H(vloop,iloop,iloop,kloop)*q(iloop)^2*q(kloop)...
-3*H(vloop,jloop,kloop,kloop)*q(jloop)*q(kloop)^2-3*H(vloop,jloop,jloop,kloop)*q(jloop)^2*q(kloop))/q(iloop)/q(iloop)/q(kloop)/6;
        end
    end
end
end





for iloop=1:Nmodes    
    for jloop=iloop:Nmodes
   Zs(:,iloop,jloop)=((Omega(iloop)+Omega(jloop))^2*M-K)\Gten(:,iloop,jloop);
   Zd(:,iloop,jloop)=((-Omega(iloop)+Omega(jloop))^2*M-K)\Gten(:,iloop,jloop);   
   a_ten(:,iloop,jloop)=1/2*(Zd(:,iloop,jloop)+Zs(:,iloop,jloop));
   b_ten(:,iloop,jloop)=1/(2*Omega(iloop)*Omega(jloop))*(Zd(:,iloop,jloop)-Zs(:,iloop,jloop));
   r_ten(:,iloop,jloop)=(Omega(jloop)-Omega(iloop))/Omega(jloop)*Zd(:,iloop,jloop)+(Omega(jloop)+Omega(iloop))/Omega(jloop)*Zs(:,iloop,jloop); 


   if iloop~=jloop
   a_ten(:,jloop,iloop)= a_ten(:,iloop,jloop);
   b_ten(:,jloop,iloop)=b_ten(:,iloop,jloop);
   r_ten(:,jloop,iloop)=(Omega(iloop)-Omega(jloop))/Omega(iloop)*Zd(:,iloop,jloop)+(Omega(iloop)+Omega(jloop))/Omega(iloop)*Zs(:,iloop,jloop); 

    end
    end
end







%% Generating tensors and running ABAQUS


[qa,qb]=generate_DNF_inp(PHI_m,q,a_ten,b_ten,dof_of_disp,dof_of_rotation,dof_per_node,Nnodes,Nmodes,BCdofs,Current_path,disp_applied,meshfile);







system('abaqus cae noGUI=inp_files/DNF_test');

%% 
Tensortype={'a','b','s','d'};

for aorb=1:2



for iloop=1:Nmodes  %master_modes_selected
    for jloop=iloop:Nmodes  %master_modes_selected

mat_name=append(Tensortype{aorb},"_ten_",num2str(iloop),num2str(jloop),"_2",".rpt");
FNL1=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);


mat_name=append(Tensortype{aorb},"_ten_",num2str(iloop),num2str(jloop),"_1",".rpt");
FNL2=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);

for vloop=1:Nmodes  %master_modes_selected

if aorb==1
        GA(vloop,iloop,jloop)=PHI(:,vloop)'*(FNL1+FNL2)/(qa(iloop,jloop))^4/2;
    else
        GB(vloop,iloop,jloop)=PHI(:,vloop)'*(FNL1+FNL2)/(qb(iloop,jloop))^4/2;
end


end
    end
end



for iloop=1:1:Nmodes
    for jloop=1:1:Nmodes
        for kloop=jloop:1:Nmodes


mat_name=append("Phi_",num2str(iloop),'_',Tensortype{aorb},"_ten_",num2str(jloop),num2str(kloop),"_1",".rpt");
FNL1=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);

mat_name=append("Phi_",num2str(iloop),'_',Tensortype{aorb},"_ten_",num2str(jloop),num2str(kloop),"_2",".rpt");
FNL2=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);

for vloop=1:1:Nmodes
if aorb==1
AA(vloop,iloop,jloop,kloop)=(PHI(:,vloop)'*(FNL1+FNL2)-2*G(vloop,iloop,iloop)*(q(iloop))^2-2*GA(vloop,jloop,kloop)*qa(jloop,kloop)^4)/(4*q(iloop)*qa(jloop,kloop)^2);
else
BB(vloop,iloop,jloop,kloop)=(PHI(:,vloop)'*(FNL1+FNL2)-2*G(vloop,iloop,iloop)*(q(iloop))^2-2*GB(vloop,jloop,kloop)*qb(jloop,kloop)^4)/(4*q(iloop)*qb(jloop,kloop)^2);
end
end
        end
    end
end


end

for vloop=1:1:Nmodes
for iloop=1:1:Nmodes
    for jloop=1:1:Nmodes
        for kloop=jloop:1:Nmodes    
AH(vloop,iloop,jloop,kloop)=2*AA(vloop,iloop,jloop,kloop)+H(vloop,iloop,jloop,kloop);
BH(vloop,iloop,jloop,kloop)=2*BB(vloop,iloop,jloop,kloop);
if jloop~=kloop
AH(vloop,iloop,kloop,jloop)=AH(vloop,iloop,jloop,kloop);
BH(vloop,iloop,kloop,jloop)=BH(vloop,iloop,jloop,kloop);
end
        end
    end
end
end



end