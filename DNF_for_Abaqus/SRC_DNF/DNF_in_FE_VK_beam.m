function [AH,BH,G,H,a_ten,b_ten,r_ten,Omega,PHI]=DNF_in_FE_VK_beam(master_modes,path)


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


Nmodes=length(master_modes);
%% compute M, K, PHI, Omega and exclude the clamped dof.

% [PHI,Omega,M,K,BCdofs,Nnodes]=Compute_MK(dof_per_node,meshfile);

% Omega=Omega(master_modes);
% PHI=PHI(:,master_modes);


%% Define Beam Properties, Create the FE model
addpath(genpath('C:\YichangShen\Matlab_code\vonKarman_3D\src'));

% ---- Material Properties ----
beam.material_E   = 210 * 1e9; %Pa
beam.material_G   = beam.material_E / (2.*(1.+0.3)); %Pa
beam.material_rho = 7850; %kg/m^3
% ---- Baseline cross section ----
% define either base and height for rectangular beam, or radius for circular
beam.cross_section_b = 1e-2; %m
beam.cross_section_h = 1e-2; %m
% Thickness coefficients for each component (not currently used, leave as [1])
beam.thick_coeff = [1];

% Read geometry and create mesh with specified maximum element length and boundary condition
[mesh] = mesh_generator(beam, 'C:\YichangShen\Matlab_code\vonKarman_3D\examples\Geometry\beam_1m_2D.mat', 1/30, [1,2], 'clamped' ,0);
% calculate beam profile properties
[beam] = profile_properties(mesh, beam);
% Form underlying linear system
[linearmodel] = LinearBeamModel(mesh, beam);







M=linearmodel.Mlin;  %mass, stiffness matrix and eigenfrequencies of the system
K=linearmodel.Klin;
Omega=linearmodel.natfreq(master_modes);
PHI=linearmodel.modeshapes(master_modes);


q=1e-3*ones(Nmodes,1);

%% generating inp files for running the STEP

% PHI_m=PHI;
% for i=1:length(BCdofs)
%     PHI_m=[PHI_m(1:BCdofs(i)-1,:);zeros(1,Nmodes);PHI_m(BCdofs(i):end,:)];
% end
% 
% for icount=1:Nmodes
%     [maxphi(icount),maxidx]=max(abs(PHI_m(thickness_dof:dof_per_node:end,icount)));
%     q(icount)=Disp_applied/maxphi(icount);
% end
% maxidx=maxidx*dof_per_node-thickness_dof;
% 
% 
% generate_STEP_inp(PHI_m,q,dof_per_node,Nnodes,Nmodes,Current_path,meshfile);
% 
% system('abaqus cae noGUI=inp_files/STEP_test');







%% Compute nonlinear coeffients



for iloop=1:Nmodes  %master_modes_selected

% FL_p(:,iloop)=K*PHI(:,iloop)*q(iloop);
% mat_name=append("Mode_",num2str(iloop),"_2",".rpt");
% FNL_p(:,iloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
% TFL_p(:,iloop)=(FNL_p(:,iloop)-FL_p(:,iloop));
% 
% 
% FL_n(:,iloop)=-K*PHI(:,iloop)*q(iloop);
% mat_name=append("Mode_",num2str(iloop),"_1",".rpt");
% FNL_n(:,iloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
% TFL_n(:,iloop)=(FNL_n(:,iloop)-FL_n(:,iloop));
% 

TFL_p=get_nonlinstiffness(mesh, beam, q(iloop)*PHI(:,iloop))*q(iloop)*PHI(:,iloop);
TFL_n=get_nonlinstiffness(mesh, beam, -q(iloop)*PHI(:,iloop))*-q(iloop)*PHI(:,iloop);
 

Gten(:,iloop,iloop)=(TFL_n+TFL_p)/q(iloop)^2/2;
Hten(:,iloop,iloop,iloop)=(TFL_p-TFL_n)/q(iloop)^3/2;

end



for vloop=1:Nmodes
        for iloop=1:Nmodes
G(vloop,iloop,iloop)=PHI(:,vloop)'*Gten(:,iloop,iloop);
H(vloop,iloop,iloop,iloop)=PHI(:,vloop)'*Hten(:,iloop,iloop,iloop);
        end
end




for iloop=1:Nmodes
    for jloop=iloop+1:Nmodes



% FL1(:,iloop,jloop)=K*(PHI(:,iloop)*q(iloop)+PHI(:,jloop)*q(jloop));
% mat_name=append("Mode_",num2str(iloop),"_Mode_",num2str(jloop),"_1.rpt");
% FNL1(:,iloop,jloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
% TFL1(:,iloop,jloop)=(FNL1(:,iloop,jloop)-FL1(:,iloop,jloop));
% 
% 
% 
% 
% FL2(:,iloop,jloop)=K*(PHI(:,iloop)*q(iloop)-PHI(:,jloop)*q(jloop));
% mat_name=append("Mode_",num2str(iloop),"_Mode_",num2str(jloop),"_2.rpt");
% FNL2(:,iloop,jloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
% TFL2(:,iloop,jloop)=(FNL2(:,iloop,jloop)-FL2(:,iloop,jloop));
% 
% 
% 
% FL3(:,iloop,jloop)=K*(-PHI(:,iloop)*q(iloop)-PHI(:,jloop)*q(jloop));
% mat_name=append("Mode_",num2str(iloop),"_Mode_",num2str(jloop),"_3.rpt");
% FNL3(:,iloop,jloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
% TFL3(:,iloop,jloop)=(FNL3(:,iloop,jloop)-FL3(:,iloop,jloop));


dispA=q(iloop)*PHI(:,iloop)+q(jloop)*PHI(:,jloop);
dispB=-q(iloop)*PHI(:,iloop)-q(jloop)*PHI(:,jloop);
dispC=q(iloop)*PHI(:,iloop)-q(jloop)*PHI(:,jloop);


TFL_1=get_nonlinstiffness(mesh, beam, dispA)*dispA;
TFL_2=get_nonlinstiffness(mesh, beam, dispB)*dispB;
TFL_3=get_nonlinstiffness(mesh, beam, dispC)*dispC;



Gten(:,iloop,jloop)=(TFL1+TFL2-2*Gten(:,iloop,iloop)*q(iloop)^2-2*Gten(:,jloop,jloop)*q(jloop)^2)/(4*q(iloop)*q(jloop));
Gten(:,jloop,iloop)=Gten(:,iloop,jloop);


for vloop=1:Nmodes

G(vloop,iloop,jloop)=PHI(:,vloop)'*Gten(:,iloop,jloop);
G(vloop,jloop,iloop)=G(vloop,iloop,jloop);

H(vloop,iloop,jloop,jloop)=((PHI(:,vloop)'*(TFL1+TFL3))-2*G(vloop,iloop,iloop)*q(iloop)^2-2*G(vloop,jloop,jloop)*q(jloop)^2-2*H(vloop,iloop,iloop,iloop)*q(iloop)^3)/(6*q(iloop)*q(jloop)^2);
H(vloop,jloop,iloop,jloop)=H(vloop,iloop,jloop,jloop);
H(vloop,jloop,jloop,iloop)=H(vloop,iloop,jloop,jloop);

H(vloop,iloop,iloop,jloop)=((PHI(:,vloop)'*(TFL2+TFL3))-2*G(vloop,iloop,iloop)*q(iloop)^2-2*G(vloop,jloop,jloop)*q(jloop)^2+2*H(vloop,jloop,jloop,jloop)*q(jloop)^3)/(-6*q(iloop)^2*q(jloop));
H(vloop,iloop,jloop,iloop)=H(vloop,iloop,iloop,jloop);
H(vloop,jloop,iloop,iloop)=H(vloop,iloop,iloop,jloop);
end


    end
end




for vloop=1:1:N
for iloop=1:1:N
    for jloop=iloop+1:1:N
        for kloop=jloop+1:1:N 

dispA=q(iloop)*PHI(:,iloop)+q(jloop)*PHI(:,jloop)+q(kloop)*PHI(:,kloop);
TFL = get_nonlinstiffness(mesh, beam, dispA)*dispA;


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
   Zs(:,iloop,jloop)=inv((Omega(iloop)+Omega(jloop))^2*M-K)*Gten(:,iloop,jloop);
   Zd(:,iloop,jloop)=inv((-Omega(iloop)+Omega(jloop))^2*M-K)*Gten(:,iloop,jloop);    
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


% [qa,qb]=generate_DNF_inp(PHI_m,q,a_ten,b_ten,dof_per_node,Nnodes,Nmodes,BCdofs,Current_path,maxidx,meshfile);
% 
% 
% system('abaqus cae noGUI=inp_files/DNF_test');




qa=1e-3*ones(Nmodes,1);
qb=1e-3*ones(Nmodes,1);



%% 
Tensortype={'a','b','s','d'};

for aorb=1:2



for iloop=1:Nmodes  %master_modes_selected
    for jloop=iloop:Nmodes  %master_modes_selected

% mat_name=append(Tensortype{aorb},"_ten_",num2str(iloop),num2str(jloop),"_2",".rpt");
% FNL1(:,iloop,jloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
% 
% 
% mat_name=append(Tensortype{aorb},"_ten_",num2str(iloop),num2str(jloop),"_1",".rpt");
% FNL2(:,iloop,jloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);


if aorb==1
dispA=qa(jloop)*qa(kloop)*a_ten(:,jloop,kloop);
dispB=-qa(jloop)*qa(kloop)*a_ten(:,jloop,kloop); 
else
dispA=qb(jloop)*qb(kloop)*b_ten(:,jloop,kloop);
dispB=-qb(jloop)*qb(kloop)*b_ten(:,jloop,kloop); 
end

FNL1=get_nonlinstiffness(mesh, beam, dispA)*dispA;
FNL2=get_nonlinstiffness(mesh, beam, dispB)*dispB;




for vloop=1:Nmodes  %master_modes_selected

if aorb==1
        GA(vloop,iloop,jloop)=PHI(:,vloop)'*(FNL1(:,iloop,jloop)+FNL2(:,iloop,jloop))/(qa(iloop))^2/(qa(jloop))^2/2;
    else
        GB(vloop,iloop,jloop)=PHI(:,vloop)'*(FNL1(:,iloop,jloop)+FNL2(:,iloop,jloop))/(qb(iloop))^2/(qb(jloop))^2/2;
end


end
    end
end



for iloop=1:1:Nmodes
    for jloop=1:1:Nmodes
        for kloop=jloop:1:Nmodes


% mat_name=append("Phi_",num2str(iloop),'_',Tensortype{aorb},"_ten_",num2str(jloop),num2str(kloop),"_1",".rpt");
% FNL1(:,iloop,jloop,kloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
% 
% mat_name=append("Phi_",num2str(iloop),'_',Tensortype{aorb},"_ten_",num2str(jloop),num2str(kloop),"_2",".rpt");

% FNL2(:,iloop,jloop,kloop)=read_rpt(mat_name,BCdofs,dof_per_node,Nnodes);
if aorb==1
dispA=q(iloop)*PHI(:,iloop)+a_ten(:,jloop,kloop)*qa(jloop)*qa(kloop);
dispB=-q(iloop)*PHI(:,iloop)-a_ten(:,jloop,kloop)*qa(jloop)*qa(kloop);
else
dispA=q(iloop)*PHI(:,iloop)+b_ten(:,jloop,kloop)*qb(jloop)*qb(kloop);
dispB=-q(iloop)*PHI(:,iloop)-b_ten(:,jloop,kloop)*qb(jloop)*qb(kloop);
end
FNL1 = get_nonlinstiffness(mesh, beam, dispA)*dispA;
FNL2= get_nonlinstiffness(mesh, beam, dispB)*dispB;


for vloop=1:1:Nmodes
if aorb==1
AA(vloop,iloop,jloop,kloop)=(PHI(:,vloop)'*(FNL1+FNL2)-2*G(vloop,iloop,iloop)*(q(iloop))^2-2*GA(vloop,jloop,kloop)*qa(jloop)^2*qa(kloop)^2)/(4*q(iloop)*qa(jloop)*qa(kloop));
else
BB(vloop,iloop,jloop,kloop)=(PHI(:,vloop)'*(FNL1+FNL2)-2*G(vloop,iloop,iloop)*(q(iloop))^2-2*GB(vloop,jloop,kloop)*qb(jloop)^2*qb(kloop)^2)/(4*q(iloop)*qb(jloop)*qb(kloop));
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