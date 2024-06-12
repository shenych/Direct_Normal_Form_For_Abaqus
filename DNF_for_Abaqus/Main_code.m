clear
clc

%% Thess are all the inputs

% Path information

path='C:/your_path';
%current path of the main code and the folder of functions.

meshfile='C:/your_path/*.inp'; 
%full path of your mesh file, the mesh should be in the *.inp format that
%ABAQUS can read.


% Modal basis information
master_modes=[];  
% which master modes you want to include in your reduced modal basis?
% It can be as many as you want, using more master modes leads to more
% accurate ROM, while the computational cost will increase accordingly.

%Selecting master modes depends on which modes in your dynamical analysis play the most
%improtant role.


% How large extent of the geogemtric nonlinearity that you want?
thickness=; 
%Thickness of your structure

disp_thickness_ratio=1; 
% Maximum displacement your ROM will be tested,  one time of thickness as default.

% Element information
dof_of_disp=; % In the type of element of the mesh, how many dofs for the displacements for each node? (DX DY DZ)
dof_of_rotation=; % and how many dofs for the rotations?  (DRX DRY DRZ)


%% The code will launch
disp_applied=disp_thickness_ratio*thickness;


warning('off', 'MATLAB:rmpath:DirNotFound');rmpath(genpath(path));
addpath(genpath(append(path,'SRC_DNF')))  %path
[AH,BH,G,H,a_ten,b_ten,r_ten,Omega,PHI]=DNF_in_FE(master_modes,dof_of_disp,dof_of_rotation,disp_applied,path,meshfile);


%% this is to save the mat file including all coefficents computed by the DNF
matfile=strrep(meshfile,'.inp','.mat');
save(matfile)

