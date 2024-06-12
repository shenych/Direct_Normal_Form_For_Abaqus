clear
clc

%% Thess are all the inputs

% Path information

path='C:/YichangShen/ROM_challenge_paper/DNF_ABAQUS_Code/Program_v4/';
%current path of the main code and the folder of functions.

meshfile='C:/YichangShen/ROM_challenge_paper/DNF_ABAQUS_Code/Program_v4/Examples/Beam_clamped_1_1.inp'; 
%full path of your mesh file, the mesh should be in the *.inp format that
%ABAQUS can read.


% Modal basis information
master_modes=[1,2];  
% which master modes you want to include in your reduced modal basis?
% It can be as many as you want, using more master modes leads to more
% accurate ROM, while the computational cost will increase accordingly.

%Selecting master modes depends on which modes in your dynamical analysis play the most
%improtant role.


% How large extent of the geogemtric nonlinearity that you want?
thickness=0.03; 
%Thickness of your structure

disp_thickness_ratio=0.3; 
% Maximum described displacement applied to structures for generating the ROMs,  
% 0.2 time of thickness as default and is suggested for some reasons.

% Element information
dof_of_disp=3; % In the type of element of the mesh, how many dofs for the displacements for each node? (DX DY DZ)
dof_of_rotation=0; % and how many dofs for the rotations?  (DRX DRY DRZ)


%% The code will launch
disp_applied=disp_thickness_ratio*thickness;


warning('off', 'MATLAB:rmpath:DirNotFound');rmpath(genpath(path));
addpath(genpath(append(path,'SRC_DNF')))  %path
[AH,BH,G,H,a_ten,b_ten,r_ten,Omega,PHI]=DNF_in_FE(master_modes,dof_of_disp,dof_of_rotation,disp_applied,path,meshfile);


%% this is to save the mat file including all coefficents computed by the DNF
matfile=strrep(meshfile,'.inp','.mat');
save(matfile)

