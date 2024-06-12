function [M,K,BCdofs,Nnodes]=Compute_MK(dof_per_node,meshfile)


file_name=append('Compute_MK.inp');

copyfile(meshfile,file_name)

fileID = fopen(file_name,'a');

text='';

text=append(text,'*STEP, name=exportmatrix',newline,...
'*MATRIX GENERATE, MASS,STIFFNESS',newline,...
'*MATRIX OUTPUT, MASS,STIFFNESS,FORMAT=COORDINATE',newline,...
'*END STEP');




fprintf(fileID,text);
fclose(fileID);


system('abaqus -job Compute_MK');


while isempty(folder_search( pwd, 'mtx'))~=0
pause(5)
end
% delete Compute_MK.com
% delete Compute_MK.sta
% delete Compute_MK.msg
% delete Compute_MK.dat
% delete Compute_MK.prt
% delete Compute_MK.odb
% delete Compute_MK.log
% delete Compute_MK_X1.sim
ListPath = folder_search( pwd, 'mtx');


MTEMP=import_matrix(ListPath{1});
% M=full(MTEMP);
M=MTEMP;

KTEMP=import_matrix(ListPath{2});
% K=full(KTEMP);
K=KTEMP;



% MTEMP=import_stiffness_matrix6(ListPath{1});
% M=full(MTEMP);
% 
% 
% KTEMP=import_stiffness_matrix6(ListPath{2});
% K=full(KTEMP);







BCdofs=find(sum(K>10^35,1));

K(:,BCdofs)=[];
K(BCdofs,:)=[];
M(:,BCdofs)=[];
M(BCdofs,:)=[];

% [PHI,BB]=eig(K,M);
% % [PHI2,BB2]=eig(inv(M)*K);
% for i=1:length(BB(1,:))
%     Omega(i)=sqrt(BB(i,i));
% end



Ndof=length(MTEMP(1,:)); %number of nodes x number of dofs for each node: 1481*6
Nnodes=Ndof/dof_per_node;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end




function ListPath = folder_search( path_target, extension_target )
list_dir = dir( path_target );
ListPath = {  };
for i = 1:1:length( list_dir )
    if isfolder( fullfile( path_target, list_dir( i ).name ) ) && ~strcmp( list_dir( i ).name, '.' ) && ~strcmp( list_dir( i ).name, '..' )
        ListPathSub = folder_search( fullfile(path_target, list_dir( i ).name), extension_target );
        if ~isempty( ListPathSub )
            ListPath = [ ListPath;ListPathSub ];
        end
    else
        idx_dot = regexp( list_dir( i ).name, '\.' );
        file_extension = list_dir( i ).name( ( idx_dot( end  ) + 1 ):end  );
        if ~isempty( find( strcmp( file_extension, extension_target ), 1 ) )
            ListPath = [ ListPath;fullfile( path_target, list_dir( i ).name ) ];
        end
    end
end
end



function [matrix] = import_matrix(mtx_file)
%============== Import Stiffness Matrix ==============%
data_raw = dlmread(mtx_file);

matrix=zeros(max(data_raw(:,1)),max(data_raw(:,2)));
for i=1:length(data_raw(:,1))

    matrix(data_raw(i,1),data_raw(i,2))=data_raw(i,3);

end

end





