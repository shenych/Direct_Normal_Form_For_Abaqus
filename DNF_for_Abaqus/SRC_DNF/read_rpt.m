function RFRM=read_rpt(file_name,BCdofs,dof_per_node,Nnodes)


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
                    RFRM(kicount,1)=data_raw(k);
                    kicount=kicount+1;
                end
            end

RFRM(BCdofs)=[];
end