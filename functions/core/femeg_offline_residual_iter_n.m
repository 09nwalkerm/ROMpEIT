function [ROM]= femeg_offline_residual_iter_n(FOM,ROM,L1,U1,tol)
%OFFLINE_RESIDUAL offline computation of the mu-independent terms of the 
%dual norm of the residual 
%   
%   [CQQ, DQQ, EQQ] = OFFLINE_RESIDUAL(FOM, V) given a FOM struct and a
%   trial basis V, computes the mu-independent terms of the 
%   dual norm of the residual. CQQ is a ROM.Qf x ROM.Qf cell array of numbers,  
%   DQQ is a ROM.Qf x ROM.Qa cell array of ROM.N x 1 vectors,  while EQQ is 
%   a ROM.Qa x ROM.Qa cell array of ROM.N x ROM.N matrices.

nv=size(ROM.V,2);

%fprintf('\n ** Computing dqq and Eqq ** \n')
if norm(FOM.Xnorm-speye(size(FOM.Xnorm,1)),'inf')~=0
    for q1=1:FOM.Qa
        [ROM.Z{nv,q1},flag_z,res_z]=pcg(FOM.Xnorm,FOM.Aq{FOM.active(q1)}*ROM.V(:,nv),tol,2000,L1,U1);
        disp(['  Done Z: ' num2str(q1) ' out of ' num2str(FOM.Qa) '. Flag: ' num2str(flag_z) '. Res: ' num2str(res_z)])
    end

    for q1=1:FOM.P
        for q2=1:FOM.P
            AA=zeros(nv);
            for ii=1:nv
                AA(ii,:)=ROM.Z{ii,q1}'*(FOM.Aq{FOM.active(q2)}*ROM.V);
            end
            ROM.Eqq{q1,q2}=AA;
        end
        for q2=1:FOM.Qf
            BB=zeros(nv,1);
            for ii=1:nv
                BB(ii)=ROM.Z{ii,q1}'*FOM.Fq{q2};
            end
            ROM.dqq{q1,q2}=BB;
        end
    end
else
    if ROM.N == 1
        for q1=1:FOM.P
                ROM.XAV{q1} = sparse(FOM.Aq{q1}*ROM.V);
            for q2=1:FOM.P
                ROM.Eqq{q1,q2}=(ROM.XAV{q1})'*(FOM.Aq{q2}*ROM.V); %tempory edit - added (1:FOM.np) to all terms to select only stiffness matrix
            end
            for q2=1:FOM.Qf
                ROM.dqq{q1,q2}=(ROM.XAV{q1})'*FOM.Fq{q2}; %same as above
            end
        end
    else
        for q1=1:FOM.P
            ROM.XAV{q1} = sparse([ROM.XAV{q1} FOM.Aq{q1}*ROM.V(:,end)]);
        end
        for q1=1:FOM.P
            for q2=1:FOM.P
                SIDE = ROM.XAV{q1}(:,1:end-1)'*ROM.XAV{q2}(:,end);
                BOTTOM = ROM.XAV{q1}(:,end)'*ROM.XAV{q2};
                ROM.Eqq{q1,q2}= [ROM.Eqq{q1,q2} SIDE; BOTTOM];
            end
            for q2=1:FOM.Qf
                BOTTOM = (ROM.XAV{q1}(:,end))'*FOM.Fq{q2};
                ROM.dqq{q1,q2}= [ROM.dqq{q1,q2}; BOTTOM];
            end
        end
    end
end
end


