function [Cqq, dqq, Eqq] = femeg_offline_residual_iter(FOM, V)
%OFFLINE_RESIDUAL offline computation of the mu-independent terms of the 
%dual norm of the residual 
%   
%   [CQQ, DQQ, EQQ] = OFFLINE_RESIDUAL(FOM, V) given a FOM struct and a
%   trial basis V, computes the mu-independent terms of the 
%   dual norm of the residual. CQQ is a ROM.Qf x ROM.Qf cell array of numbers,  
%   DQQ is a ROM.Qf x ROM.Qa cell array of ROM.N x 1 vectors,  while EQQ is 
%   a ROM.Qa x ROM.Qa cell array of ROM.N x ROM.N matrices.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

[L1,U1]=ilu(FOM.Xnorm);
tol=5e-5;

fprintf('\n     Compute Fq''*(X^(-1) Fq) and Fq''*(X^(-1) Aq V) terms   \n')
for q1 = 1 : FOM.Qf
       
    [t,flag_t,res_t]=pcg(FOM.Xnorm,FOM.Fq{q1},tol,2000,L1,U1);
%     disp(['          Done t: ' num2str(q1) ' out of ' num2str(FOM.Qf) '. Flag: ' num2str(flag_t) '. Res: ' num2str(res_t)])
    
    for q2 = 1: FOM.Qf
        Cqq{q1,q2} = t'*FOM.Fq{q2};
    end

end

fprintf('\n     Compute V''Aq''*(X^(-1) Fq) and V''Aq''*(X^(-1) (Aq * V)) terms \n')
for q1 = 1 : FOM.Qa
    
    Z=zeros(size(FOM.Xnorm,1),size(V,2));
    for kk=1:size(V,2)
        [Z(:,kk),flag_z,res_z]=pcg(FOM.Xnorm,FOM.Aq{q1}*V(:,kk),tol,2000,L1,U1);
%         disp(['            Done Z: ' num2str(kk) ' out of ' num2str(size(V,2)) '. Flag: ' num2str(flag_z) '. Res: ' num2str(res_z)])
    end
    
    disp(['          Done Eqq and dqq: ' num2str(q1) ' out of ' num2str(FOM.Qa)])

    for q2 = 1 : FOM.Qa
        Eqq{q1,q2} = Z'*(FOM.Aq{q2}*V);
    end

    for q2 = 1 : FOM.Qf
        dqq{q1,q2} = Z'*FOM.Fq{q2};
    end
end



end


