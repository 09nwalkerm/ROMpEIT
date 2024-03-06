function [ deltaN ] = femeg_ROM_error_estimate_eit(ROM, uN, mu, FOM, varargin)
%ERROR_ESTIMATE compute an estimate of the RB error for a given parameter
%value
%  
%   [ DELTAN ] = ERROR_ESTIMATE(ROM, UN, MU) given a ROM struct and the RB 
%   solution UN corresponding to the parameter value MU, computes an
%   estimate of the error between the reduced and full order solutions.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

N = length(uN);
%ROM.Eqq{1,1}

active = ROM.active;
non_active = ROM.non_active;

n_mu = zeros(1,ROM.P);
n_mu(active) = mu;
n_mu(non_active) = ROM.mu_min(non_active);

[ theta_a, theta_f ] = femeg_ROM_Theta_EIT( n_mu, ROM );

INDX = 1:N;


%% evaluate a lower bound or an approximation to the stability factor
if isfield(ROM, 'stabFactor')
    beta   = RBF_OnlineInterpolation(ROM, mu);
else
    beta   = 1;
end

%% evaluate dual norm of the residual
res_aa = 0;
res_af = 0;
res_ff = 0;

for q1 = 1 : ROM.Qf
    for q2 = 1 : ROM.Qf
        res_ff = res_ff + theta_f(q1)' * theta_f(q2) * ROM.Cqq{q1,q2}; %same as below
    end
end


for q1 = 1 : ROM.P
    for q2 = 1 : ROM.P
        res_aa = res_aa + theta_a(q1)' * theta_a(q2) * uN'*(ROM.Eqq{q1,q2}(INDX,INDX)*uN); %removed (INDX,INDX) and added (1:FOM.np) to solution terms
    end

    for q2 = 1 : ROM.Qf
        res_af = res_af + theta_a(q1)' * theta_f(q2)  * uN'*ROM.dqq{q1,q2}; %... same as above
                        %+ theta_a(q1)  * theta_f(q2)' * ((ROM.dqq{q1,q2}(INDX))'*uN);
                        %+ theta_a(q1) * theta_f(q2)' * (ROM.dqq{2,q2,q1}(INDX)*uN);

    end
end


res = res_ff - 2*(res_af) + res_aa;

deltaN = sqrt(res);

end
