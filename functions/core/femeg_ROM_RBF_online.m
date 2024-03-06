function [betaa]=femeg_ROM_RBF_online(mu_test,FOM)

% hacer eficiente!!!


fi=@(x)exp(-x.^2) ;%x.^2 .* log(x);
%fi=@(x)x.^2.*log(x);

nP=size(FOM.mu_train,1);
nT=size(FOM.mu_train,2);
%mu_test=mu_test(1:FOM.P);

for kk=1:nP
    fi_tot(kk)=fi(norm(mu_test-FOM.mu_train(kk,:)));
end

betaa=exp(FOM.RRBF(end)+mu_test*FOM.RRBF(end-nT:end-1)+fi_tot*FOM.RRBF(1:end-nT-1));