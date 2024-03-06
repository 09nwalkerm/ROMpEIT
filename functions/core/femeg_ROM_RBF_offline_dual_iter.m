function [betaa,res]=femeg_ROM_RBF_offline_dual_iter(FOM,n_coef)

% FOM
% mu_train: samples for training
warning('off','all')
%M_mu_mu=femeg_ROM_ensemble_mu_dual(FOM,FOM.mu_train(n_coef,:));
M_mu_mu = FOM.muAssemble(FOM.mu_train(n_coef,:));

opt.UPDATEP = 'no';opt.INNERIT=8;opt.MAXITERATION=100;opt.UPDATEM='no';opt.DISP=0;
[betaa,~,res]= bleigifp(0.5*(M_mu_mu+M_mu_mu'),FOM.Xnorm,2,opt);

betaa=real(betaa(2));res=res(2);

warning('on','all')
