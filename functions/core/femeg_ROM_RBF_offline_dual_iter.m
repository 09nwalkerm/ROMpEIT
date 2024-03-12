function [betaa,res]=femeg_ROM_RBF_offline_dual_iter(FOM,n_coef)

% FOM
% mu_train: samples for training
warning('off','all')
%M_mu_mu=femeg_ROM_ensemble_mu_dual(FOM,FOM.mu_train(n_coef,:));
M_mu_mu = FOM.muAssemble(FOM.mu_train(n_coef,:));

opt.UPDATEP = 'no';opt.INNERIT=8;opt.MAXITERATION=100;opt.UPDATEM='no';opt.DISP=0;
try
    %FOM.logger.debug('femeg_ROM_RBF_offline_dual_iter','Using bleigifp function to calculate beta')
    [betaa,~,res]= bleigifp(0.5*(M_mu_mu+M_mu_mu'),FOM.Xnorm,2,opt);
catch ME
    if (strcmp(ME.identifier,'MATLAB:UndefinedFunction'))
        FOM.logger.warn('femeg_ROM_RBF_offline_dual_iter',...
            'Could not find BLEIGIFP function (download recommended) so using eigs instead.')
        betaa = eigs(M_mu_mu,2,'smallestabs');
        res=0;
    else
        rethrow(ME)
    end
end

betaa=real(betaa(2));%res=res(2);
FOM.logger.info('femeg_ROM_RBF_offline_dual_iter',['Beta value is ' num2str(betaa)])
warning('on','all')
