function cluster_job_BETA()

coef_n=str2double(getenv('SLURM_ARRAY_TASK_ID'));
top = getenv("ROMEG_TOP");
fprintf('\n ************************* COEFFICIENT NUMBER  %d \n', coef_n);

FOM = FOMClass.loadFOM(top,false);
setenv("ROMEG_LOGS",[top '/Results/logs'])
FOM=FOM.startLogger(['BETA_' num2str(coef_n)]);

betaa=femeg_ROM_RBF_offline_dual_iter(FOM,coef_n);

save([top '/Results/ROM/other/betaa_' num2str(coef_n)],'betaa')
FOM.logger.info('cluster_job_BETA',['Saved beta value to ' top '/Results/ROM/other/betaa_' num2str(coef_n)])
FOM.logger.info('cluster_job_BETA','Closing MATLAB instance.')
end