function cluster_job_BETA()

coef_n=str2double(getenv('SLURM_ARRAY_TASK_ID'));
top = getenv("ROMEG_TOP");
fprintf('\n ************************* COEFFICIENT NUMBER  %d \n', coef_n);

FOM = FOMClass.loadFOM(top,false);

fprintf('Starting beta calculation...\n\n')
betaa=femeg_ROM_RBF_offline_dual_iter(FOM,coef_n);

disp(betaa)

save([top '/Results/ROM/other/betaa_' num2str(coef_n)],'betaa')

end