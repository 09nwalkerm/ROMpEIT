function cluster_job_inv_TRAD()

    coef_n=str2double(getenv('SLURM_ARRAY_TASK_ID'));
    fprintf('\n ************************* COEFFICIENT NUMBER  %d \n', coef_n);
    top = getenv("ROMEG_TOP");
    setenv("ROMEG_LOGS",[top '/Results/logs'])
    load([top '/Results/inverse/TRAD/prep.mat'], 'inv')
    inv.log_tag = ['invTRAD_' num2str(coef_n)];
    inv = inv.startLogger();
    inv.runInverse(coef_n);
end