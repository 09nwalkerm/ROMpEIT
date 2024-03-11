function cluster_job_inv_TRAD()

    coef_n=str2double(getenv('SLURM_ARRAY_TASK_ID'));
    fprintf('\n ************************* COEFFICIENT NUMBER  %d \n', coef_n);
    top = getenv("ROMEG_TOP");
    load([top '/Results/inverse/TRAD/prep.mat'], 'inv')
    setenv("ROMEG_LOGS",[top '/Results/logs'])
    inv.startLogger(['invTRAD_' num2str(coef_n)])
    inv.runInverse(coef_n);
end