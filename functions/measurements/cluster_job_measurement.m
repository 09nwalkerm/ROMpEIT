function cluster_job_measurement()

    coef_n=str2double(getenv('SLURM_ARRAY_TASK_ID'));
    fprintf('\n ************************* COEFFICIENT NUMBER  %d \n', coef_n);
    top = getenv("ROMEG_TOP");
    load([top '/Results/measurements/prep.mat'], 'Data')
    setenv("ROMEG_LOGS",[top '/Results/logs'])
    Data=Data.startLogger(['Data_' num2str(coef_n)]);
    Data = Data.genData(coef_n);
    Data.saveData(coef_n);
end