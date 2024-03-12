function cluster_job_measurement()

    coef_n=str2double(getenv('SLURM_ARRAY_TASK_ID'));
    fprintf('\n ************************* COEFFICIENT NUMBER  %d \n', coef_n);
    top = getenv("ROMEG_TOP");
    setenv("ROMEG_LOGS",[top '/Results/logs'])    
    load([top '/Results/measurements/prep.mat'], 'Data')
    Data.log_tag = ['Data_' num2str(coef_n)];
    Data=Data.startLogger();
    Data = Data.genData(coef_n);
    Data.saveData(coef_n);
end