function cluster_job_ROM()

electrode=str2double(getenv('SLURM_ARRAY_TASK_ID'));
top = getenv("ROMEG_TOP");
fprintf('\n ************************* ELECTRODE NUMBER  %d \n \n \n', electrode);

ROM=ROMClass('electrode',electrode,'top',top);
ROM=ROM.buildROM();
ROM.saveLF(electrode)

end