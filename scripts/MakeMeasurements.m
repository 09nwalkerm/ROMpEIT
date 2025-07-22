%% Make EEG-EIT Measurements

model = '/home/c1616132/Documents/PhD/colin_sources/head_model_CEM.mat';

mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];
GM = 5;

sample_num = 11;
dp_num = 1:6;

current_injection = 2e-05;

data = getenv("ROMEG_DATA");

for ii = 1:length(sample_num)

    load([data '/Result' num2str(sample_num(ii)) '/Results/prep_SEP.mat'],'EEG')
    conds = [EEG.conds 5];
    disp(conds)
    GenMeasurements('model',model,'noise',0.82e-6,'mu_min',mu_min,'mu_max',mu_max,...
        'sample_num',sample_num(ii),'current',current_injection,...
        'debug',true,'Cluster',true,'synth_cond',conds)
    
end