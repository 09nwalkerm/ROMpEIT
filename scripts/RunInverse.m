%% RunInverse

% A script to run the inverse problem given an RBModel exists in the path
% specified by the evironment variable $ROMEG_DATA

% Optional: create synthetic measurements before running the inverse
% problem.

%% Prelims

%clear

tree = getenv("ROMEG");
%model = [tree '/models/Real/head_model.mat'];
%model2 = '/cubric/data/c1616132/DataDump/ROMEG/Experiments/Real/head_model.mat';
model = '/home/c1616132/Documents/PhD/real_data/Cardiff_bundle/test_models/head_model_107.mat';
%model1 = '/home/c1616132/Documents/PhD/fat/model/test_mesh_fatless.mat';
%model2 = '/home/c1616132/Documents/PhD/fat/model/test_mesh.mat';
%model3 = '/home/c1616132/Documents/PhD/fat/model/test_mesh2.mat';
%model4 = '/home/c1616132/Documents/PhD/fat/model/test_mesh_marrow.mat';sq

mu_min = [0.30,0.02,0.002,0.013,0.268,0.092,1.49,5];
mu_max = [0.40,0.2,0.01,0.033,0.508,0.177,1.794,5];

%mu_min = [0.303,0.002,0.268,0.092,1.450,5];
%mu_max = [0.303,0.02,0.508,0.177,1.794,5];

%mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
%mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

%mu_min = [0.303,0.002,0.013,0.0173,0.0079,1.450,0.092,0.268,5];
%mu_max = [0.444,0.009,0.043,0.0173,0.0079,1.794,0.177,0.508,5];

% Anisotropic values in the skull
%mu_min = [0.303,0.002,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,0.043,1.794,0.508,0.177,5];

% For merged skull and spongiform bone
%mu_min = [0.303,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,1.794,0.508,0.177,5];

num_samples=2;
num_start=2;

c = [1 2 3 4 5 6];
%synth = [0.4,0.01,0.01,1.6,0.33,0.12,5];

current_injection = 2.1e-05;

%% Make measurements
% Type `help GenMeasurements` for examples, help and options

% GenMeasurements('model',model4,'noise',0.82e-6,'mu_min',mu_min,'mu_max',mu_max,...
%     'sample_num',num_start:num_samples,'current',current_injection,...
%     'debug',true,'Cluster',true,'pre_stiff',true)

%% Run the inverse problem

% noise_levels=1e-6:2e-6:1.5e-5;
% 
% for i=1:length(noise_levels)
%     GenInverse('model',model,'ROM',true,'sample_num',num_start:num_samples,...
%         'active_layers',c,'use_noise',true,...%,'new_sinks',true,...
%         'tag',['noise' num2str(i)],'ground','avg','ref_sink',133,...
%         'current',current_injection,'fix_conds',true,'snaps',true,...
%         'Cluster',true,'noise',noise_levels(i))
% end

% GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%     'active_layers',c,'use_noise',true,'use_sinks',true,'fix_conds',true,...
%     'tag','','noise',0.82e-6,'debug',true,'ref_sink',257,'snaps',true)

% GenInverse('model',model,'ROM',true,'sample_num',num_start:num_samples,...
%     'active_layers',[1 2 6],'new_sinks',true,'fix_conds',[1.79,0.33,0.2,1.5,1.5,3.14],'real',true,...
%     'tag','marrow_fixed_test2','debug',true,'ref_sink',129,'ground',129,...
%     'current',current_injection,'weighted',true,...
%     'lb',[0.137 0.0008 0.001],'ub',[0.6196 0.013 0.8])%,'simultaneous',true)

% GenInverse('model',model,'ROM',true,'sample_num',num_start:num_samples,...
%     'active_layers',[1 2 3 4 5 6],'new_sinks',true,'fix_conds',true,'real',true,...
%     'tag','marrow_test4','debug',true,'ref_sink',129,'ground',129,...
%     'current',current_injection,'weighted',true,...
%     'lb',[0.137 0.0008 1.38 0.06 0.065 0.001],...
%     'ub',[0.6196,0.0130,1.8,0.7391,0.228,0.30])%,'simultaneous',true)

GenInverse('model',model,'ROM',true,'sample_num',num_start:num_samples,...
    'active_layers',[1 2 3 4 5 6],'new_sinks',true,'fix_conds',true,'real',true,...
    'tag','ref_test6','debug',true,'ref_sink',129,'ground','avg_ref',...
    'current',current_injection,'weighted',true,...
    'lb',[0.137 0.0008 1.38 0.06 0.065 0.001],...
    'ub',[0.6196,0.0130,1.8,0.7391,0.228,0.30])













