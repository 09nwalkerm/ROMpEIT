%% RunInverse

% A script to run the inverse problem given an RBModel exists in the path
% specified by the evironment variable $ROMEG_DATA

% Optional: create synthetic measurements before running the inverse
% problem.

%% Prelims

clear

tree = getenv("ROMEG");
%model = [tree '/models/Real/head_model.mat'];
%model2 = '/cubric/data/c1616132/DataDump/ROMEG/Experiments/Real/head_model.mat';
model = '/home/c1616132/Documents/PhD/real_data/Cardiff_bundle/test_models/head_model_ST.mat';
%model1 = '/home/c1616132/Documents/PhD/fat/model/test_mesh_fatless.mat';
%model2 = '/home/c1616132/Documents/PhD/fat/model/test_mesh.mat';
%model3 = '/home/c1616132/Documents/PhD/fat/model/test_mesh2.mat';
%model4 = '/home/c1616132/Documents/PhD/fat/model/test_mesh_marrow.mat';

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

num_samples=1;
num_start=1;

c = [1 2 6];
%synth = [0.4,0.01,0.01,1.6,0.33,0.12,5];

current_injection = 2.3923e-05;

%% Make measurements
% Type `help GenMeasurements` for examples, help and options

% GenMeasurements('model',model4,'noise',0.82e-6,'mu_min',mu_min,'mu_max',mu_max,...
%     'sample_num',num_start:num_samples,'current',current_injection,...
%     'debug',true,'Cluster',true,'pre_stiff',true)

%% Run the inverse problem

%OrderedModelClass.patterns('model',model3,'num_sinks',1,'new_sinks',true,'electrodes',1:163)

% GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%     'active_layers',c,'use_noise',true,...
%     'tag','1-133-known','noise',0.82e-6,'debug',true)

% GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%     'active_layers',c,'use_noise',true,'use_sinks',true,'fix_conds',true,...
%     'tag','','noise',0.82e-6,'debug',true,'ref_sink',257,'snaps',true)

GenInverse('model',model,'ROM',true,'sample_num',num_start:num_samples,...
    'active_layers',c,'new_sinks',true,'fix_conds',[1.79,0.33,0.2,1.5,1.5,3.14],'real',true,...
    'tag','sim3_fixed_test','debug',true,'ref_sink',257,'ground',257,...
    'current',current_injection,'weighted',true,'simultaneous',true)
