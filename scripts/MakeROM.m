%% Generate a RBModel

% Script to determine which layers are best to estimate by the reduction in
% overall estimation error.

%% Prelims
tree = getenv("ROMEG");
model = [tree '/models/Spherical/head_model.mat']; % path to model to use
%model = '/home/c1616132/Documents/PhD/real_data/Cardiff_bundle/test_models/head_model_108I.mat';
%model = '/home/c1616132/Documents/PhD/fat/model/test_mesh.mat';

%mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
%mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

%scalp,comp_skull,CSF,GM,WM,marrow,eyes,gel,z
%mu_min = [0.303,0.002,1.450,0.268,0.092,0.013,1.5,1.5,5];
%mu_max = [0.444,0.009,1.794,0.508,0.177,0.043,1.5,1.5,5];
%mu_min = [0.1359,0.0008,1.0,0.0598,0.0652,0.0010,1.5,1.5,3.14];
%mu_max = [0.6196,0.0131,1.7945,0.7391,0.2283,0.043,1.5,1.5,3.14];

%mu_min = [0.303,0.02,0.002,0.268,0.092,1.450,5];
%mu_max = [0.444,0.06,0.02,0.508,0.177,1.794,5];

%mu_min = [0.1359,0.002,0.268,0.092,1.450,5];
%mu_max = [0.6196,0.02,0.508,0.177,1.794,5];

%mu_min = [0.1359,0.02,0.002,0.268,0.092,1.450,5];
%mu_max = [0.6196,0.2,0.02,0.508,0.177,1.794,5];

mu_min = [0.303,0.002,1.450,0.1,5];
mu_max = [0.444,0.043,1.794,0.508,5];

current_injection = 2.0e-05;

%% Set up files

OMC = OrderedModelClass();
OMC = OMC.checkPaths('type','ROM');

%% Create sink file
% Type `help OrderedModelClass.patterns` for help, options and examples

OrderedModelClass.patterns('model',model,'electrodes',1:120,'out',120)

%% Run functions to make RB model
% Type `help GenRBModel` for help, options and examples

GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'tolGREEDY',5e-10,'Nmax',40,...
    'use_sinks',true,'nic',150,'debug',true,...
    'Cluster',true,'current',current_injection)