%% Generate a RBModel

% Script to determine which layers are best to estimate by the reduction in
% overall estimation error.

%% Prelims
tree = getenv("ROMEG");
%model = [tree '/models/Spherical/head_model.mat']; % path to model to use
%model = '/home/c1616132/Documents/PhD/new_colin/CEM_semi_complete_with_marrow_133_merged.mat';
%model = '/home/c1616132/Documents/PhD/real_data/Cardiff_bundle/test_models/head_model_AM.mat';
model='/home/c1616132/Documents/PhD/colin_sources/head_model_CEM.mat';

%mu_min = [0.1359,0.0008,0.3,0.06,0.065,0.0010,1.5,1.5,3.14];
%mu_max = [0.6196,0.0131,2.0,0.7391,0.228,0.800,1.5,1.5,3.14];

%mu_min = [0.1359,0.0008,0.01,1.3881,0.0598,0.0652,5];
%mu_max = [0.6196,0.0131,0.0430,1.7945,0.7391,0.2283,5];

mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

%mu_min = [0.303,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,1.794,0.508,0.177,5];

current_injection = 2e-05;

%% Set up files

%OMC = OrderedModelClass();
%OMC = OMC.checkPaths('type','ROM');

%% Create sink file
% Type `help OrderedModelClass.patterns` for help, options and examples

%OrderedModelClass.patterns('model',model,'electrodes',1:129,'out',129,'pre_stiff',true)

%% Run functions to make RB model
% Type `help GenRBModel` for help, options and examples

% GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'tolGREEDY',5e-10,'Nmax',40,...
%     'use_sinks',true,'nic',150,'debug',true,...
%     'Cluster',true,'current',current_injection)

GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'tolGREEDY',5e-10,'Nmax',150,...
    'use_sinks',true,'nic',150,'debug',true,...
    'Cluster',true,'current',current_injection,'use_FOM',true)