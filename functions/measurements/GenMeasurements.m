function GenMeasurements(varargin)
%
%   GenMeasurements(name1,value1,name2,value2...)
%
% Description:
%   Function to generate synthetic measurements from a head model and some
%   provided conductivities (see arguments below). The default current is
%   0.02mA but can be specified with an argument. Please ensure an
%   injection pattern data file is available for the function to run
%   properly. For help with this please see 'help
%   OrderedModelClass.patterns'.
%
% Arguments:
%   synth_cond: (essential) the conductivities of the tissues in
%                         the head model
%   model: path to the head model
%   sinks_path: path the file with the injection patterns
%   sinks: the injection pattern data in a 2D array
%   current: injection current
%   top: path to the top of the ROMEG data tree
%   num_samples: number of conductivity samples
%   sample_num: array of sample numbers to generate (e.g 10:20)
%   mu_min: (essential for num_sample) array of the minimum conductivity values
%   mu_max: (essential for num_sample) array of the maximum conductivity values
%   noise: how much noise do you want to include?
%   ROM: would you like these to be calculated using ROM?
%   use_sinks: tell this function that ROM used 'use_sinks'
%   anis_tan: array of tissue numbers to tag with tangential
%                         conductivity
%   anis_rad: array of tissue numbers to tag with radial
%                         conductivity.
%   angles: true if model contains theta values (default: false)
%   ratio: The ratio between tan and rad conds, e.g, 5
%                         means tan is 5 times rad. A value of 1 will give
%                         an isotropic conductivity.
%   new_sinks: (boolean) use a new set of sinks called new_sinks.mat
%   redo: Redo the measurements using existing prep.mat
%                         file
%   RB_path: path for the RBModel if NOT in $ROMEG_DATA
%   debug: (boolean) turn debug mode on
%   Cluster: should these measurements be made on the cluster?
%
%
% Examples:
%   Basic isotropic measurements with no noise:
%
%   GenMeasurements('model',model,'num_samples',num_samples,'mu_min',mu_min,...
%   'mu_max',mu_max,'noise',0)
%
%   Skull anisotropic measurements with noise and a ratio of tangential to 
%   radial conductivities where the angles are already in the head_model.mat
%   file:
%
%   GenMeasurements('model',model,'num_samples',num_samples,'mu_min',mu_min,...
%   'mu_max',mu_max,'anis_tan',2,'anis_rad',3,'angles',true,'ratio',3,'noise',0.82e-6)
    
    params = [];
    params_S = struct();
    paramslist = [{'synth_cond'},{'model'},{'sinks_path'},{'sinks'},{'current'}, ...
        {'top'},{'num_samples'},{'mu_min'},{'mu_max'},{'noise'},{'ROM'},{'Cluster'}, ...
        {'use_sinks'},{'anis_tan'},{'anis_rad'},{'angles'},{'ratio'},{'new_sinks'}, ...
        {'redo'},{'sample_num'},{'RB_path'},{'debug'}];

    if ~isempty(varargin)
        for i = 1:2:length(varargin) % work for a list of name-value pairs
            if ischar(varargin{i}) && ~isempty(find(strcmp(paramslist,varargin{i})))
                params = [params varargin(i) varargin(i+1)];
                params_S.(varargin{i}) = varargin{i+1};
            end
        end
    end

    if ~isfield(params_S,'model')
        error('Please provide path to head model.')
    end

    if ~isfield(params_S,'current')
        params = [params {'current'} {0.020e-3}];
    end

    if isfield(params_S,'synth_cond') %%~isfield(params_S,'num_samples') && ~isfield(params_S,'sample_num')
        %params_S.num_samples = 1;
        synth_cond = params_S.synth_cond;
    else
        if isfield(params_S,'num_samples')
            params_S.samp_len = params_S.num_samples;
        elseif isfield(params_S,'sample_num')
            params_S.samp_len = length(params_S.sample_num);
        end
        if ~isfield(params_S,'ratio')
            mu_cube=lhsdesign(params_S.samp_len,length(params_S.mu_min)); % normalized design
            synth_cond=bsxfun(@plus,params_S.mu_min,bsxfun(@times,mu_cube,(params_S.mu_max-params_S.mu_min)));
        else
            params_S.mu_min(params_S.anis_tan) = params_S.mu_min(params_S.anis_rad)*params_S.ratio;
            mu_cube=lhsdesign(params_S.samp_len,length(params_S.mu_min)); % normalized design
            synth_cond=bsxfun(@plus,params_S.mu_min,bsxfun(@times,mu_cube,(params_S.mu_max-params_S.mu_min)));
            synth_cond(:,params_S.anis_rad) = synth_cond(:,params_S.anis_tan)/params_S.ratio;
        end
    end
    
    if ~isfield(params_S,'num_samples') && ~isfield(params_S,'sample_num')
        samples = 1;
    elseif isfield(params_S,'sample_num')
        samples = params_S.sample_num;
    elseif isfield(params_S,'num_samples')
        samples = 1:params_S.num_samples;
    end

    Data = MeasurementClass(params);
    Data = Data.processModel();
    %Data.SL = length(unique(Data.f(:,end)))-1;
    
    for i = 1:length(samples)
        
        Data.synth_cond = synth_cond(i,:);
        
        if ~isfield(params_S,'RB_path')
            Data = Data.checkPaths('type','measurement','num',samples(i));
        else
            Data = Data.checkPaths('type','measurement','num',samples(i),'RB_path',params_S.RB_path);
        end
        
        Data = Data.loadSinks();
        
        if ~isfield(params_S,'ROM')
            
            if isfield(params_S,'redo')
                Data_tmp = Data;
                load([Data.top '/Results/measurements/prep.mat'],'Data')
                Data.sinks = Data_tmp.sinks;
                Data.num_sinks = Data_tmp.num_sinks;
                Data.num_patterns = Data_tmp.num_patterns;
                Data.top = Data_tmp.top;
                clear Data_tmp
                disp('Redo Warning: Redo==true: Using previously made prep.mat file to make measurements.')
                disp(['Redo Warning: Head model path: ' Data.model])
                disp(['Redo Warning: top path: ' Data.top])
                disp(['Redo Warning: Conductivities: ' num2str(Data.synth_cond)])
                disp(['Redo Warning: switched sinks to local ROM version of size ' num2str(size(Data.sinks))])
            end
            
            if isfield(params_S,'Cluster') && params_S.Cluster
                Data.savePrep();

                setenv("num",num2str(Data.num_patterns));

                [status,cmdout] = system('sbatch --array 1-$num -o $ROMEG_TOP/Results/slurm_logs/SYNTH_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/SYNTH_%a_%j.err --job-name SYNTH $ROMEG/functions/cluster/cluster_job.sh measurement');
                disp(cmdout)

                OrderedModelClass.wait('SYNTH',20);
            else
                for ii = 1:Data.num_patterns
                    Data = Data.genData(ii);
                    Data.saveData(ii);
                end
            end
        else
            
            Data = Data.loadLF();
            
            for ii=1:size(Data.sinks,1)
                Data = Data.genDatafromROM(ii);
                Data.saveData(ii);
            end

        end
    
        disp('Finished making synthetic measurements. Saved to Results/measurements.')
    end
end