function [FOM,RBModel] = GenRBModel(varargin)
%
% GenRBModel(name1, value1, name2, value2 ...)
%
% Description:
%   GenRBModel is the main function in the ROMEG toolbox. It will generate
%   a FOM.mat object and RBModel.mat object in the Results/ folder. These
%   are the full order and reduced order head models respectively. Simply
%   running this function will generate example models but the function
%   comes with many name-value pair options that can be used to tailor the
%   head model. See below for the arguments to this function. After running
%   this function, you are encouraged to run the inverse problem. For
%   guidance with this type: `help GenInverse`.
%
% Requirements:
%   - A small sized cluster (50 cores ~ 1 large cpu node)
%   - MATLAB R2021a or later
%   - Statistics and Machine Learning Toolbox
%
% Arguments:
%   - model      - path to head model you want to use. Must contain p,t,f.
%   - mu_min     - array of minimum conductivities in order of head tissue.
%   - mu_max     - array of maximum conductivities in order of head tissue.
%   - nic        - number of beta interpolation points.
%   - anis_tan   - array of tissue numbers to tag with tangential
%                  conductivity.
%   - anis_rad   - array of tissue numbers to tag with radial
%                  conductivity.
%   - angles     - true if model contains theta values (default: false)
%   - current    - injection current used.
%   - tolGREEDY  - The minumum error tolerance the greedy algorithm
%                  will stop at.
%   - Nmax       - The maximum number of snapshots. If specified
%                  with tolGREEDY, the algorithm will stop at which 
%                  ever value is reached first. Recommended < 300.
%   - verbose    - true or false: if true returns additional values in
%                  RBModel and prints debugging information. This produces
%                  much more data and therefore only runs for one electrode
%                  locally, for 2 snapshots with only one stability factor 
%                  calculation.
%   - use_sinks  - (boolean) would you like to use the sinks.mat file for
%                  injection patterns
%   - complim    - the limit on the number of nodes used at one time on the
%                  cluster
%   - use_FOM    - is there an existing FOM to look for instead of making a
%                  new one.
%
% Examples:
%   
%   This will create an RBModel with the max and min conductivites
%   specified, with appropriate stopping conditions. The nuber of beta
%   interpolation points will be 100 and ROM will be trained according to
%   the injection pattern file provided. There will alse be a limit of 60
%   jobs allowed to run in parallel.
%
%   GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'tolGREEDY',...
%       5e-10,'Nmax',200,'use_sinks',true,'nic',100,'complim',60)
%
%   This will create an RBModel similar to above but where the second layer
%   is isotropic, meaning parameters 2 and 3 are the tangential and radial
%   conductivities for this layer and where the model contains the 
%   necessary angle values for the computation.
%
%   GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'tolGREEDY',...
%       5e-8,'Nmax',100,'use_sinks',true,'nic',80,'anis_tan',3,...
%       'anis_rad',2,'angles',true,'complim',60)
%
%

    %% Preliminaries
    
    paramsFOM = [];
    paramsROM = [];
    ROMlist = [{'electrode'},{'current'},{'tolGREEDY'},{'Nmax'},{'verbose'},{'top'},{'split'},{'use_sinks'}];
    FOMlist = [{'model'},{'mu_min'},{'mu_max'},{'nic'},{'anis_tan'},{'anis_rad'},{'angles'},{'verbose'},{'top'},{'complim'},{'use_FOM'}];

    if ~isempty(varargin)
        for i = 1:2:length(varargin) % work for a list of name-value pairs
            if ischar(varargin{i}) && ~isempty(find(strcmp(ROMlist,varargin{i})))
                paramsROM = [paramsROM varargin(i) varargin(i+1)];
            end
            if ischar(varargin{i}) && ~isempty(find(strcmp(FOMlist,varargin{i})))
                paramsFOM = [paramsFOM varargin(i) varargin(i+1)];
            end
        end
    end

    obj = OrderedModelClass();
    obj = obj.checkPaths('type','ROM');
    
    % Get the top of the ROMEG tree
    top = getenv("ROMEG_TOP");
    addpath(genpath([top '/Functions']));
    paramsROM = [paramsROM {'top'} {top}];
    paramsFOM = [paramsFOM {'top'} {top}];

    % path to head model
    if isempty(find(strcmp('model',paramsFOM)))
        disp("Using example Spherical head model.")
        tree = getenv("ROMEG");
        paramsFOM=[{'model'} {strcat(tree,'/Models/Spherical/head_model.mat')} paramsFOM];
    end
    
    if isempty(find(strcmp('angles',paramsFOM)))
        paramsFOM = [paramsFOM {'angles'} {false}];
    end
    
    %% Step 1: Compute FOM

%     if isfile([top '/Results/ROM/FOM.mat'])
%         prompt = "A FOM.mat file already exists. \nWould you like to use this instead of calculating a new one? (y/n): ";
%         txt = input(prompt,"s");
%         if isempty(txt)
%             txt = 'y';
%         end
%         if txt == "y" || txt == "Y"
    if ~isempty(find(strcmp('use_FOM',paramsFOM)))
        
            disp('Using FOM.mat from Results/ folder')
            load([top '/Results/ROM/FOM.mat'],'FOM')
            if ~isempty(find(strcmp('use_sinks',paramsROM)))
                OMC = OrderedModelClass();
                OMC = OMC.checkPaths('type','ROM');
                OMC = OMC.loadSinks();
                setenv('SE',num2str(OMC.num_patterns));
            else
                setenv('SE',num2str(FOM.L-1));
                OMC.num_patterns = FOM.L-1;
            end
    else
            FOM = FOMClass(paramsFOM);
            FOM = FOM.buildFOM();
            FOM = FOM.saveROMparams(paramsROM);
            FOM.saveFOM(top);
            setenv('SL',num2str(FOM.nic));
            if ~isempty(find(strcmp('use_sinks',paramsROM)))
                OMC = OrderedModelClass();
                OMC = OMC.checkPaths('type','ROM');
                OMC = OMC.loadSinks();
                setenv('SE',num2str(OMC.num_patterns));
            else
                setenv('SE',num2str(FOM.L-1));
                OMC.num_patterns = FOM.L-1;
            end
            %% Step 2 
            % Compute interpolation coefficients (CLUSTER)
            % and collect back into FOM Class Object
            
            % Run job in parallel
%             if find(strcmp('verbose',paramsFOM))
%                 FOMClass.verboseBeta(top,FOM);
%             else
                if isempty(find(strcmp('complim',paramsFOM)))
                    !sbatch --array 1-$SL -o $ROMEG_TOP/Results/slurm_logs/BETA_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/BETA_%a_%j.err --job-name BETA $ROMEG/Functions/Cluster/cluster_job.sh BETA
                else
                    setenv("COMPLIM",num2str(FOM.complim))
                    !sbatch --array 1-$SL%$COMPLIM -o $ROMEG_TOP/Results/slurm_logs/BETA_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/BETA_%a_%j.err --job-name BETA $ROMEG/Functions/Cluster/cluster_job.sh BETA
                end
                
                % Wait until the job is done
                OrderedModelClass.wait('BETA',40);
                
                FOM=FOM.readBeta(top,FOM.nic); 
                FOM=femeg_ROM_RRBF(FOM); % Compute RRBF function
                FOM.saveFOM(top);
%             end
    end

    %% Step 3
    % Compute ROM for each electrode (CLUSTER)
    % and collect back into ROM Class Object
%     if find(strcmp('verbose',paramsROM))
%         ROM=ROMClass('electrode',1,'top',top,'verbose',true);
%         ROM = ROM.buildROM();
%         ROM.saveLF(1);
%         RBModel = RBModelClass('verbose',true);
%         RBModel = RBModel.readLF(top,1);
%         RBModel.saveROM(top);
%     else
        % Run job in parallel
        if isempty(find(strcmp('complim',paramsFOM)))
            !sbatch --array 1-$SE -o $ROMEG_TOP/Results/slurm_logs/ROM_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/ROM_%a_%j.err --job-name ROM $ROMEG/Functions/Cluster/cluster_job.sh ROM
        else
            setenv("COMPLIM",num2str(FOM.complim))
            !sbatch --array 1-$SE%$COMPLIM -o $ROMEG_TOP/Results/slurm_logs/ROM_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/ROM_%a_%j.err --job-name ROM $ROMEG/Functions/Cluster/cluster_job.sh ROM
        end
        
        % Wait until the job is done
        OrderedModelClass.wait('ROM',180);
        
        RBModel = RBModelClass();
        RBModel = RBModel.readLF(top,OMC.num_patterns);
        RBModel.saveROM(top);
%     end
end