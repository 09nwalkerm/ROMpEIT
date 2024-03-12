classdef InverseClass < OrderedModelClass

    properties
        estimate        % array of estimated conductivities   
        te              % index for conductivities to estimate
        u               % measurements (artificial measurements)
        lf              % index for conductivities not to estimate
        cond_lf         % conductivities not to estimate
        eL              % number of electrodes
        verbose         % boolean - verbose mode for debugging
        x0              % starting point of the optimisation
        data_path       % path to measurment data
        synth_cond      % conductivity set to make the synthetic measurements with
        el_in           % injection electrode
        injection       % inection electrode pattern number in sink file
        ROM             % is this being used for ROM inverse
        TRAD            % is the traditional inverse method being run as well?
        lb              % lower conductivity bound
        ub              % upper conductivity bound
        current         % injection current
        num_samples     % number of samples to run the inverse problem for
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        nonlcon = [];
        fix_conds       % fix the non active conductivity layers to centre
        active_layers   % when selected by the user, can determine the layers actually used for inverse
        sensitivity
        complim
        use_noise       % use the measurements with noise
        noise           % for adding noise to the combined measurements

    end

    methods

        function runInverse(obj,injection)

            obj.injection = injection;
            
            if ~isempty(obj.use_sinks) && obj.use_sinks
                if ~isempty(obj.simultaneous) && obj.simultaneous
                    obj = obj.loadMultiMeasurements(unique(obj.sinks(:,1))');
                else
                    obj = obj.loadMeasurements(obj.injection);
                end
            elseif ~isempty(obj.new_sinks) && obj.new_sinks
                if ~isempty(obj.simultaneous) && obj.simultaneous
                    sink_nums=unique(obj.sinks)';sink_nums=sink_nums(1:end);
                    obj = obj.loadMultiMeasurements(sink_nums);
                    measurements=obj.u;obj.u=[];
                    for jj = 1:size(obj.sinks,1)
                        obj = obj.combineMeasurements(jj,measurements);
                    end
                else
                    obj = obj.loadMultiMeasurements(obj.sinks(obj.injection,:));
                    measurements=obj.u;obj.u=[];
                    obj = obj.combineMeasurements(obj.injection,measurements);
                end
            end
            
            obj = obj.setUp();

            obj = obj.opt();
            
            obj.saveInv();
        end

        function obj = loadMeasurements(obj,num)
            
            if ~isempty(obj.data_path)
                path = obj.data_path;
                disp('Loading measurement data from path given')
            else
                path = [obj.top '/Results/measurements'];
            end
            
            try
                if ~isempty(obj.use_noise) && obj.use_noise
                    load([path '/pattern_' num2str(num) '_noise.mat'],'Data')
                    disp(['Loading noisey synthetic pattern ' num2str(num) ' solution...'])
                else
                    load([path '/pattern_' num2str(num) '.mat'],'Data')
                    disp(['Loading synthetic pattern ' num2str(num) ' solution...'])
                end
                obj.u{num} = Data.u;
                obj.synth_cond = Data.synth_cond;
            catch
                error('Measurement data missing, please either supply real data or generate synthetic data using GenMeasurements.')
            end
        end
        
        function obj = combineMeasurements(obj,num,measurements)
            
            el_out = obj.sinks(num,2:end);
            el_in = obj.sinks(num,1);
            u_comb = measurements{el_in};

            for ii=el_out
                u_comb = u_comb - measurements{ii}/length(el_out);
            end
            
            obj.u{num} = u_comb;
            
            if ~isempty(obj.use_noise) && obj.use_noise
                obj = obj.addNoise(num);
            end
            
            disp(['Combined measurements for pattern ' num2str(num)])
        end
        
        function obj = loadMultiMeasurements(obj,num)
            
            noise = obj.use_noise;
            obj.use_noise = [];
            for ii=num
                obj = obj.loadMeasurements(ii);
            end
            obj.use_noise = noise;
        end
        
        function obj = addNoise(obj,num)

            len = length(obj.u{num}); %un must be column vector
            
            r = normrnd(0,obj.noise,len,1);
            
            M = obj.u{num} + r;

            obj.u{num} = M;
        end
    end
end




































