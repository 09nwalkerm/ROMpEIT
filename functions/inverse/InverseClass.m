classdef InverseClass < OrderedModelClass

    properties
        estimate        % array of estimated conductivities
        estimates       % 2D-array of estimates from all injection patterns
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
        pattern         % electrode pattern number in sink file
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
        ref_sink        % the reference electrode used as a sink for all pairs.
        iter
        real
    end

    methods

        function runInverse(obj,pattern)

            obj.pattern = pattern;

            if ~isempty(obj.new_sinks) && obj.new_sinks && isempty(obj.real)
                if ~isempty(obj.simultaneous) && obj.simultaneous
                    sink_nums=unique(obj.sinks)';sink_nums=sink_nums(1:end);    
                    noise = obj.use_noise;
                    obj.use_noise = [];
                    obj = obj.loadMultiMeasurements(sink_nums);
                    obj.use_noise = noise;
                    measurements=obj.u;obj.u=[];
                    for jj = 1:size(obj.sinks,1)
                        obj = obj.combineMeasurements(jj,measurements);
                    end
                else
                    noise = obj.use_noise;
                    obj.use_noise = [];
                    obj = obj.loadMultiMeasurements(obj.sinks(obj.pattern,:));
                    obj.use_noise = noise;
                    measurements=obj.u;obj.u=[];
                    obj = obj.combineMeasurements(obj.pattern,measurements);
                end
            elseif ~isempty(obj.real) && obj.real && obj.new_sinks %&& ~isempty(obj.simultaneous)
                obj = obj.loadMultiMeasurements(1:size(obj.sinks,1));
            else
                if ~isempty(obj.simultaneous) && obj.simultaneous
                    obj = obj.loadMultiMeasurements(1:size(obj.sinks,1));
                else
                    obj = obj.loadMeasurements(obj.pattern);
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
                    load([path '/pattern_' num2str(num) '_noise.mat'])
                    obj.logger.debug('loadMeasurements',['Loading noisey synthetic pattern ' num2str(num) ' solution...'])
                else
                    load([path '/pattern_' num2str(num) '.mat'])
                    obj.logger.debug('loadMeasurements',['Loading synthetic pattern ' num2str(num) ' solution...'])
                end
                
                if isempty(obj.real) || ~obj.real
                    obj.synth_cond = Data.synth_cond;
                    obj.u{num} = Data.u;
                else
                    obj.u{num} = u;
                end
            catch
                if num==obj.ref_sink
                    obj.logger.debug('loadMeasurements',['Measurements for reference electrode ' num2str(obj.ref_sink) ' do not exist'])
                    return
                else
                    obj.logger.error('loadMeasurements',['Measurements not found for pattern ' num2str(num)])
                    error('Measurement data missing, please either supply real data or generate synthetic data using GenMeasurements.')
                end
            end
        end
        
        function obj = combineMeasurements(obj,num,measurements)
            
            el_out = obj.sinks(num,2:end);
            el_in = obj.sinks(num,1);
            u_comb = measurements{el_in};
            obj.logger.debug('combineMeasurements',['Starting measurement combo with ' num2str(el_in)])
            obj.logger.debug('combineMeasurements',['Starting measurements have electrode values ' num2str(u_comb(end-132:end-122)') '...' num2str(u_comb(end)')])

            for ii=el_out
                if ~(ii == obj.ref_sink)
                    u_comb = u_comb - measurements{ii}/length(el_out);
                    obj.logger.debug('combineMeasurements',['Linearly combining measurement ' num2str(ii)])
                else
                    obj.logger.debug('combineMeasurements',['Skipping ref electrode ' num2str(ii)])
                end
            end
            obj.logger.debug('combineMeasurements',['Combined measurements have values ' num2str(u_comb(end-132:end-122)') '...' num2str(u_comb(end)')])
            
            obj.u{num} = u_comb;
            
            if ~isempty(obj.use_noise) && obj.use_noise
                obj = obj.addNoise(num);
            end
            
            obj.logger.info('combineMeasurements',['Combined measurements for pattern ' num2str(num)])
        end
        
        function obj = loadMultiMeasurements(obj,num)

            for ii=num
                obj = obj.loadMeasurements(ii);
            end
        end
        
        function obj = addNoise(obj,num)

            len = length(obj.u{num}); %un must be column vector
            
            r = normrnd(0,obj.noise,len,1);
            
            M = obj.u{num} + r;

            obj.u{num} = M;
        end
    end
end




































