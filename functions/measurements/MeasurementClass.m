classdef MeasurementClass < OrderedModelClass

    properties
        SL
        synth_cond
        injection
        el_in
        u
        current
        num_samples
        mu_min
        mu_max
        noise
        ROM
        LF
        CBM
        ratio
        redo
    end

    methods
        function obj = MeasurementClass(varargin)
            
            obj@OrderedModelClass(varargin);
        end

        function obj = genData(obj,injection)
            
            obj.el_in = obj.sinks(injection,1);

            p = obj.p; t = obj.t; f= obj.f;
            obj.logger.debug('genData',['size of t ' num2str(size(t))])
            
            if ~isempty(obj.anis_tan)
                theta = obj.theta;
                % Construct change of basis matrix, 2=theta, 1=psi
                tissue = min([obj.anis_tan,obj.anis_rad]); % only support for one anisotropic layer at a time
                obj.logger.info('genData',['tissue for anisotropy ' num2str(tissue)])
%                 tissue = [];
%                 for ii = 1:length(obj.anis_tan)
%                     tissue_tmp = obj.anis_tan(ii)-(ii-1);
%                     tissue = [tissue tissue_tmp];
%                 end
                t2 = t(t(:,5)==tissue,:);
                oo = zeros(length(t2(:,1)),1);
                obj.logger.info('genData','Constructing CBM matrix for anisotropic layers')
                obj.CBM = [cos(theta(:,1)).*cos(theta(:,2)) sin(theta(:,1)).*cos(theta(:,2)) -sin(theta(:,2)) -sin(theta(:,1)) cos(theta(:,1)) oo cos(theta(:,1)).*sin(theta(:,2)) sin(theta(:,1)).*sin(theta(:,2)) cos(theta(:,2))];
            end
            
            layers = length(obj.synth_cond)-1; % -1 to account for the impedance (Z) number.

            D = zeros(size(t,1),6);
            for i=1:layers
                if isempty(obj.anis_tan)
                    if isempty(t(t(:,5)==i))
                        obj.logger.error('genData','Number of conductivities given and layers in model does not match. Please adjust.')
                        error('Number of conductivities given and layers in model does not match. Please adjust.')
                    end
                    Ind = (t(:,5)==i);D(Ind,[1,4,6])=obj.synth_cond(i);
                    obj.logger.info('genData',['Calculating stiffness matrix for layer ' num2str(i)])
                else
                    is_anit = obj.anis_tan(obj.anis_tan==i);
                    is_anir = obj.anis_rad(obj.anis_rad==i);
                    tt = i;
                    lower = min([obj.anis_tan,obj.anis_rad]);
                    for ii=1:length(lower)
                        if lower(ii) < i
                            tt = i - 1;
                        else
                            tt = i;
                        end
                    end
                    if (isempty(is_anit)) && (isempty(is_anir))
                        obj.logger.info('genData','Calculating values for isotropic layer')
                        D(t(:,5)==tt,[1,4,6])=obj.synth_cond(i);
                    elseif isempty(is_anit)
                        obj.logger.info('genData','Calculating radial values for stiffness')
                        D_tmp = zeros(size(t,1),6);
                        D_tmp(t(:,5)==tt,6)=obj.synth_cond(i);
                        D_tmp(t(:,5)==tt,:) = change_basis(D_tmp(t(:,5)==tt,:),obj.CBM,"rad");
                        D = D + D_tmp;
                    else
                        obj.logger.info('genData','Calculating tangential values for stiffness')
                        D_tmp = zeros(size(t,1),6);
                        D_tmp(t(:,5)==tt,[1,4])=obj.synth_cond(i);
                        D_tmp(t(:,5)==tt,:) = change_basis(D_tmp(t(:,5)==tt,:),obj.CBM,"tan");
                        D = D + D_tmp;
                    end
                end
            end

            %%
            np=size(p,1);
            eL = length(unique(f(:,end)))-1;
            
            %S=EIT_stiffness(p,t,D);
            S=sparse([femeg_stiffness(p,t,D),zeros(np,eL);zeros(eL,np+eL)]);
            [A,B,C] = femeg_stiffness_cem(p,t,f,obj.synth_cond(end)); %Adjusted from 5 to more general obj.synth_cond(end)
            S_cem = [A, -B; -B', C];
            S_cem = S_cem + S;
            [L,U] = ilu(S_cem);
            
            % Injection nodes
            el_in=obj.el_in; I_in=obj.current;%p_in=unique(f(f(:,4)==el_in));I_in=1/size(p_in,1);
            el_out=obj.sinks(injection,2:end); I_out=(-1*I_in)/length(el_out);%p_out=unique(f(f(:,4)==el_out));I_out=-1/size(p_out,1);
            %b = sparse([p_in;p_out],1,[I_in;I_out],size(p,1),1);
            b = spalloc(size(p,1)+eL,1,length(el_out)+1);
            b(np+el_in) = I_in;
            b(el_out+np) = I_out;
            
            [u,flag] = pcg(S_cem,b,1e-10,4000,L,U);
            
            %u2 = u(end-(eL-1):end);
            
            obj.u = u;

        end
        
        function obj = genDatafromROM(obj,injection)
            
            obj.injection = injection; 
            
            if obj.use_sinks
                obj.el_in = injection;
            else
                obj.el_in = obj.sinks(obj.injection,1);
            end

            el_in = obj.el_in; el_out = obj.sinks(obj.injection,2:end);
            
            mu_a = obj.synth_cond;
            
            n_mu=length(mu_a); % number of parameters
            M_mu=obj.LF{el_in}.ANq{n_mu}/mu_a(n_mu); % Z
            
            for kk=1:n_mu-1
                M_mu=M_mu+mu_a(kk)*obj.LF{el_in}.ANq{kk}; % Stiffness
            end
            
            zN=M_mu\obj.LF{el_in}.FNq{1};
            zNh = obj.LF{el_in}.V*zN;
            zNh1 = zNh;
            
            if ~obj.use_sinks
                for jj=1:length(el_out)
                    M_mu=obj.LF{el_out(jj)}.ANq{n_mu}/mu_a(n_mu); % Z

                    for kk=1:n_mu-1
                        M_mu=M_mu+mu_a(kk)*obj.LF{el_out(jj)}.ANq{kk}; % Stiffness
                    end

                    zN=M_mu\obj.LF{el_out(jj)}.FNq{1};
                    zNh = obj.LF{el_out(jj)}.V*zN;
                    zNh2{jj} = zNh;
                end

                zNh3 = zeros(size(zNh1));
                for kk=1:length(el_out)
                    zNh = (zNh1-zNh2{kk})/length(el_out);
                    zNh3 = zNh3 + zNh;
                end
            else
                zNh3 = zNh1;
            end
        
            obj.u = zNh3;
        end
        
        function obj = loadLF(obj)
            obj.logger.info('loadLF','Loading RBModel from Results folder')
            load([obj.top '/Results/ROM/RBModel.mat'],'RBModel')
            if ~RBModel.LF{1}.use_sinks
                obj.LF = {};
                N_list = [];
                for i=[obj.el_in obj.sinks(obj.injection,2:end)]
                    obj.LF{i} = RBModel.LF{i};
                    N_list = [N_list RBModel.LF{i}.N];
                end
            else
                obj.LF = RBModel.LF;
                obj.use_sinks = RBModel.LF{1}.use_sinks;
            end
        end
        
        function savePrep(obj)
            Data = obj;
            save([obj.top '/Results/measurements/prep.mat'], 'Data')
            obj.logger.debug('savePrep',['Prep.mat saved: ' obj.top '/Results/measurements/prep.mat'])
        end

        function saveData(obj,num)
%             Data = obj;
%             
%             % Clean up
%             Data.p = [];
%             Data.t = [];
%             Data.f = [];
%             Data.Ind_E = [];
            Data = struct();
            Data.u = obj.u;
            Data.synth_cond = obj.synth_cond;
            
            if isempty(obj.noise)
                save([obj.top '/Results/measurements/pattern_' num2str(num) '.mat'], 'Data')
            else
                save([obj.top '/Results/measurements/pattern_' num2str(num) '.mat'], 'Data')
                obj = obj.addNoise();
                Data.u = obj.u;
                save([obj.top '/Results/measurements/pattern_' num2str(num) '_noise.mat'], 'Data')
            end
                
        end

        function obj = addNoise(obj)

            len = length(obj.u); %un must be column vector
            
            r = normrnd(0,obj.noise,len,1);
            
            M = obj.u + r;

            obj.u = M;
        end

    end
end