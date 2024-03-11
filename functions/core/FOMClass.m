classdef FOMClass < OrderedModelClass
    properties
        mu_min
        mu_max
        active
        non_active
        Qa
        Qf
        P
        Pp
        L
        np
        nt
        Aq
        ANq
        Fq
        FNq
        nic
        Xnorm
        mu_train
        betaa
        RRBF
        CBM
        paramsROM
        complim
    end
    methods

        function obj = FOMClass(varargin)
        %
        % FOMClass(p,t,f,name1, value1, name2, value2 ...)
        %
        % Description:
        %   A class object for the construction of the full order model.
        %   Needed for the construction of the Reduced order model in made
        %   by using ROMClass. Please use the GenRBModel function to
        %   generate this class object although feel free to use this for
        %   debugging purposes. See below for arguments.
        %
        % Args:
        %   model: path to head model you want to use. Must contain p,t,f.
        %   mu_min: array of minimum conductivities in order of head tissue.
        %   mu_max: array of maximum conductivities in order of head tissue.
        %   nic: number of beta interpolation points.
        %   anis_tan: array of tissue numbers to tag with tangential
        %                conductivity.
        %   anis_rad: array of tissue numbers to tag with radial
        %                conductivity.
        %   current: injection current used.
        %
            
            if isa(varargin{1},'cell')
                varargin = varargin{1};
            end

%             for i=1:2:length(varargin)
%                 if ischar(varargin{i}) % check if is character
%                     obj.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
%                 end
%             end

            obj = obj.processArgs(varargin);
            obj = obj.getTOP();

        end

        function obj = buildFOM(obj)
        % FOM(p,t,f,mu_min,value,mu_max,value,anis_tan,value,anis_rad,value,nic,value)
            obj = obj.processModel();
        
            if ~isempty(obj.anis_rad)
                if isempty(obj.theta) || ~obj.angles
                    obj = obj.computeAngles(obj.t);
                else
                    obj = obj.computeAngles(obj.t,obj.theta);
                end
            end
        
            obj = obj.assembleFOM();

            obj = obj.computeStiff();
            
            obj = obj.computeXnorm();

            obj = obj.computeMuTrain();

            disp('Done')
        end
        function obj = assembleFOM(obj)

            p = obj.p; t = obj.t; f= obj.f;

            obj.np = length(p);
            obj.nt = length(t);
            obj.L=length(unique(f(:,end)))-1;
            obj.Pp=length(unique(t(:,5)));

            if isempty(obj.nic)
                obj.nic = 20;
                fprintf("Using default number of 20 for Stability Factor Interpolation points. \n")
            end

            if isempty(obj.mu_min) && isempty(obj.mu_max)
                % Use defult values
                if obj.Pp == 4
                    fprintf("Using default conducitivity parameters. \n" + ...
                        "Please specify your own in future using the name-value pair mu_min and mu_max. \n")
                    obj.mu_max=[0.65,0.04,1.71,0.33,5];
                    obj.mu_min=[0.15,0.001,1.71,0.33,5];
                elseif obj.Pp == 6
                    % Range for parameters: [scalp,skull-tan,skull-rad,csf,wm,gm,Z]
                    fprintf("Using default conducitivity parameters. \n" + ...
                        "Please specify your own in future using the name-value pair mu_min and mu_max. \n")
                    obj.mu_max=[0.66,0.06,0.06,1.71,0.22,0.47,5];
                    obj.mu_min=[0.15,0.001,0.001,1.71,0.22,0.47,5];
                end
            end

            layers = obj.mu_max - obj.mu_min;
            obj.active = find(logical(layers));
            obj.non_active = find(~logical(layers));

            obj.Qa=length(obj.active);%Pp+npz; % tissue layers + 1 for Z
            obj.Qf=1;
            obj.P=length(obj.mu_max);%Pp+npz; % total number of parameters

%             if length(obj.mu_min) == obj.active(end)
%                 % Z
%                 [A,B,C]=femeg_stiffness_cem(p,t,f,1);
%                 obj.Aq{length(obj.mu_min)}=[A,-B;-B',C];
%                 clear A B C Dn
%                 active = obj.active(1:end-1);
%                 obj = obj.computeStiff(active,obj.non_active,p,t);
%             else
%                 % Z
%                 [A,B,C]=femeg_stiffness_cem(p,t,f,1);%mu_min(end));
%                 obj.Aq{length(obj.mu_min)}=[A,-B;-B',C];
%                 clear A B C Dn
%                 non_active = obj.non_active(1:end-1);
%                 obj = obj.computeStiff(obj.active,non_active,p,t);
%             end

        end

        function obj = computeAngles(obj,t,varargin)
            if nargin < 3
                obj = obj.computeTheta();
            end
            tissue = min([obj.anis_tan,obj.anis_rad]); % only support for one anisotropic layer at a time
%             tissue = [];
%             for ii = 1:length(obj.anis_tan)
%                 tissue_tmp = obj.anis_tan(ii)-(ii-1);
%                 tissue = [tissue tissue_tmp];
%             end
            t2 = t(t(:,5)==tissue);

            if ~isempty(obj.anis_tan)
                theta = obj.theta;
                % Construct change of basis matrix, 2=theta, 1=psi
                oo = zeros(length(t2(:,1)),1);
                obj.CBM = [cos(theta(:,1)).*cos(theta(:,2)) sin(theta(:,1)).*cos(theta(:,2)) -sin(theta(:,2)) -sin(theta(:,1)) cos(theta(:,1)) oo cos(theta(:,1)).*sin(theta(:,2)) sin(theta(:,1)).*sin(theta(:,2)) cos(theta(:,2))];
            end
        end

        function obj = computeTheta(obj)
            disp(['Computing theta angle for each tetrahedra. /n ' ...
                'This could take a few minutes.'])
            tree = getenv("ROMEG");
            obj.theta = angles_for_mesh(obj.p,obj.t,obj.f,[tree '/Models']);
        end

        function obj = computeStiff(obj)
            
            p = obj.p; t = obj.t; f = obj.f;

            [A,B,C]=femeg_stiffness_cem(p,t,f,1);
            obj.Aq{length(obj.mu_min)}=[A,-B;-B',C];
            clear A B C Dn

            if length(obj.mu_min) == obj.active(end)
                active = obj.active(1:end-1);
                non_active = obj.non_active; % removed (1:end-1);
            else
                active = obj.active;
                non_active = obj.non_active(1:end-1);
            end

            %%%%% Stiffness matrices
            disp('Computing stiffness matrices...')
            for kk=active%1:Pp
                if isempty(obj.anis_tan)
                    Dn=zeros(obj.nt,6);Dn(t(:,5)==kk,[1,4,6])=1;
                    obj.Aq{kk}=sparse([femeg_stiffness(p,t,Dn),zeros(obj.np,obj.L);zeros(obj.L,obj.np+obj.L)]);
                    disp(['Done ' num2str(kk) ' from ' num2str(obj.Pp)])
                else
                    is_anit = obj.anis_tan(obj.anis_tan==kk);
                    is_anir = obj.anis_rad(obj.anis_rad==kk);
                    tt = kk;
                    lower = min([obj.anis_tan,obj.anis_rad]);
                    for ii=1:length(lower)
                        if lower(ii) < kk
                            tt = kk - 1;
                        else
                            tt = kk;
                        end
                    end
            
                    if (isempty(is_anit)) && (isempty(is_anir))
                        Dn=zeros(obj.nt,6);Dn(t(:,5)==tt,[1,4,6])=1;
                        obj.Aq{kk}=sparse([femeg_stiffness(p,t,Dn),zeros(obj.np,obj.L);zeros(obj.L,obj.np+obj.L)]);
                        disp(['Done ' num2str(kk) ' from ' num2str(length(obj.mu_min)) ' in tissue ' num2str(tt)])
                    elseif isempty(is_anit)
                        Dn=zeros(obj.nt,6);Dn(t(:,5)==tt,6)=1;
                        Dn(t(:,5)==tt,:) = change_basis(Dn(t(:,5)==tt,:),obj.CBM,"rad");
                        obj.Aq{kk}=sparse([femeg_stiffness(p,t,Dn),zeros(obj.np,obj.L);zeros(obj.L,obj.np+obj.L)]);
                        disp(['Done ' num2str(kk) ' from ' num2str(length(obj.mu_min)) ' in tissue ' num2str(tt)])
                    else
                        Dn=zeros(obj.nt,6);Dn(t(:,5)==tt,[1,4])=1;
                        Dn(t(:,5)==tt,:) = change_basis(Dn(t(:,5)==tt,:),obj.CBM,"tan");
                        obj.Aq{kk}=sparse([femeg_stiffness(p,t,Dn),zeros(obj.np,obj.L);zeros(obj.L,obj.np+obj.L)]);
                        disp(['Done ' num2str(kk) ' from ' num2str(length(obj.mu_min)) ' in tissue ' num2str(tt)])
                    end
                end
            end
            
            for kk=non_active%1:Pp
                if isempty(obj.anis_tan)
                    Dn=zeros(obj.nt,6);Dn(t(:,5)==kk,[1,4,6])=1;
                    obj.Aq{kk}=sparse([femeg_stiffness(p,t,Dn),zeros(obj.np,obj.L);zeros(obj.L,obj.np+obj.L)]);
                    disp(['Done ' num2str(kk) ' from ' num2str(obj.Pp)])
                else
                    is_anit = obj.anis_tan(obj.anis_tan==kk);
                    is_anir = obj.anis_rad(obj.anis_rad==kk);
                    tt = kk;
                    lower = min([obj.anis_tan,obj.anis_rad]);
                    for ii=1:length(lower)
                        if obj.lower(ii) < kk
                            tt = kk - 1;
                        else
                            tt = kk;
                        end
                    end
                        
                    if (isempty(is_anit)) && (isempty(is_anir))
                        Dn=zeros(obj.nt,6);Dn(t(:,5)==tt,[1,4,6])=1;
                        obj.Aq{kk}=sparse([femeg_stiffness(p,t,Dn),zeros(obj.np,obj.L);zeros(obj.L,obj.np+obj.L)]);
                        disp(['Done ' num2str(kk) ' from ' num2str(length(obj.mu_min)) ' in tissue ' num2str(tt)])
                    elseif isempty(is_anit)
                        Dn=zeros(obj.nt,6);Dn(t(:,5)==tt,6)=1;
                        Dn(t(:,5)==tt,:) = change_basis(Dn(t(:,5)==tt,:),obj.CBM,"rad");
                        obj.Aq{kk}=sparse([femeg_stiffness(p,t,Dn),zeros(obj.np,obj.L);zeros(obj.L,obj.np+obj.L)]);
                        disp(['Done ' num2str(kk) ' from ' num2str(length(obj.mu_min)) ' in tissue ' num2str(tt)])
                    else
                        Dn=zeros(obj.nt,6);Dn(t(:,5)==tt,[1,4])=1;
                        Dn(t(:,5)==tt,:) = change_basis(Dn(t(:,5)==tt,:),obj.CBM,"tan");
                        obj.Aq{kk}=sparse([femeg_stiffness(p,t,Dn),zeros(obj.np,obj.L);zeros(obj.L,obj.np+obj.L)]);
                        disp(['Done ' num2str(kk) ' from ' num2str(length(obj.mu_min)) ' in tissue ' num2str(tt)])
                    end
                end
            end
        end

        function obj = computeXnorm(obj)
            disp('Computing norm...')
            obj.Xnorm= speye(obj.np+obj.L);%blkdiag(femeg_ROM_Xnorm_EEG(p,t),speye(L)); % eq 3.25 from Somersalo, with c=0 / speye(np+L);%
        end

        function obj = computeMuTrain(obj)
            % Compute a positive RBF interpolant of the stability factor
            disp('Generating random RBF conductivity interpolation set...')
            mu_cube=lhsdesign(obj.nic,length(obj.active)); % normalized design
            obj.mu_train=bsxfun(@plus,obj.mu_min(obj.active),bsxfun(@times,mu_cube,(obj.mu_max(obj.active)-obj.mu_min(obj.active))));
        end

        function obj = readBeta(obj)
            betaa_tot=zeros(obj.nic,1);
            for kk=1:obj.nic
                load([obj.top '/Results/ROM/other/betaa_' num2str(kk)], 'betaa')
                betaa_tot(kk)=betaa;
            end
            obj.betaa=betaa_tot;
        end

        function saveFOM(obj)

            FOM = obj;
            save([obj.top '/Results/ROM/FOM.mat'], 'FOM')
            fprintf("\n \n To see FOM.mat, please check the Results folder in the ROMEG tree. \n \n")

        end

        function obj = loadpreFOM(top)
            load([top '/Results/ROM/FOM.mat'], 'FOM')
            obj = FOM;
            disp("Loaded pre-made FOM Class from Results/.")
        end

        function obj = saveROMparams(obj,paramsROM)
            obj.paramsROM = paramsROM;
        end
    end

    methods (Static)
        function obj = loadFOM(top)
        % 
        % loadFOM(top)
        % 
        % Description:
        %   Load a previously constructed Full Order Model object from the
        %   Results/ folder.
        %
        % Arguments:
        %   - top - path to top of ROMEG tree


            load([top '/Results/ROM/FOM.mat'],'FOM')
            obj = FOM;
            disp("Loaded pre-made FOM Class object from Results/.")
        end

    end
end
