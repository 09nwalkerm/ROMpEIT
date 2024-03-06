classdef RBModelClass < OrderedModelClass
    properties
        LF
        verbose
        use_sinks
    end
    methods
        function obj = RBModelClass(varargin)
            
            obj = obj.processArgs(varargin);
        end

        function obj = readLF(obj,top,n)

            for kk=1:n
                try
                    load([top '/Results/ROM/other/LF_EIT_' num2str(kk) '.mat'], 'ROM')
                    if size(ROM.V,1) > ROM.L
                        ROM.V = ROM.V(end-(ROM.L-1):end,:);
                    end
                catch
                    disp(['Could not load LF_EIT_' num2str(kk) '.mat'])
                    continue
                end
                obj.LF{kk}=ROM;
                clear ROM
                disp(['LF_EIT_' num2str(kk) ' loaded of ' num2str(n) '.'])
            end
            
        end

        function saveROM(obj,top)
        %
        % saveROM(top)
        %
        % Description:
        %   A function to save the reduce order model as RBModel.mat in the
        %   Results/ folder.
        %
        % Arguments:
        %   - top - Path to the top of the ROMEG tree
        %

            RBModel = obj;
            if obj.verbose
                save([top '/Results/verbose/RBModel.mat'],'RBModel')
            else
                save([top '/Results/ROM/RBModel.mat'],'RBModel')
                fprintf("\n \n The reduced order model has been saved in the Results/ folder as RBModel. \n \n")
                fprintf(" \n To run the inverse problem for this model: type `help GenInverse` into the command line. \n \n")
            end
        end

        function runInverse(obj,varargin)
            GenInverse('ROM',true)
        end
        
    end
end