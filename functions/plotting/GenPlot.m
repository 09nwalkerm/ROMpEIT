function GenPlot(varargin)
%
%   GenPlot(name1,value1,name2,value2,...)
%
% Description:
%   A function to load and plot the processed data from the specified
%   ROMEG_DATA repository of data.
%
% Arguments:
%   sense           - plot the sensitivity analysis
%   layer_names     - list of names in a cell
%   plot            - type of plot, "Bar", "Box"
%   num_samples     - number of samples to read from
%   anis            - isotropic layer to be compared to
%                     anisotropic counter parts
%   snap            - plot the snap shot comparison data
%   bound           - plot the bound data
%
% Example:
%   After running AllPlot.m script you may wish to make all the plots like
%   so:
%
%   GenPlot('sense',true,'snap',true,'bound',true)
%
%   Or to simply see the sense plot as a Bar graph instead of Box (default):
%
%   GenPlot('sense',true,'plot','Bar')
%

    senseparams = [];
    params_S = [];
    paramslist = [{'sense'},{'snap'},{'bound'}];
    senselist = [{'layer_names'},{'plot'},{'num_samples'},{'anis'}];

    if ~isempty(varargin)
        for i = 1:2:length(varargin) % work for a list of name-value pairs
            if ischar(varargin{i}) && ~isempty(find(strcmp(senselist,varargin{i})))
                senseparams = [senseparams varargin(i) varargin(i+1)];
            end
            if ischar(varargin{i}) && ~isempty(find(strcmp(paramslist,varargin{i})))
                params_S.(varargin{i}) = varargin{i+1};
            end
        end
    end

    OrderedModelClass.changePath('ROM'); top = getenv("ROMEG_TOP");
    
    if isfield(params_S,'sense') && (params_S.sense==true)
        load([top '/Results/sense_data.mat'],'sense')
        sense.plotSensitivity(senseparams)
    end
    
    if isfield(params_S,'snap') && (params_S.snap==true)
        load([top '/Results/snap_data.mat'],'snap')
        snap.plotSnapshots()
    end
    
    if isfield(params_S,'bound') && (params_S.bound==true)
        load([top '/Results/bound_data.mat'],'bound')
        bound.plotBound()
    end

end