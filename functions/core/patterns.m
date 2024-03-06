function [sinks] = patterns(p,f,num_sinks,varargin)

electrodes = 1:length(unique(f(:,4)))-1;

elec_centers = zeros(size(length(unique(f(:,4))),1),3);
for ii = 1:length(electrodes)
    nodes = f(f(:,4)==electrodes(ii),1:3);
    nodes = unique(nodes);
    elec_centers(ii,:) = [mean(p(nodes,1)) mean(p(nodes,2)) mean(p(nodes,3))];
end

sinks = zeros(size(elec_centers,1),num_sinks);

for ii = 1:length(unique(f(:,4)))-2
    pos = elec_centers(ii,:);

    elec_sinks = elec_centers;
    lengths = zeros(size(elec_sinks,1),3);

    for jj = 1:size(elec_sinks,1)-1
        lengths(jj,:) = elec_sinks(jj,:) - pos(:,:);
        if nargin > 3
            if (elec_sinks(jj,3) < varargin{1}), lengths(jj,:) = []; end
        end
    end

    norms = vecnorm(lengths');

    [a,ind] = maxk(norms,num_sinks);

    sinks(ii,:) = ind;
end
    
% elseif isequal(model,'Spherical')
%     elec_centers = zeros(size(length(unique(f(:,4))),1),3);
%     for ii = electrodes
%         nodes = f(f(:,4)==electrodes(ii),1:3);
%         nodes = unique(nodes);
%         elec_centers(ii,:) = [mean(p(nodes,1)) mean(p(nodes,2)) mean(p(nodes,3))];
%     end
% 
%     sinks = zeros(size(elec_centers,1),num_sinks);
% 
%     for ii = 1:length(unique(f(:,4)))-2
%         pos = elec_centers(ii,:);
% 
%         elec_sinks = elec_centers;
%         lengths = zeros(size(elec_sinks,1),3);
% 
%         for jj = 1:size(elec_sinks,1)-1
%             lengths(jj,:) = elec_sinks(jj,:) - pos(:,:);
%         end
% 
%         norms = vecnorm(lengths');
% 
%         [a,ind] = maxk(norms,num_sinks);
% 
%         sinks(ii,:) = ind;
%     end