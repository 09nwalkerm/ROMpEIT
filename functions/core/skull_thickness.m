function [thickness] = skull_thickness(p,t)%,direct)

    t2 = t(t(:,5)==2,1:4);
    t3 = t(t(:,5)==4,1:4);
%     thickness = zeros(size(t2,1),1);
% 
%     disp('Creating centroids...')
%     centroids = meshcentroid(p,t2);
%     centroids2 = meshcentroid(p,t3);
%     
%     tic
%     disp('Starting pool of workers...')
%     parpool(12)
%     parfor i=1:size(centroids,1)   
%         PQ = dsearchn(centroids2,centroids(i,:));
%         points = centroids2(PQ,:);
%         v = centroids(i,:) - points;
%         thickness(i) = vecnorm(v,2,2);
%     end
%     disp('Finshed thickness calc.')
%     toc
    
    
    %indx = unique(t2(:,1:4));
    indx2 = unique(t3(:,1:4));
    %pp = p(indx,:);
    ppp = p(indx2,:);
    thickness = zeros(length(p(:,1)),1);
    
    tic
    disp('Starting pool of workers...')
    parpool(12)
    parfor i=1:size(p,1)   
        PQ = dsearchn(ppp,p(i,:));
        point = ppp(PQ,:);
        v = p(i,:) - point;
        thickness(i) = vecnorm(v,2,2);
    end
    disp('Finshed thickness calc.')
    toc
end