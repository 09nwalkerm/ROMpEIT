function [theta,centroids,centroids2] = angles_for_mesh(p,t,f)%,direct)

tissue = 2; %skull
t2 = t(t(:,5)==tissue);
centroids = meshcentroid(p,t2);
centroids2 = meshcentroid(p,f(:,1:3));
theta = zeros(length(t2(:,1)),2);

%if center == 1

    tic
    for ii = 1:length(t2(:,1))
        v = gpuArray(centroids2) - gpuArray(centroids(ii,:));
        v = gather(v);
        norms = vecnorm(v,2,2);
        [~,index] = min(norms);
        theta2 = acos(v(index,3)/norm(v(index,:)));
        theta1 = acos(v(index,1)/norm(v(index,[1,2])));
        if (v(index,2) < 0)
            theta1 = -1*theta1;
        end
        
        theta(ii,:) = [theta1 theta2];
    end
    toc;
    
    %save([direct '/head_model_theta2.mat'],'p','t','f','theta','centroids','centroids2')
    
% elseif center == 0
%     
%     tic
%     for ii = 1:length(t2(:,1))
%         v = gpuArray(centroids2) - gpuArray(centroids(ii,:));
%         v = gather(v);
%         norms = vecnorm(v,2,2);
%         [~,index] = min(norms);
%         theta2 = acos(v(index,3)/norm(v(index,:)));
%         theta1 = acos(v(index,1)/norm(v(index,[1,2])));
%         if (v(index,2) < 0)
%             theta1 = -1*theta1;
%         end
%         
%         theta(ii,:) = [theta1 theta2];
%     end
%     toc;
%     
%     save([direct '/head_model_theta.mat'],'p','t','f','theta','centroids','centroids2')
% end
