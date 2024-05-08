function [Ksk,D,d] = makeSkullCond(ptot,ttot,estimates,patterns,elecs)
    % to make a scalable sitff, estimates must be normalised
    
    logger = log4m.getLogger();
    
    t=ttot(ttot(:,5)~=10,:);
    [p,t(:,1:4)]=removeisolatednode(ptot,t(:,1:4)); % Remove isolated nodes
    
    t1 = t(t(:,5)==2,1:4);
    t2 = t(t(:,5)==3,1:4);
    %t8 = t(t(:,5)==8,1:4);

    logger.info('makeCondMatrix','Creating centroids...')
    %centroids8 = meshcentroid(p,t8);
    centroids2 = meshcentroid(p,t1);
    centroids3 = meshcentroid(p,t2);
    logger.info('makeCondMatrix','Centroids made... finding electrode centres')

    elec_centers = ptot(elecs(:,1),:);

    skull1 = zeros(size(estimates,1),3);
    csf1 = zeros(size(estimates,1),3);
    skull2 = zeros(size(estimates,1),3);
    logger.info('makeCondMatrix','Starting loops...')

    for i = 1:size(estimates,1)
        % closest skull centroid to electrode
        %indx = unique(t2(:,1:4));
        %pp = p(indx,:);
        
        PQ = dsearchn(centroids2,elec_centers(patterns(i,1),:));
        skull1(i,1:3) = centroids2(PQ,1:3);
        %skull1(i,1:3) = elec_centers(patterns(i,1),:);%centroids2(PQ,1:3);
        

        % closest csf centroids3 to skull centroid
        PQ = dsearchn(centroids3,skull1(i,:));
        csf1(i,1:3) = centroids3(PQ,1:3);

        % closest skull to csf
        PQ = dsearchn(centroids2,csf1(i,:));
        skull2(i,1:3) = centroids2(PQ,1:3);
        %disp(['Finished electrode ' num2str(i)])
    end

    logger.info('makeCondMatrix','Making inteprolant ...')
    % interp
    
    estimates1=estimates(:,1)/mean(estimates(:,1));
    
    F1 = scatteredInterpolant([skull1; skull2],[estimates1; estimates1],'natural','linear');
    VQ = F1(centroids2);
    D1= VQ;
    D1(isnan(D1))=1;
    D1(D1<0.05)= 0.0001;
    D1(D1>1.8)=0.0001;
    
    D = zeros(size(t,1),6);
    D(t(:,5)==1,[1,4,6])=[D1 D1 D1];
    logger.info('makeCondMatrix','Calculating values for isotropic layer')

    np = size(p,1);
    L=length(elecs);
    Ksk = sparse([femeg_stiffness(p,t,D),zeros(np,L);zeros(L,np+L)]);
    
    d = zeros(size(p,1),1);
    t_p = unique(t(t(:,5)==2,1:4));
    VQ2 = F1(p(t_p,1:3));
    d(t_p) = VQ2;
    
    d(isnan(d))=1;
    d(d<0.05)= 0.0001;
    d(d>1.9)=0.0001;
    logger.info('makeCondMatrix','Done')
end