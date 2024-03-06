function D = makeCondMatrix(p,t,f,estimates,theta)
    
    % closest centroid for the skull
    % closest csf centroid to that piece of skull
    % closest skull centroid to that csf
    % collect 266 vlaues and inteprolate across the skull
    logger = log4m.getLogger();
    
    if (size(estimates,1) == 1)
        synth_cond = estimates;
        D = zeros(size(t,1),6);
        layers = length(synth_cond);
        for i=1:layers
            D(t(:,5)==i,[1,4,6])=synth_cond(i);
        end
        logger.info('makeCondMatrix','Made conductivity tensors')
        return
    end
        
    t2 = t(t(:,5)==2,1:4);
    t3 = t(t(:,5)==3,1:4);

    logger.info('makeCondMatrix','Creating centroids...')
    centroids2 = meshcentroid(p,t2);
    centroids3 = meshcentroid(p,t3);
    logger.info('makeCondMatrix','Centroids made... finding electrode centres')

    elec_centers = zeros(size(estimates,1),3);

    for ii = 1:size(estimates,1)
        nodes = f(f(:,4)==ii,1:3);
        nodes = unique(nodes);
        elec_centers(ii,:) = [mean(p(nodes,1)) mean(p(nodes,2)) mean(p(nodes,3))];
    end
    logger.info('makeCondMatrix','Electrode centers found...')

    skull1 = zeros(size(estimates,1),3);
    csf1 = zeros(size(estimates,1),3);
    skull2 = zeros(size(estimates,1),3);
    logger.info('makeCondMatrix','Starting loops...')

    for i = 1:size(estimates,1)
        % closest skull centroid to electrode
        %indx = unique(t2(:,1:4));
        %pp = p(indx,:);
        PQ = dsearchn(centroids2,elec_centers(i,:));
        skull1(i,1:3) = centroids2(PQ,1:3);

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
    F1 = scatteredInterpolant([skull1; skull2],[estimates(:,2); estimates(:,2)]);
    VQ = F1(centroids2);
    D1= VQ;

    F2 = scatteredInterpolant([skull1; skull2],[estimates(:,3); estimates(:,3)]);
    VQ = F2(centroids2);
    D2= VQ;

    obj = struct();
    if size(estimates,2)==length(unique(t(:,5)))
        obj.anis_tan = [];
        obj.anis_rad = [];
    else
        obj.anis_tan = 3;
        obj.anis_rad = 2;
        
        oo = zeros(length(t2(:,1)),1);
        logger.info('makeCondMatrix','Constructing CBM matrix for anisotropic layers')
        obj.CBM = [cos(theta(:,1)).*cos(theta(:,2)) sin(theta(:,1)).*cos(theta(:,2)) -sin(theta(:,2)) -sin(theta(:,1)) cos(theta(:,1)) oo cos(theta(:,1)).*sin(theta(:,2)) sin(theta(:,1)).*sin(theta(:,2)) cos(theta(:,2))];

    end


    D = zeros(size(t,1),6);
    layers = length(estimates(1,:));
    for i=1:layers
        is_anit = obj.anis_tan(obj.anis_tan==i);
        is_anir = obj.anis_rad(obj.anis_rad==i);
        tt = i;
        for ii=1:length(obj.anis_rad)
            if obj.anis_rad(ii) < i
                tt = i - 1;
            else
                tt = i;
            end
        end
        if (isempty(is_anit)) && (isempty(is_anir))
            logger.info('makeCondMatrix','Calculating values for isotropic layer')
            if tt==2
                D(t(:,5)==tt,[1,4,6])=[D1 D1 D1];
            else
                D(t(:,5)==tt,[1,4,6])=mean(estimates(:,i));
            end
        elseif isempty(is_anit)
            logger.info('makeCondMatrix','Calculating radial values for stiffness')
            D_tmp = zeros(size(t,1),6);
            D_tmp(t(:,5)==tt,6)=D1;
            D_tmp(t(:,5)==tt,:) = change_basis(D_tmp(t(:,5)==tt,:),obj.CBM,"rad");
            D = D + D_tmp;
        else
            logger.info('makeCondMatrix','Calculating tangential values for stiffness')
            D_tmp = zeros(size(t,1),6);
            D_tmp(t(:,5)==tt,[1,4])=[D2 D2];
            D_tmp(t(:,5)==tt,:) = change_basis(D_tmp(t(:,5)==tt,:),obj.CBM,"tan");
            D = D + D_tmp;
        end
    end
end