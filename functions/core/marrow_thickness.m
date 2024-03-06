function [ratio] = marrow_thickness(p,t,f)%,direct)
    
    %t2 = t(t(:,5)==2,1:4);
    t3 = t(t(:,5)==4,1:4);
    
    elec_centers = zeros(size(unique(f(:,4)),1)-1,3);
    for ii = 1:length(elec_centers)
        nodes = f(f(:,4)==ii,1:3);
        nodes = unique(nodes);
        elec_centers(ii,:) = [mean(p(nodes,1)) mean(p(nodes,2)) mean(p(nodes,3))];
    end
    
    %indx = unique(t2(:,1:4));
    indx2 = unique(t3(:,1:4));
    %pp = p(indx,:);
    ppp = p(indx2,:);
    ratio = zeros(size(elec_centers,1),1);
    
    tic
    disp('Starting pool of workers...')
    %parpool(12)
    
    for i=1:size(elec_centers,1)
        % Searching for closest CSF tetrahedron index
        PQ = dsearchn(ppp,elec_centers(i,:));
        point = ppp(PQ,:);
        v = elec_centers(i,:) - point;
        % creating points along the vector
        lx = linspace(elec_centers(i,1),point(1),100);
        ly = linspace(elec_centers(i,2),point(2),100);
        lz = linspace(elec_centers(i,3),point(3),100);
        % searching for the tetrehdron around each point
        INDX = zeros(1,100);
        parfor j=1:100
            INDX(j) = tsearchn(p,t(:,1:4),[lx(:,j),ly(:,j),lz(:,j)]);
        end
        if isempty(find(isnan(INDX)))
            layers = t(INDX,5);
            count = accumarray(layers,ones(1,100));
            %disp(count)
            if size(count,1) > 2
                ratio(i,:) = count(3)/(count(2) + count(3));
            end
            disp(['Finished electrode ' num2str(i)])
        else
            disp(['Finished electrode ' num2str(i) '. Set to zero due to NaN.'])
        end
        disp(['The ratio for this electrode is ' num2str(ratio(i,:))])
    end
    disp('Finshed thickness calc.')
    toc
end