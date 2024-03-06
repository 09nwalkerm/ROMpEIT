classdef PlottingClass < OrderedModelClass
    properties
        results_ROM
        results_TRAD = {};
        conds
        num_samples
        TRAD_RE
        ROM_RE
        tri
        axis
        cut
        fill
        expr
        recursion
        VQ
        d
        cond_map
        folder
        elec_err
        tissue
        amp_dense
        LIC
        map_tissue = 2
        thickness
        cmap = 'parula'
        estimates
        sources
    end
    methods

        function plotHeadModel(obj,varargin)
        % 
        %   plotHeadModel(name1,value1,name2,value2...)
        %
        % Arguments:
        %   model       - path to model
        %   axis        - axis to cut along if desired e.g. 'x' or 'y'
        %   cut         - how far along would you like to cut
        %   fill        - how would you like to colour in the tetrahedrons
        %   solution    - solution to the FEM equations.
        %   electrodes  - (boolean) would you like the electrodes plotted 
        %                 on the scalp?
        %   cond_map    - (boolean) for mapping the anisotropy in the
        %                 skull which should be the second layer
        %   sources     - is the cond_map made with sources
        %   sample_num  - sample number to load estimates for.
        %   folder      - (character vector) name of inverse sub-folder, leave
        %                 out to indicate main folder
        %   elec_err    - plot the error of estimation on each electrode,
        %                 give layer number.
        %   tissue      - which tissue to display (a number)
        %   map_tissue  - tissue to interpolate onto skull
        %   amp_dense   - plot the current density - must come with fill
        %                 argument and conds argument
        %   conds       - conductivity of head model for current density
        %   LIC         - plot vector field
        %   thickness   - boolean dispay thickness of skull as colour
        %   cmap        - matlab colormap option
        %   estimates   - the estimates on the electrodes to plot
        %   
        %   
        %   Examples:
        %       
        %       plotting.plotHeadModel('model',model,'electrodes',true,...
        %           'sample_num',4,'folder','inverse_12345','elec_err',3)
        %
        %       plotting.plotHeadModel('model',model,'cond_map',true,...
        %           'sample_num',4,)
        %

            obj = obj.processArgs(varargin);
            obj = obj.processModel();
        
            p = obj.p; t = obj.t; f = obj.f;
            
            if ~isempty(obj.axis) && isempty(obj.cut)
                error('Please provide cut name value pair')
            elseif ~isempty(obj.cut) && isempty(obj.axis)
                error('Please provide axis name value pair')
            end
            
            if ~isempty(obj.axis)
                switch obj.axis
                    case 'x'
                        obj.expr = ['p(:,1)>' num2str(obj.cut)];
                    case 'y'
                        obj.expr = ['p(:,2)>' num2str(obj.cut)];
                    case 'z'
                        obj.expr = ['p(:,3)>' num2str(obj.cut)];
                end
            end
            
            if ~isempty(obj.fill) 
                if isvector(obj.fill) && isempty(obj.amp_dense)
                    d = obj.fill;
                elseif ~isempty(obj.amp_dense) && ~isempty(obj.fill) && ~isempty(obj.conds)
                    D = zeros(size(t,1),6);
                    for i=1:size(obj.conds,2)
                        Ind = (t(:,5)==i);D(Ind,[1,4,6])=obj.conds(i);
                    end
                    J = EIT_current_density(p,t,D,obj.fill);
                    d = J;
                    if ~isempty(obj.LIC)

                        m=min(p);M=max(p);h=4e2;step=(M(1)-m(1))/h;
                        s1=m(1):step:M(1);
                        s3=m(3):step:M(3);
                        [px3,pz3]=meshgrid(s1,s3);

                        pmt=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:)+p(t(:,4),:))/4; % compute centroid

                        F2_m=scatteredInterpolant(pmt(:,1),pmt(:,2),pmt(:,3),J(:,1));
                        qy_m=F2_m(px3,0.1*ones(size(px3)),pz3);
                        T=tsearchn(p,t(:,1:4),[px3(:),0.1*ones(numel(px3),1),pz3(:)]);
                        qy_m(isnan(T))=0;
                        clear F2_m

                        F3_m=scatteredInterpolant(pmt(:,1),pmt(:,2),pmt(:,3),J(:,3));
                        qz_m=F3_m(px3,0.1*ones(size(px3)),pz3);
                        qz_m(isnan(T))=0;
                        clear F3_m

                        figure()
                        plotvfield(qy_m',qz_m',2.2,2,jet,1,'ef',1);
                        return
                    end
                else
                    switch obj.fill
                        case 'x'
                            d = p(:,1);
                        case 'y'
                            d = p(:,2);
                        case 'z'
                            d = p(:,3);
                    end
                end
            elseif ~isempty(obj.thickness) && (obj.tissue==2)
                tree = getenv("ROMEG");
                load([tree '/Models/Real/thickness.mat'],'thickness')
                d = thickness;
            else
                d = p(:,1);
            end
            
            if ~isempty(obj.elec_err)
                OrderedModelClass.changePath(['Result' num2str(obj.sample_num)]); obj.top = getenv("ROMEG_TOP");
                if ~isempty(obj.folder)
                    disp('Loading estimates and synthetic conductivities')
                    load([obj.top '/Results/inverse/ROM/' obj.folder '/estimates.mat'],'estimates')
                    load([obj.top '/Results/measurements/prep.mat'],'Data')
                    layers=find(logical(Data.mu_max-Data.mu_min));
                    obj.conds = Data.synth_cond(layers);
                    obj.ROM_RE = abs(estimates-obj.conds)./obj.conds;
                    for ii=1:size(obj.ROM_RE,2)
                        obj.ROM_RE(:,ii) = normalize(obj.ROM_RE(:,ii),"range");
                    end
                end
            end
            
            if ~isempty(obj.cond_map) && (obj.cond_map==true)
                obj = obj.plotMap();
                d = obj.d;
                t = obj.t;
            end
            
            if ~isempty(obj.tissue)
                t = obj.t(obj.t(:,5)==obj.tissue,1:4);
            end
            
            figure()
            
            tri1=surftri(p,t);
            icol=.1*ones(1,3);
            t=t(:,1:4);tor=t;

            if size(t)<1,return,end

            if ~isempty(obj.expr)
                incl=find(eval(obj.expr));
                t2=t(any(ismember(t(:,1:4),incl),2),:);
                tri1=tri1(any(ismember(tri1,incl),2),:);
                tri2=surftri(p,t2);
                tri2=setdiff(tri2,tri1,'rows');
            else
                tri2=tri1;
            end
                
            if size(d,1)==size(p,1)
                s.FaceVertexCData=d(:);s.Faces=tri1;
            elseif size(d,1)==size(t,1)
%                 ty=tsearchn(p,double(tor),(p(tri2(:,1),:)+p(tri2(:,2),:)+p(tri2(:,3),:))/3);
%                 s.FaceVertexCData=d(~isnan(ty));
%                 s.Faces=tri2;
%                 ty(~isnan(ty)==1)=ty(1);
%                 s.FaceVertexCData=d(~isnan(ty));
                s.CData = d(:);
            else
                warning('Error');return;
            end

            if ~isempty(obj.expr) && size(d,1)==size(p,1)
                h=trimesh(tri2,p(:,1),p(:,2),p(:,3),s.FaceVertexCData);
                set(h,'facecolor',icol,'edgecolor','k');
                hold on
            end

            s.Vertices=p;
            V=patch(s);

            if size(d(:),1)==size(p,1)
                shading interp
            else
                set(V,'FaceColor','flat','FaceAlpha','1','EdgeAlpha','0.3'); 
            end

            if obj.electrodes
                if nargin>3 && ~isempty(obj.expr)
                    incl=find(eval(obj.expr));
                    f2=f(any(ismember(f(:,1:3),incl),2),:);
                else
                    f2=f;
                end

                array = unique(f2(:,4));
                array = array(2:end);

                for ii=array'
                    el_node_ind=f2(:,4) == ii;
                    el_node = f2(el_node_ind,1:3);
                    for kk=1:length(el_node(:,1))
                        verts = [p(el_node(kk,1),:); p(el_node(kk,2),:); p(el_node(kk,3),:)];
                        faces = [1 2 3];
                        if ~isempty(obj.elec_err) && ~(ii==array(end))
                            color = [1 obj.ROM_RE(ii,obj.elec_err) obj.ROM_RE(ii,obj.elec_err)];
                            patch('Faces',faces,'Vertices',verts,'FaceColor',color,'EdgeColor','none')
                        else
                            patch('Faces',faces,'Vertices',verts,'FaceColor','red','EdgeColor','none')
                        end
                    end
                    index = round(length(el_node(:,1))/2);
                    text_node = el_node(index,1);
                    text_node_pos = p(text_node,:) + p(text_node,:)/20;
                    text(text_node_pos(1),text_node_pos(2),text_node_pos(3),num2str(ii,'%.0f'),'FontSize' ,10)
                end
            end
            axis image
            axis off
            %disp(obj.cmap)
            colormap(obj.cmap)
            
            if ~isempty(obj.cond_map)
                colorbar('FontSize',14)
            end

            function tri=surftri(p,t)
            %SURFTRI Find surface triangles from tetrahedra mesh
            %   TRI=SURFTRI(P,T)

            %   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.

            % Form all faces, non-duplicates are surface triangles

            faces=[t(:,[1,2,3]);
                   t(:,[1,2,4]);
                   t(:,[1,3,4]);
                   t(:,[2,3,4])];
            node4=[t(:,4);t(:,3);t(:,2);t(:,1)];
            faces=sort(faces,2);
            [~,ix,jx]=unique(faces,'rows');
            vec=histc(jx,1:max(jx));
            qx=find(vec==1);
            tri=faces(ix(qx),:);
            node4=node4(ix(qx));

            % Orientation
            v1=p(tri(:,2),:)-p(tri(:,1),:);
            v2=p(tri(:,3),:)-p(tri(:,1),:);
            v3=p(node4,:)-p(tri(:,1),:);
            ix=find(dot(cross(v1,v2,2),v3,2)>0);
            tri(ix,[2,3])=tri(ix,[3,2]);
            end
        end
        
        function obj = plotMap(obj)
            
            OrderedModelClass.changePath(['Result' num2str(obj.sample_num)]); obj.top = getenv("ROMEG_TOP");
            
            %load sinks
            obj = obj.loadSinks();
            electrodes = obj.sinks(:,1);
            
            %load estimates
            
            if isempty(obj.sources) || obj.sources == false
                if isempty(obj.estimates)
                    if isempty(obj.folder)
                        for ii=1:length(electrodes)
                            load([obj.top '/Results/inverse/ROM/inv_' num2str(ii) '.mat'],'inv')
                            estimates = [estimates; inv.estimate];
                        end
                        obj.estimates = estimates(:,obj.map_tissue);
                    else
                        load([obj.top '/Results/inverse/ROM/' obj.folder '/estimates.mat'],'estimates')
                        obj.estimates = estimates(:,obj.map_tissue);
                    end
                end

                t2 = obj.t(obj.t(:,5)==obj.map_tissue,1:4);

                elec_centers = zeros(size(obj.estimates,1),3);
                for ii = 1:length(electrodes)
                    nodes = obj.f(obj.f(:,4)==electrodes(ii),1:3);
                    nodes = unique(nodes);
                    elec_centers(ii,:) = [mean(obj.p(nodes,1)) mean(obj.p(nodes,2)) mean(obj.p(nodes,3))];
                end
            else
                %tbd
            end

            indx = unique(t2(:,1:4));
            pp = obj.p(indx,:);
            PQ = dsearchn(pp,elec_centers);
            F = scatteredInterpolant(pp(PQ,:),obj.estimates);
            obj.VQ = F(obj.p);
            obj.d = obj.VQ;
            obj.t=t2;
        end
    end
end