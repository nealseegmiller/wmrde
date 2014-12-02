function h = drawSurfaces(surfaces,ax)
%ax, axis handle

if nargin < 2
    ax = gca;
end

%PLOT
h = zeros(length(surfaces),1);


for i = 1:length(surfaces)
    %TODO, if PlaneSurf
    
    if isa(surfaces{i},'TriMeshSurf')
        V = surfaces{i}.vertices';
        I = surfaces{i}.indices;
        C = .8*ones(1,3);
%         Ce = .75*ones(1,3);
        Ce = 'k';
        alpha = 1;
        h(i) = patch('Vertices',V,'Faces',I,'FaceColor',C,'FaceAlpha',alpha);
        set(h(i),'EdgeColor',Ce,'Clipping','off');
    end
    
    
    if isa(surfaces{i},'GridSurf')
        
        V=surfaces{i}.Z.Values;
        
        if 0
            
            h(i) = surf(ax, surfaces{i}.X,surfaces{i}.Y,V);
            set(h(i),'EdgeColor','none');
            
            gs=surfaces{i}.dx; %grid spacing (meters per grid index)
            if gs~=surfaces{i}.dy
                disp('WARNING, equal grid spacing in X,Y is recommended')
            end
            ppm = 300;
            apply_texture(h(i), ppm, 'desert_sand.jpg', V, gs);
            
            
        else
            colormap copper
%             colormap gray
            
            cmax = max(V(:));
            cmin = min(V(:));
            
            %DEBUGGING, so not too dark or light
            cmin_ = cmin; %backup
            cmax_ = cmax;
            cmin = cmax_ - 1.25*(cmax_-cmin_); 
            cmax = cmin_ + 1.25*(cmax_-cmin_);
            
            caxis([cmin cmax]); %make_vrml_hgt changes this?
            
            
%             h(i)=mesh(ax, surfaces{i}.X, surfaces{i}.Y, V);
            h(i)=surf(ax, surfaces{i}.X, surfaces{i}.Y, V, 'EdgeColor', 'none');
        end
        
        set(h(i),'Clipping','off')
    end
end

end


function apply_texture(h,ppm,im_filename,V,gs)
    %resolution
    %ppm:   pixels per meter
    
    im=load_texture_image(im_filename,ppm); %get scaled texture image

    C = double(im)/255;

    %tile and crop image to match size of Z
    %assumes image is seamless


    ppi=round(ppm*gs); %pixels per grid index

    [m,n]=size(V);

    C = repmat(C,[ceil(m*ppi/size(C,1)),ceil(n*ppi/size(C,2)),1]); %tile
    C = C(1:m*ppi,1:n*ppi,:); %crop

    if 0
        %SCALE INTENSITY ACCORDING TO HEIGHT
        %helps to see geometry if no light

        %interpolate Z
        S = change_res(V,[m n]*ppi);

        %scale
        range_S = range(S(:));
        median_S = median(S(:));
        S = 1.0/range_S*(S - median_S) + 1;
        S = repmat(S,[1 1 3]);

        C = S.*C;
        C(C < 0) = 0;
        C(C > 1) = 1;
    end

    set(h,'CData',C,'FaceColor','texturemap');


end
