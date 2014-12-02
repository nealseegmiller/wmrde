function [h_line,h_surf] = drawEllipse(C,m,nv,conf,marg)
%draw ellipse or ellipsoid
%INPUTS
%C:         (2x2 or 3x3) covariance matrix
%m:         (2x1 or 3x1) mean
%nv:        number of vertices per ellipse
%(optional)
%conf:      confidence (0 < conf < 1), .95 default
%marg:      {1,2,3} dimension to marginalize out or []. Only if C is 3x3
%OUTPUTS
%h_line:    (scalar or 3x1) handle line(s)
%h_surf:    (scalar) handle surf

%If C is 3x3 & marg is empty, ellipsoid is drawn as 3 ellipses
%If C is 3x3 & marg not empty, single ellipse is drawn for marginal distribution

%normcdf(s,0,1) - normcdf(-s,0,1) == chi2cdf(s^2,1)


if nargin < 5
    marg = [];
    if nargin < 4
        conf=.95;
        if nargin < 3
            nv=36*2;
        end
    end
end

dim=size(C,1);
dof=size(C,1); %degrees of freedom for Chi square distribution

if dim==3 && ~isempty(marg)
    %marginalize out one dimension
    %for multivariate normal distribution, simply drop out the dimension from mean and covariance
    %http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Marginal_distributions
    I = 1:3~=marg;
    C = C(I,I);
    m_marg = m(marg);
    m = m(I);
    dim=2;
end

[V,D] = eig(C);
sig = sqrt(chi2inv(conf,dof));
d = sqrt(diag(D))*sig; %semi diameters

%repmat
D = d*ones(1,nv+1);
M = m*ones(1,nv+1);

if dim == 2
    %2D ellipse

    if isempty(marg)
        %unit circle
        pts = pointsOnCircle(1,nv,true);
        %scale
        pts = pts .* D;
        %rotate and add mean
        pts = V*pts + M;

        h_line = plot(pts(1,:), pts(2,:));
    else
        %unit circle
        pts = pointsOnCircle(1,nv,true,marg);
        %scale
        pts(I,:) = pts(I,:) .* D;
        %rotate and add mean
        pts(I,:) = V*pts(I,:) + M;
        pts(marg,:) = m_marg;

        h_line = plot3(pts(1,:), pts(2,:), pts(3,:));
    end
    h_surf = [];
    
elseif dim == 3
    %3D ellipsoid
    
    %plot three 2D ellipses
    h_line = zeros(3,1);

    for i = 1:3;
        %unit circle
        pts = pointsOnCircle(1,nv,true,i);
        %scale
        I = 1:3~=i;
        pts(I,:) = pts(I,:) .* D(I,:);
        %rotate and add mean
        pts = V*pts + M;

        h_line(i) = plot3(pts(1,:), pts(2,:), pts(3,:));
    end

    %%
    
    [X,Y,Z] = ellipsoid(0,0,0,d(1),d(2),d(3),nv);
    pts = [X(:), Y(:), Z(:)]';
    pts = V*pts + m*ones(1,size(pts,2)); %rotate and add mean
    X(:) = pts(1,:);
    Y(:) = pts(2,:);
    Z(:) = pts(3,:);
    h_surf = surf(X,Y,Z,'FaceAlpha',.2,'FaceColor',[1 1 1]/2,'EdgeColor','none');
    
end
    
    