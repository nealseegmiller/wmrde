function [V,I] = rock(nv,roundness,L,dsplit,npert)
%make a random rock TriMeshSurf object

% http://markjstock.org/rocktools/

%create_cubic_nodes
V = rand(3,nv) - .5;

%sphericalize_nodes

%DEBUGGING
% V0 = V; 
% I0 = convhull(V(1,:), V(2,:), V(3,:));

for i = 1:round(roundness*50);
    for j = 1:nv    %loop over vertices
        D = V(:,j)*ones(1,nv) - V;
        dist = sqrt(sum(D.^2,1));
        diff = .0075./dist - dist*.01333;
        
        tmp = (ones(3,1)*diff) .* D;
        tmp(:,j) = [];
        tmp = sum(tmp,2);
        V(:,j) = V(:,j) + tmp;
    end
end

%scale
V = V .* (L*ones(1,nv));

I = convhull(V(1,:), V(2,:), V(3,:));

E = TriangleEdgeLengths(V,I);

Ie = [1 2; 2 3; 3 1];


%DEBUGGING
% V0 = V;
% I0 = I;

for iter = 1:dsplit
    nt_ = size(I,1);
    npert = npert/2;
    
    for sno = 1:nt_ %split number

        %which edge is longest?
        [E_,eno] = max(E,[],2);
        [~,tno] = max(E_);
        eno = eno(tno);

        i1 = I(tno,Ie(eno,1));
        i2 = I(tno,Ie(eno,2));

        %which two triangles contain edge?
        tno = [];
        N = []; %normal
        i3 = []; %the index that does not belong to the edge
        
        nt = size(I,1);
        for i = 1:nt
            I_ = I(i,:);
            if any(i1 == I_) && any(i2 == I_)
                tno(end+1) = i;
                N(:,end+1) = TriangleNormal(V(:,I_));
                i3(end+1) = I_(I_ ~= i1 & I_ ~= i2);
            end
        end

        %add vertex
        V(:,end+1) = (V(:,i1) + V(:,i2))/2 + rand*npert*N(:,1) + rand*npert*N(:,2);
        nv = nv+1;

        %add triangles
        Inew = [i3(1) i1 nv;
                i3(1) nv i2;
                i3(2) nv i1;
                i3(2) i2 nv];

        %TODO, necessary?
        N_ = TriangleNormal(V(:,Inew(1,:)));
        if N_'*N(:,1) < 0
            Inew(:,[2 3]) = Inew(:,[3 2]); %make clockwise
        end

        I = [I; Inew];

        E = [E; TriangleEdgeLengths(V,Inew)];

        %remove triangles
        I(tno,:) = [];
        E(tno,:) = [];

    end
end

if 0
    %%
    %PLOT

    set(figure, 'name', 'rock')
    hold on
    
    C = .8*ones(1,3);
    alpha = 1;
    h = patch('Vertices',V','Faces',I,'FaceColor',C,'FaceAlpha',alpha);
    set(h,'EdgeColor','none');
  
    axis equal
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    make_legible(14)

    h_light = light('Position',[1 0 1],'Parent',gca);
    
%     %DEBUGGING
%     C = 'none';
%     h = patch('Vertices',V0','Faces',I0,'FaceColor',C,'FaceAlpha',alpha);
%     set(h,'EdgeColor','r');
    

end

end

function E = TriangleEdgeLengths(V,I)

    Ie = [1 2; 2 3; 3 1];
    
    for i = 1:3
        I1 = I(:,Ie(i,1));
        I2 = I(:,Ie(i,2));
        E(:,i) = sqrt(sum((V(:,I2) - V(:,I1)).^2,1));
    end

end

function N = TriangleNormal(V)
    N = cross3(V(:,2)-V(:,1),V(:,3)-V(:,1));
    N = N/norm(N);
    
end

