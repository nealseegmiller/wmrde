function surface = ramp()

%set the following parameters:
%w:     width of ramp
%h:     height of ramp
%a:     ramp angle
%lf:    length of flat section (may be zero)
%x,y:   location of ramp center
%yaw:   ramp orientation

intom=.0254;

if 0
    %Zoe rover ramp with flat
    w = 24*intom;
    h = 16*intom;
    a = asin(h/(27*intom));
    lf = 24*intom;
    x = -1;
    y = -1;
    yaw = 0*pi/180;

elseif 1
    %Zoe rover ramp with no flat
    w = 24*intom;
    h = 16.5*intom;
    a = asin(h/(30*intom));
    lf = 0;
    x = 3;
    y = -.9;
    yaw = 0*pi/180;

elseif 0
    %for LandTamer
    w = 2;
    h = 2;
    a = 20*pi/180;
    lf = 0;
    x = 10;
    y = 0;
    yaw = 0*pi/180;
elseif 0
    %for Rocky 7
    w = .25;
    h = .15;
    a = 30*pi/180;
    lf = 0;
    x = .65;
    y = -.2;
    yaw = 0*pi/180;
end


ls = h/tan(a); %length of sloped surface
HT = poseToHT([0 0 yaw],[x y 0]);

surface = TriMeshSurf();

if lf > 0
    V = [
         ls+lf/2  w/2 0;
         ls+lf/2 -w/2 0;
         lf/2     w/2 h;
         lf/2    -w/2 h;
        -lf/2     w/2 h;
        -lf/2    -w/2 h;
        -ls-lf/2  w/2 0;
        -ls-lf/2 -w/2 0
    ]';
    V = applyHT(HT,V);
    
    addVertices(surface,V);
    addQuad(surface,[1 3 4 2]);
    addQuad(surface,[3 5 6 4]);
    addQuad(surface,[5 7 8 6]);
else
    V = [
         ls  w/2 0;
         ls -w/2 0;
         0   w/2 h;
         0  -w/2 h;
        -ls  w/2 0;
        -ls -w/2 0
    ]';
    V = applyHT(HT,V);
    
    addVertices(surface,V);
    addQuad(surface,[1 3 4 2]);
    addQuad(surface,[3 5 6 4]);
    
end

return




