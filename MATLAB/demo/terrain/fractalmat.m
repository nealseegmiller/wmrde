function Z = fractalmat(lx,ly,ppm,size0,ni,hmult,hmin,sig)
%INPUTS
%lx:        length in x (m)
%ly:        length in y (m)
%ppm:       points per meter
%size0:     size of zm input to generate_brownian_mesh(), output is size0*(2^ni) + 1?
%ni:        num. iterations of Brownian motion
%hmult:     height multiplier
%hmin:      min. height
%sig:       sigma for gaussian filter (m), [] for no filter

if nargin < 8
    sig = [];
    if nargin < 7
        hmin = -Inf;
        if nargin < 6
            hmult = 1.0;
        end
    end
end

L = max(lx,ly);

Z = generate_brownian_mesh(ni,zeros(size0)); %returns size ?

N = L*ppm + 1;
Z = changeRes(Z,N*[1 1]);

Z = Z*hmult;
Z(Z < hmin) = hmin;

%Gaussian filter
if ~isempty(sig)
    hsize = round_to_odd(10*sig*ppm);
    H = fspecial('gaussian',hsize,sig*ppm);
    Z = conv2(Z,H,'same');
end

%crop if necessary
if lx ~= ly
    Nx = lx*ppm+1;
    Ny = ly*ppm+1;
    Z = Z(1:Nx,1:Ny);
end


end
