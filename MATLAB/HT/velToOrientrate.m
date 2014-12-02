function [orientrate,T,dTdorient] = velToOrientrate(varargin)
orient = varargin{1};

if nargout > 2
    if length(orient) == 3
        [orientrate,T,dTdorient] = velToEulerrate(varargin{:});
    elseif length(orient) == 4
        [orientrate,T,dTdorient] = velToQuatrate(varargin{:});
    end
else
    if length(orient) == 3
        [orientrate,T] = velToEulerrate(varargin{:});
    elseif length(orient) == 4
        [orientrate,T] = velToQuatrate(varargin{:});
    end
end