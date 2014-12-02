function [angvel,T] = orientrateToVel(varargin)

orient=varargin{1};
if numel(orient)==3
    [angvel,T] = eulerrateToVel(varargin{:});
elseif numel(orient)==4
    [angvel,T] = quatrateToVel(varargin{:});
end
