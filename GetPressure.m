%%%% Get Pressure using lagrangian particle trajectories

function [P] = GetPressure(p1,x0,y0,z0,x,y,z)

n = numel(x0);

X = x;
Y = y;
Z = z;

p_interp = nan([n 1]);
parfor j = 1:n
    p_interp(j) = interp3(X,Y,Z,p1,x0(j),y0(j),z0(j));
end

P = p_interp;