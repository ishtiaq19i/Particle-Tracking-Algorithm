
function [dudx_eL, dudy_eL, dudz_eL, dvdx_eL, dvdy_eL, dvdz_eL, dwdx_eL, dwdy_eL, dwdz_eL] = GetVelocityGradientEulcondLag(x0, y0, z0, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)

tic
x0_nodal = round((165/40).*x0); %in pixel coordinates
y0_nodal = round((165/40).*y0);
z0_nodal = round((165/40).*z0);

if x0_nodal > 330
    x0_nodal = 330;
elseif x0_nodal < 1
    x0_nodal = 1;
end

if y0_nodal > 165
    y0_nodal = 165;
elseif y0_nodal < 1
    y0_nodal = 1;
end

if z0_nodal > 165
    z0_nodal = 165;
elseif z0_nodal < 1
    z0_nodal = 1;
end

n = numel(x0);
dudx_eL = nan(n,1);
dudy_eL = nan(n,1);
dudz_eL = nan(n,1);
dvdx_eL = nan(n,1);
dvdy_eL = nan(n,1);
dvdz_eL = nan(n,1);
dwdx_eL = nan(n,1);
dwdy_eL = nan(n,1);
dwdz_eL = nan(n,1);

parfor i=1:numel(x0)
    
    if isnan(y0_nodal(i)) || isnan(x0_nodal(i)) || isnan(z0_nodal(i))
        dudx_eL_temp = nan;
        dudy_eL_temp = nan;
        dudz_eL_temp = nan;
        dvdx_eL_temp = nan;
        dvdy_eL_temp = nan;
        dvdz_eL_temp = nan;
        dwdx_eL_temp = nan;
        dwdy_eL_temp = nan;
        dwdz_eL_temp = nan;
    else
        dudx_eL_temp = dudx(y0_nodal(i),x0_nodal(i),z0_nodal(i));
        dudy_eL_temp = dudy(y0_nodal(i),x0_nodal(i),z0_nodal(i));
        dudz_eL_temp = dudz(y0_nodal(i),x0_nodal(i),z0_nodal(i));
        dvdx_eL_temp = dvdx(y0_nodal(i),x0_nodal(i),z0_nodal(i));
        dvdy_eL_temp = dvdy(y0_nodal(i),x0_nodal(i),z0_nodal(i));
        dvdz_eL_temp = dvdz(y0_nodal(i),x0_nodal(i),z0_nodal(i));
        dwdx_eL_temp = dwdx(y0_nodal(i),x0_nodal(i),z0_nodal(i));
        dwdy_eL_temp = dwdy(y0_nodal(i),x0_nodal(i),z0_nodal(i));
        dwdz_eL_temp = dwdz(y0_nodal(i),x0_nodal(i),z0_nodal(i));
    end

        dudx_eL(i) = dudx_eL_temp;
        dudy_eL(i) = dudy_eL_temp;
        dudz_eL(i) = dudz_eL_temp;
        dvdx_eL(i) = dvdx_eL_temp;
        dvdy_eL(i) = dvdy_eL_temp;
        dvdz_eL(i) = dvdz_eL_temp;
        dwdx_eL(i) = dwdx_eL_temp;
        dwdx_eL(i) = dwdy_eL_temp;
        dwdz_eL(i) = dwdz_eL_temp;
end
toc
end