x0 = [7.27272727272727
7.27272727272727
4.84848484848485
2.42424242424242
7.27272727272727
12.1212121212121
2.42424242424242
7.27272727272727
12.1212121212121
4.84848484848485];

y0 = x0;
z0 = x0;

load('dudx.mat')
load('dudy.mat')
load('dudz.mat')
load('dvdx.mat')
load('dvdy.mat')
load('dvdz.mat')
load('dwdx.mat')
load('dwdy.mat')
load('dwdz.mat')
load('vof.mat')

vis_kin = 1;
% [dudx_eL, dudy_eL, dudz_eL, dvdx_eL, dvdy_eL, dvdz_eL, dwdx_eL, dwdy_eL, dwdz_eL] = GetVelocityGradientEulcondLag(x0, y0, z0, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,vof);
Nx = numel(dudx(1,:,1));
Ny = numel(dudx(:,1,1));
Nz = numel(dudx(:,1,1));
e = nan([1 Nx*Ny*Nz]);

[e] = EnergyDissipation(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vis_kin,e);
