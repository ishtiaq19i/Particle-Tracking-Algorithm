function[e] = EnergyDissipation(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vis_kin,e)

% dudx free of bubbles
Nx = numel(dudx(1,:,1));
Ny = numel(dudx(:,1,1));
Nz = numel(dudx(:,1,1));

dudx = dudx(:);
dudy = dudy(:);
dudz = dudz(:);
dvdx = dvdx(:);
dvdy = dvdy(:);
dvdz = dvdz(:);
dwdx = dwdx(:);
dwdy = dwdy(:);
dwdz = dwdz(:);

p = numel(dvdx);

parfor i = 1:p % number of points
    % strain rate components
    s_xx = 0.5*(dudx(i) + dudx(i));
    s_xy = 0.5*(dudy(i) + dvdx(i));
    s_xz = 0.5*(dudz(i) + dwdx(i));
    s_yy = 0.5*(dvdy(i) + dvdy(i));
    s_yx = 0.5*(dvdx(i) + dudy(i));
    s_yz = 0.5*(dvdz(i) + dwdy(i));
    s_zz = 0.5*(dwdz(i) + dwdz(i));
    s_zx = 0.5*(dwdx(i) + dudz(i));
    s_zy = 0.5*(dwdy(i) + dvdz(i));

    % strain rate
    s_ij = [s_xx s_xy s_xz; s_yx  s_yy s_yz; s_zx s_zy s_zz];
    % Lagrangian Mechanical Energy Dissipation in mm2/s3
    e(i) = 2*vis_kin*(sum(s_ij.^2,'all'));

end

e = mean(e,'omitnan'); %average energy dissipation
end