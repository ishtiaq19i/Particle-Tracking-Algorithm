%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finite_diff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes gradients of F(x,y,z) on a Cartesian grid (x,y,z). 
% The user has the option to specify the FD order used in the computation of 
% gradients at interior points. Gradients at the boundaries are always 
% computed with a 1st-order forward/backward difference. 
%
% 14-July-2018
% 1-August-2018 (Added 4th order central differencing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dkdx, dkdy, dkdz] = Finite_diff_3Ddz(k, method, dx, dy, dz)

[nrow,ncol,N] = size(k);

switch method
    case 1
        %Forward finite differencing
        disp('1st order differencing')
        dkdx = diff(k,1,2)/(dx); dkdx = [dkdx NaN(nrow,1,N)]; %columnwise derivative
        dkdx(:,end,:) = (k(:,end,:)-k(:,end-1,:))/(dx);       %Backward difference
        
        dkdy = diff(k,1,1)/(dy); dkdy = [dkdy;NaN(1,ncol,N)]; %row-wise derivative
        dkdy(end,:,:) = (k(end,:,:)-k(end-1,:,:))/(dy);       %Backward difference
        
        dkdz = diff(k,1,3)/(dz); dkdz = cat(3,dkdz,NaN(nrow,ncol,1)); %depth-wise derivative
        dkdz(:,:,end) = (k(:,:,end)-k(:,:,end-1))./dz(:,:,end);       %Backward difference
        
    case 2
        %Central finite differencing
        disp('2nd order CD')
        A1 = [NaN(nrow,1,N) k]; A2 = [k NaN(nrow,1,N)];
        dkdx = (diff(A1,1,2) + diff(A2,1,2))./(2*dx);           %columnwise derivative
        dkdx(:,1,:) = (k(:,2,:)-k(:,1,:))./dx(:,1,:);           %Forward difference
        dkdx(:,end,:) = (k(:,end,:)-k(:,end-1,:))./dx(:,end,:); %Backward difference
        
        A1 = [NaN(1,ncol,N); k]; A2 = [k; NaN(1,ncol,N)];
        dkdy = (diff(A1,1,1) + diff(A2,1,1))./(2*dy);           %row-wise derivative
        dkdy(1,:,:) = (k(2,:,:)-k(1,:,:))./dy(1,:,:);           %Forward difference
        dkdy(end,:,:) = (k(end,:,:)-k(end-1,:,:))./dy(end,:,:); %Backward difference
        
        A1 = cat(3,NaN(nrow,ncol,1), k); A2 = cat(3,k,NaN(nrow,ncol,1));
        dkdz = (diff(A1,1,3) + diff(A2,1,3))./(2*dz);           %depth-wise derivative
        dkdz(:,:,1) = (k(:,:,2)-k(:,:,1))./dz(:,:,1);           %Forward difference
        dkdz(:,:,end) = (k(:,:,end)-k(:,:,end-1))./dz(:,:,end); %Backward difference
        
%         A1 = cat(3,NaN(nrow,ncol,1), k);
%         dkdz = diff(A1,1,3)./(dz);                                    %depth-wise derivative
%         dkdz(:,:,end) = (k(:,:,end)-k(:,:,end-1))./dz(:,:,end);       %Backward difference
        
    case 3
        %4th order central differencing 
        disp('4th order CD')
        A1 = [NaN(nrow,1,N) k]; A2 = [k NaN(nrow,1,N)];
        A3 = [NaN(nrow,2,N) k(:,1:ncol-2,:)]; A4 = [k(:,3:ncol,:) NaN(nrow,2,N)];
        dkdx = (2/3)*(diff(A1,1,2) + diff(A2,1,2))./(dx) - (1/12)*(A4-A3)./(dx); %columnwise derivative
        dkdx(:,2,:) = (k(:,3,:)-k(:,1,:))./(2*dx(:,2,:));                        %Central difference
        dkdx(:,end-1,:) = (k(:,end,:)-k(:,end-2,:))./(2*dx(:,end-1,:));          %Central difference
        dkdx(:,1,:) = (k(:,2,:)-k(:,1,:))./(dx(:,1,:));                          %Forward difference
        dkdx(:,end,:) = (k(:,end,:)-k(:,end-1,:))./(dx(:,end,:));                %Backward difference
        
        A1 = [NaN(1,ncol,N); k]; A2 = [k; NaN(1,ncol,N)];
        A3 = [NaN(2,ncol,N); k(1:nrow-2,:,:)]; A4 = [k(3:nrow,:,:); NaN(2,ncol,N)];
        dkdy = (2/3)*(diff(A1,1,1) + diff(A2,1,1))./(dy) - (1/12)*(A4-A3)./(dy); %row-wise derivative
        dkdy(2,:,:) = (k(3,:,:)-k(1,:,:))./(2*dy(2,:,:));                        %Central difference
        dkdy(end-1,:,:) = (k(end,:,:)-k(end-2,:,:))./(2*dy(end-1,:,:));          %Central difference
        dkdy(1,:,:) = (k(2,:,:)-k(1,:,:))./(dy(1,:,:));                          %Forward difference
        dkdy(end,:,:) = (k(end,:,:)-k(end-1,:,:))./(dy(end,:,:));                %Backward difference
        
        A1 = cat(3,NaN(nrow,ncol,1), k); A2 = cat(3,k,NaN(nrow,ncol,1));
        A3 = cat(3,NaN(nrow,ncol,2), k(:,:,1:N-2)); A4 = cat(3,k(:,:,3:N),NaN(nrow,ncol,2));
        dkdz = (2/3)*(diff(A1,1,3) + diff(A2,1,3))./(dz) - (1/12)*(A4-A3)./(dz); %depth-wise derivative
        dkdz(:,:,2) = (k(:,:,3)-k(:,:,1))./(2*dz(:,:,2));                        %Central difference
        dkdz(:,:,end-1) = (k(:,:,end)-k(:,:,end-2))./(2*dz(:,:,end-1));          %Central difference
        dkdz(:,:,1) = (k(:,:,2)-k(:,:,1))./(dz(:,:,1));                          %Forward difference
        dkdz(:,:,end) = (k(:,:,end)-k(:,:,end-1))./(dz(:,:,end));                %Backward difference

end

end

