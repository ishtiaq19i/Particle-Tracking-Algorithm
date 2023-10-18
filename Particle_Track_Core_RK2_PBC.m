%% Particle Track Code

% x,y,z,u1,v1,w1,x0,y0,z0,n,r,delta_t,method_interp are
% needed as input
% x - grid points in x
% y - grid points in y
% z - grid points in z
% u1 - velocity in x
% v1 - velocity in y
% w1 - velocity in z
% x0, y0, z0 - initial points
% n - number of particles
% r - the number of time steps
% delta_t - interval between consecutive timesteps
% method_interp - interpolation method used to interpolate velocities at
% non-nodal points where velocity data is unavailable

function[x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,rows_outofDomain_yz, u_last, v_last, w_last] = Particle_Track_Core_RK2_PBC(x,y,z,u1_t1,u1_t2,v1_t1,v1_t2,w1_t1,w1_t2,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof_t2,rows_outofDomain_yz,nth_time,number_of_timesteps,u_last,v_last,w_last,vof_t1)
X = x;
Y = y;
Z = z;
Lx = max(x);
Ly = max(y);
Lz = max(z);

parfor j=1:n

    % interpolate velocities at initial points
    u_interp1 = interp3(X,Y,Z,u1_t1,x0(j),y0(j),z0(j)); %linear interpolatio at bubble-free location
    v_interp1 = interp3(X,Y,Z,v1_t1,x0(j),y0(j),z0(j));
    w_interp1 = interp3(X,Y,Z,w1_t1,x0(j),y0(j),z0(j));

    u_trajectory1 = u_interp1; % velocity responsible for the trajectory
	v_trajectory1 = v_interp1;
    w_trajectory1 = w_interp1;
    
    % PREDICTOR STEP
	x_track = x0(j) + u_interp1*delta_t;
    y_track = y0(j) + v_interp1*delta_t;
    z_track = z0(j) + w_interp1*delta_t;
    
    % Periodic Boundary Conditions

    x_bound = nan; %nan values if they are not outside of boundary and warping does not take place
    y_bound = nan;
    z_bound = nan;
    x_bound2 = nan; %nan values if they are not outside of boundary and warping does not take place
    y_bound2 = nan;
    z_bound2 = nan;

    if x_track > max(x) || x_track < min(x) || y_track > max(y) || y_track < min(x) || z_track > max(z) || z_track < min(x)
        
        if x_track > max(x)
            x_bound = round(Lx*330/80); %check to see if the particle wraps through a bubble
            x_track = x_track - (Lx-min(x));
            x0(j) = x0(j) - (Lx-min(x));
            y_bound = round(y_track*330/80);
            z_bound = round(z_track*330/80);
            x_bound2 = round(min(x)*330/80);
            y_bound2 = y_bound;
            z_bound2 = z_bound;
        elseif x_track < min(x)
            x_bound = round(min(x)*330/80); %check to see if the particle wraps through a bubble
            x_track = x_track + (Lx-min(x));
            x0(j) = x0(j) + (Lx-min(x));
            y_bound = round(y_track*330/80);
            z_bound = round(z_track*330/80);
            x_bound2 = round(Lx*330/80);
            y_bound2 = y_bound;
            z_bound2 = z_bound;
        end
 
        if y_track > max(y)
            y_bound = round(Ly*330/80); %check to see if the particle wraps through a bubble
            y_track = y_track - (Ly-min(y));
            y0(j) = y0(j) - (Ly-min(y));
            rows_outofDomain_yz(j) = 1;
            x_bound = round(x_track*330/80);
            z_bound = round(z_track*330/80);
            y_bound2 = round(min(y)*330/80);
            x_bound2 = x_bound;
            z_bound2 = z_bound;
        elseif y_track < min(y)
            y_bound = round(min(y)*330/80); %check to see if the particle wraps through a bubble
            y_track = y_track + (Ly-min(y));
            y0(j) = y0(j) + (Ly-min(y));
            rows_outofDomain_yz(j) = 1;
            x_bound = round(x_track*330/80);
            z_bound = round(z_track*330/80);
            y_bound2 = round(Ly*330/80);
            x_bound2 = x_bound;
            z_bound2 = z_bound;
        end
 
        if z_track > max(z)
            z_bound = round(Lz*330/80); %check to see if the particle wraps through a bubble
            z_track = z_track - (Lz-min(z));
            z0(j) = z0(j) - (Lz-min(z));
            rows_outofDomain_yz(j) = 1; 
            x_bound = round(x_track*330/80);
            y_bound = round(y_track*330/80);
            z_bound2 = round(min(z)*330/80);
            x_bound2 = x_bound;
            y_bound2 = y_bound;
        elseif z_track < min(z)
            z_bound = round(min(z)*330/80); %check to see if the particle wraps through a bubble
            z_track = z_track + (Lz-min(z));
            z0(j) = z0(j) + (Lz-min(z));
            rows_outofDomain_yz(j) = 1;  
            x_bound = round(x_track*330/80);
            y_bound = round(y_track*330/80);
            z_bound2 = round(Lz*330/80);
            x_bound2 = x_bound;
            y_bound2 = y_bound;
        end
%         x_track = nan; % if PBC does not exist and domain is finite
%         y_track = nan;
%         z_track = nan;
    else
    end
        
    x_track_temp = x_track; %positions do not know if bubbles at the boundary overlapped with them on their way during warp
    y_track_temp = y_track;
    z_track_temp = z_track;

    % checking whether particles move through bubbles at the boundaries of
    % the domain to warp back to the other side during rk2

    %x_bound is non-nan only when warp happens

    if isnan(y_bound) && isnan(x_bound) && isnan(z_bound)
%         disp('Warping free of bubble intervention at boundary (during RK2)');    
    else
        if vof_t1(y_bound,x_bound,z_bound) == 1 || vof_t2(y_bound,x_bound,z_bound) == 1 || vof_t1(y_bound2,x_bound2,z_bound2) == 1 || vof_t2(y_bound2,x_bound2,z_bound2) == 1
            % Bubble Termination at boundary given x_bound is non nan

        x_track_temp = nan;
        y_track_temp = nan;
        z_track_temp = nan;
        disp('Warping CANNOT happen through the bubbles at the boundaries (during RK2)');
        else
        end
    end
    
    x_track = x_track_temp; %positions have cleared away from bubbles at the boundary
    y_track = y_track_temp;
    z_track = z_track_temp;

    % vof mask - if warping does happen, check if the new warped locations
    % overlap with bubbles
    x_trajectory_nodal = round((165/40)*x_track); %in pixel coordinates
    y_trajectory_nodal = round((165/40)*y_track);
    z_trajectory_nodal = round((165/40)*z_track);

    % indices should be within the domain
    if x_trajectory_nodal > 330
        x_trajectory_nodal = 330;
    elseif x_trajectory_nodal < 1
        x_trajectory_nodal = 1;
    end
    if y_trajectory_nodal > 165
        y_trajectory_nodal = 165;
    elseif y_trajectory_nodal < 1
        y_trajectory_nodal = 1;
    end
    if z_trajectory_nodal > 165
        z_trajectory_nodal = 165;
    elseif z_trajectory_nodal < 1
        z_trajectory_nodal = 1;
    end

    % check to see if the predicted position is safe from bubbles
    if isnan(y_trajectory_nodal) || isnan(x_trajectory_nodal) || isnan(z_trajectory_nodal)
        x_track = nan;
        y_track = nan;
        z_track = nan;
    else
        if vof_t2(y_trajectory_nodal,x_trajectory_nodal,z_trajectory_nodal) == 1 % Bubble Termination

        x_track = nan;
        y_track = nan;
        z_track = nan;
%         disp('Bubble-terminated Particle(s) during RK2');
        else
        end
    end

    x_trajectory_temp = x_track; % saving trajectory values to prevent zero-ing of data due to using parfor loop
    y_trajectory_temp = y_track;
    z_trajectory_temp = z_track;
    
    % velocity interpolation at the predicted step free of bubbles
    u_interp = interp3(X,Y,Z,u1_t2,x_trajectory_temp,y_trajectory_temp,z_trajectory_temp);
    v_interp = interp3(X,Y,Z,v1_t2,x_trajectory_temp,y_trajectory_temp,z_trajectory_temp);
    w_interp = interp3(X,Y,Z,w1_t2,x_trajectory_temp,y_trajectory_temp,z_trajectory_temp);
   
    % CORRECTOR STEP
    x_trajectory_pseudo = x0(j) + (u_trajectory1+u_interp)*delta_t/2;
    y_trajectory_pseudo = y0(j) + (v_trajectory1+v_interp)*delta_t/2;
    z_trajectory_pseudo = z0(j) + (w_trajectory1+w_interp)*delta_t/2;

    % Periodic Boundary Conditions after RK2

    x_track = x_trajectory_pseudo;
    y_track = y_trajectory_pseudo;
    z_track = z_trajectory_pseudo;

    x_bound = nan; %nan values if they are not outside of boundary and warping does not take place
    y_bound = nan;
    z_bound = nan;
    x_bound2 = nan; %nan values if they are not outside of boundary and warping does not take place
    y_bound2 = nan;
    z_bound2 = nan;

    if x_track > max(x) || x_track < min(x) || y_track > max(y) || y_track < min(x) || z_track > max(z) || z_track < min(x)
        
        if x_track > max(x)
            x_bound = round(Lx*330/80); %check to see if the particle wraps through a bubble
            x_track = x_track - (Lx-min(x));
            y_bound = round(y_track*330/80);
            z_bound = round(z_track*330/80);
            x_bound2 = round(min(x)*330/80);
            y_bound2 = y_bound;
            z_bound2 = z_bound;
        elseif x_track < min(x)
            x_bound = round(min(x)*330/80); %check to see if the particle wraps through a bubble
            x_track = x_track + (Lx-min(x));
            y_bound = round(y_track*330/80);
            z_bound = round(z_track*330/80);
            x_bound2 = round(Lx*330/80);
            y_bound2 = y_bound;
            z_bound2 = z_bound;
        end
 
        if y_track > max(y)
            y_bound = round(Ly*330/80); %check to see if the particle wraps through a bubble
            y_track = y_track - (Ly-min(y));
            rows_outofDomain_yz(j) = 1;
            x_bound = round(x_track*330/80);
            z_bound = round(z_track*330/80);
            y_bound2 = round(min(y)*330/80);
            x_bound2 = x_bound;
            z_bound2 = z_bound;
        elseif y_track < min(y)
            y_bound = round(min(y)*330/80); %check to see if the particle wraps through a bubble
            y_track = y_track + (Ly-min(y));
            rows_outofDomain_yz(j) = 1;
            x_bound = round(x_track*330/80);
            z_bound = round(z_track*330/80);
            y_bound2 = round(Ly*330/80);
            x_bound2 = x_bound;
            z_bound2 = z_bound;
        end
 
        if z_track > max(z)
            z_bound = round(Lz*330/80); %check to see if the particle wraps through a bubble
            z_track = z_track - (Lz-min(z));
            rows_outofDomain_yz(j) = 1; 
            x_bound = round(x_track*330/80);
            y_bound = round(y_track*330/80);
            z_bound2 = round(min(z)*330/80);
            x_bound2 = x_bound;
            y_bound2 = y_bound;
        elseif z_track < min(z)
            z_bound = round(min(z)*330/80); %check to see if the particle wraps through a bubble
            z_track = z_track + (Lz-min(z));
            rows_outofDomain_yz(j) = 1;  
            x_bound = round(x_track*330/80);
            y_bound = round(y_track*330/80);
            z_bound2 = round(Lz*330/80);
            x_bound2 = x_bound;
            y_bound2 = y_bound;
        end
%         x_track = nan; % if PBC does not exist and domain is finite
%         y_track = nan;
%         z_track = nan;
    else
    end

    x_trajectory_pseudo = x_track;
    y_trajectory_pseudo = y_track;
    z_trajectory_pseudo = z_track;
   
    % checking whether particles move through bubbles at the boundaries of
    % the domain to warp back to the other side after RK2
    % The function should output positions free of bubble intervention

    % vof_t1 and vof_t2 are used because warping is instantaneous thus
    % final positions exist in two timesteps
    if isnan(y_bound) && isnan(x_bound) && isnan(z_bound)
%         disp('Warping free of bubble intervention at boundary (after RK2)');    
    else
        if vof_t1(y_bound,x_bound,z_bound) == 1 || vof_t2(y_bound,x_bound,z_bound) == 1 || vof_t1(y_bound2,x_bound2,z_bound2) == 1 || vof_t2(y_bound2,x_bound2,z_bound2) == 1% Bubble Termination at boundary given x_bound is non nan

        x_trajectory_pseudo = nan;
        y_trajectory_pseudo = nan;
        z_trajectory_pseudo = nan;
        disp('Warping CANNOT happen through the bubbles at the boundaries (after RK2)');
        else
        end
    end

    % After warp checking at boundary, the final positions must check
    % for bubble intervention
    x_track_nodal = round((165/40).*x_trajectory_pseudo); %in pixel coordinates
    y_track_nodal = round((165/40).*y_trajectory_pseudo);
    z_track_nodal = round((165/40).*z_trajectory_pseudo);
    if x_track_nodal > 330
        x_track_nodal = 330;
    elseif x_track_nodal < 1
        x_track_nodal = 1;
    end
    if y_track_nodal > 165
        y_track_nodal = 165;
    elseif y_track_nodal < 1
        y_track_nodal = 1;
    end
    if z_track_nodal > 165
        z_track_nodal = 165;
    elseif z_track_nodal < 1
        z_track_nodal = 1;
    end

% if particles survive bubbles at boundaries...
    if isnan(y_track_nodal) || isnan(x_track_nodal) || isnan(z_track_nodal)
        x_trajectory_pseudo = nan;
        y_trajectory_pseudo = nan;
        z_trajectory_pseudo = nan;
    else
        if vof_t2(y_track_nodal,x_track_nodal,z_track_nodal) == 1 % Bubble Termination

        x_trajectory_pseudo = nan;
        y_trajectory_pseudo = nan;
        z_trajectory_pseudo = nan;
%         disp('Bubble-terminated Particle(s) after RK2:');
        else
        end
    end

    x_trajectory(j) = x_trajectory_pseudo; %final position after several checks
    y_trajectory(j) = y_trajectory_pseudo;
    z_trajectory(j) = z_trajectory_pseudo;

    u_trajectory(j) = u_trajectory1;
    v_trajectory(j) = v_trajectory1;
    w_trajectory(j) = w_trajectory1;

    % interpolating lagrangian velocity for the last point
    u_interp2 = interp3(X,Y,Z,u1_t2,x_trajectory(j),y_trajectory(j),z_trajectory(j));
    v_interp2 = interp3(X,Y,Z,v1_t2,x_trajectory(j),y_trajectory(j),z_trajectory(j));
    w_interp2 = interp3(X,Y,Z,w1_t2,x_trajectory(j),y_trajectory(j),z_trajectory(j));

% recording the velocity data on the last timestep
    if nth_time < number_of_timesteps
        u_last_temp = u_interp2;
        v_last_temp = v_interp2;
        w_last_temp = w_interp2;
    else
        u_last_temp = nan;
        v_last_temp = nan;
        w_last_temp = nan;
    end
    u_last(j) = u_last_temp;
    v_last(j) = v_last_temp;
    w_last(j) = w_last_temp;
end % end parfor

end

