
function[x_track_int, y_track_int, z_track_int, u_track_int, v_track_int, w_track_int,rows_outofDomain_yz, u_last, v_last, w_last] = Particle_Track_Intermediate_PBC(x,y,z,u1,v1,w1,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof,rows_outofDomain_yz,nth_time,number_of_timesteps,u_last,v_last,w_last)

u1_t1 = u1(:,:,:,1); 
v1_t1 = v1(:,:,:,1);
w1_t1 = w1(:,:,:,1);

u1_t2 = u1(:,:,:,2);
v1_t2 = v1(:,:,:,2);
w1_t2 = w1(:,:,:,2);

vof_t1 = vof(:,:,:,1);
vof_t2 = vof(:,:,:,2);

% Change to a specific directory
newDirectory = 'C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\Particle_Track_Core_RK2_PBC';
cd(newDirectory);

% Check the current working directory
currentDirectory = pwd;
disp(['Current working directory: ' currentDirectory])
        
tic

[x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,rows_outofDomain_yz, u_last, v_last, w_last] = Particle_Track_Core_RK2_PBC_mex(x,y,z,u1_t1,u1_t2,v1_t1,v1_t2,w1_t1,w1_t2,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof_t2,rows_outofDomain_yz,nth_time,number_of_timesteps,u_last,v_last,w_last,vof_t1);

toc

x_track_int = x_trajectory;
y_track_int = y_trajectory;
z_track_int = z_trajectory;
u_track_int = u_trajectory;
v_track_int = v_trajectory;
w_track_int = w_trajectory;

% Change to a specific directory
newDirectory = 'C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions';
cd(newDirectory);

%Check the current working directory
currentDirectory = pwd;
disp(['Current working directory: ' currentDirectory]);

end
