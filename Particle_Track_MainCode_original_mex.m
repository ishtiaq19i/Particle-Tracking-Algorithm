tic % measuring the code's elapsed time; starting time

% Please include ReadDNSData3, Particle_Trake_Intermediate,
% Particle_Track_Core_RK2_032923_mex together with this script in the same
% folder. Please note that the directories need to adjust accordingly to
% the location of the data and mex function

%% Input Values

clearvars -except x_trajectory y_trajectory z_trajectory
% Use approprate directory where the datafiles are saved
directory = 'F:\Bubbles_DNS\2_05'; % Path can be changed

timeSkip = 1; % default 1 for consecutive data files (2nd timestep - 1st timestep)
delta_t = 0.001*timeSkip;
n = 100000; % NUMBER OF PARTICLES
% n = numel(x_trajectory(:,1)); %keeping particle number consistent
n = single(n);
StartTime = 550; %Starting Timestep
EndTime = 3000; % End TIMESTEP + 1 for RK2
number_of_timesteps = (EndTime - StartTime)/timeSkip + 1; % including the first and last timesteop

% Fluid Properties
vis_kin = 1; %mm2/s

%% Importing VOF timestep 1
% NECESSARY STEP TO ENSURE PARTICLES ARE NOT SEEDED IN BUBBLE LOCATIONS
% to make sure coordinates inputted don't overlap with the bubbles

a = StartTime;
if a <= 9
    myfilename = sprintf('Matlab_stag00%d.dat', a);
    
elseif a > 9 && a <= 99
    myfilename = sprintf('Matlab_stag0%d.dat', a);
else
    myfilename = sprintf('Matlab_stag%d.dat', a);
end
filename = fullfile(directory, myfilename); % LOCATION OF THE DATA IS CHANGED HERE
[Data] = ReadDNSData4(filename);

vof_1(:,:,:) = Data.vof; %logical
u_1(:,:,:) = Data.u; %m/s
v_1(:,:,:) = Data.v; %m/s
w_1(:,:,:) = Data.w; %m/s
p_1(:,:,:) = Data.p; %Pascal

u_1 = u_1*1000; %mm/s
v_1 = v_1*1000; %mm/s
w_1 = w_1*1000; %mm/s

Nx = numel(u_1(1,:,1));
Ny = numel(u_1(:,1,1));
Nz = numel(u_1(:,1,1));
%% Creating the GRID and Seeding Particles

x = double((1:1:Nx));
y = double((1:1:Ny));
z = double((1:1:Nz));
x = transpose(x);
y = transpose(y);
z = transpose(z);
x_range = [min(x), max(x)];
y_range = [min(y), max(y)];
z_range = [min(z), max(z)];
rows_outofDomain_yz = zeros([n 1]); %rows that correspond to particles that have reached the domain in y and z directions

% Initialize an array to store the generated coordinates
coordinates = zeros(n, 3);

% Generate random 3D coordinates
for i = 1:n
    x_coordinate = randi(x_range);
    y_coordinate = randi(y_range);
    z_coordinate = randi(z_range);
    coordinates(i, :) = [x_coordinate, y_coordinate, z_coordinate];
end
positions = coordinates;
x0 = positions(:,1);
y0 = positions(:,2);
z0 = positions(:,3);

parfor i = 1:numel(x0)
    while vof_1(y0(i),x0(i),z0(i)) == 1
        disp(i)
        x0(i) = randi(x_range);
        y0(i) = randi(y_range);
        z0(i) = randi(z_range);
    end
end

x0 = x0*(40/165); %from indices to real world coordinates in mm
y0 = y0*(40/165);
z0 = z0*(40/165);

x = x*(40/165); %from indices to real world coordinates in mm
y = y*(40/165);
z = z*(40/165);

%%%%%%% You can load custom x0 y0 z0 data over here but make sure they do
%%%%%%% NOT overlap with bubbles from the first timestep simulated
x0 = x_trajectory(:,end); %starting positions for the current set should be the same as the ending positions for the last set
y0 = y_trajectory(:,end);
z0 = z_trajectory(:,end);

x0_original = x0;
y0_original = y0;
z0_original = z0;

clearvars x_trajectory y_trajectory z_trajectory
%% Eulerian Velocity Gradient Invariant Conditioned upon Lagrangian Particle Position

method = 2; %% Specify the method of Differencing vlocity over space
% method = 1; Forward Finite Differencing
% method = 2; 2nd order Central Finite Differencing
% method = 3; 4th order central finite differencing
dx = x(2) - x(1); %mm
dy = dx;
dz = dx;

k = u_1; % eulerian velocity
[dkdx, dkdy, dkdz] = Finite_diff_3Ddz(k, method, dx, dy, dz);
dudx = dkdx; %3d matrix
dudy = dkdy;
dudz = dkdz;
clearvars dkdx dkdy dkdz k
k = v_1;
[dkdx, dkdy, dkdz] = Finite_diff_3Ddz(k, method, dx, dy, dz);
dvdx = dkdx;
dvdy = dkdy;
dvdz = dkdz;
clearvars dkdx dkdy dkdz k
k = w_1;
[dkdx, dkdy, dkdz] = Finite_diff_3Ddz(k, method, dx, dy, dz);
dwdx = dkdx;
dwdy = dkdy;
dwdz = dkdz;
clearvars dkdx dkdy dkdz k

% dudx filter bubbles
dudx = dudx(:);
dudy = dudy(:);
dudz = dudz(:);
dvdx = dvdx(:);
dvdy = dvdy(:);
dvdz = dvdz(:);
dwdx = dwdx(:);
dwdy = dwdy(:);
dwdz = dwdz(:);
vof = vof_1(:);
for i = 1:numel(vof)
    if vof(i) == 1
        dudx(i) = nan;
        dudy(i) = nan;
        dudz(i) = nan;
        dvdx(i) = nan;
        dvdy(i) = nan;
        dvdz(i) = nan;
        dwdx(i) = nan;
        dwdy(i) = nan;
        dwdz(i) = nan;
    else
    end
end
dudx = reshape(dudx,Ny,Nx,Nz);
dudy = reshape(dudy,Ny,Nx,Nz);
dudz = reshape(dudz,Ny,Nx,Nz);
dvdx = reshape(dvdx,Ny,Nx,Nz);
dvdy = reshape(dvdy,Ny,Nx,Nz);
dvdz = reshape(dvdz,Ny,Nx,Nz);
dwdx = reshape(dwdx,Ny,Nx,Nz);
dwdy = reshape(dwdy,Ny,Nx,Nz);
dwdz = reshape(dwdz,Ny,Nx,Nz);

cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\GetVelocityGradientEulcondLag')
[dudx_eL, dudy_eL, dudz_eL, dvdx_eL, dvdy_eL, dvdz_eL, dwdx_eL, dwdy_eL, dwdz_eL] = GetVelocityGradientEulcondLag_mex(x0, y0, z0, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);

dudx_eL_first = dudx_eL; % should contain nan values at bubble locations
dudy_eL_first = dudy_eL;
dudz_eL_first = dudz_eL;
dvdx_eL_first = dvdx_eL;
dvdy_eL_first = dvdy_eL;
dvdz_eL_first = dvdz_eL;
dwdx_eL_first = dwdx_eL;
dwdy_eL_first = dwdy_eL;
dwdz_eL_first = dwdz_eL;

%% Eulerian Energy Dissipation
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions');

e = nan([1 Nx*Ny*Nz]);
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\EnergyDissipation')
[e] = EnergyDissipation_mex(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vis_kin,e); %free of bubbles
e_first = e; % average eulerian energy dissipation for the first timestep devoid of bubbles
clearvars dudx_eL dudy_eL dudz_eL dvdx_eL dvdy_eL dvdz_eL dwdx_eL dwdy_eL dwdz_eL dudx_eL dudy_eL dudz_eL dvdx_eL dvdy_eL dvdz_eL dwdx_eL dwdy_eL dwdz_eL e vof

%% Lagrangian Pressure Computation

cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\GetPressure');
p1 = p_1;
[P] = GetPressure_mex(p1,x0,y0,z0,x,y,z); %collecting pressure values at points free of bubbles
P1 = P; % Lagrangian pressure of first timestep

clearvars p_1 p1 P
%% Eulerian U_RMS first timestep

cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\U_RMS1');
[u_rms, v_rms, w_rms] = U_RMS1_mex(u_1,v_1,w_1,vof_1); 
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions');

u_rms_first = u_rms;
v_rms_first = v_rms;
w_rms_first = w_rms;

clearvars u_rms v_rms w_rms u_1 v_1 w_1 vof_1
%% Reading Data

q = StartTime; %StartTime

x_track = nan(n,1); %n = number of particles; row
y_track = nan(n,1);
z_track = nan(n,1);

u_track = nan(n,1); %n = number of particles; row
v_track = nan(n,1);
w_track = nan(n,1);

%GLOBAL INITIALISATION 

x_trajectory = nan(n,1); %n = number of particles; row
y_trajectory = nan(n,1);
z_trajectory = nan(n,1);

u_trajectory = nan(n,1);
v_trajectory = nan(n,1);
w_trajectory = nan(n,1);

u1 = zeros(165,330,165,2);
v1 = zeros(165,330,165,2);
w1 = zeros(165,330,165,2);
vof = zeros(165,330,165,2);
p1 = zeros(165,330,165); % just reading pressure of the second timestep of the two

temp_u1 = zeros(size(u1));
temp_v1 = zeros(size(temp_u1));
temp_w1 = zeros(size(temp_u1));
temp_vof = zeros(size(temp_u1));

u_last = nan(n,1);
v_last = nan(n,1);
w_last = nan(n,1);

%eulerian u_rms
u_rms1 = nan(1,number_of_timesteps-1); %parse with values of the second timestep
v_rms1 = nan(1,number_of_timesteps-1);
w_rms1 = nan(1,number_of_timesteps-1);

P_track = nan(n,number_of_timesteps-1);

dudx_eL_second = nan(n,number_of_timesteps-1); %parse with values of the second timestep
dudy_eL_second = nan(n,number_of_timesteps-1);
dudz_eL_second = nan(n,number_of_timesteps-1);
dvdx_eL_second = nan(n,number_of_timesteps-1);
dvdy_eL_second = nan(n,number_of_timesteps-1);
dvdz_eL_second = nan(n,number_of_timesteps-1);
dwdx_eL_second = nan(n,number_of_timesteps-1);
dwdy_eL_second = nan(n,number_of_timesteps-1);
dwdz_eL_second = nan(n,number_of_timesteps-1);
e_second = nan(1,number_of_timesteps-1); %average energy dissipation rate

% Main Loop
for b = 1:number_of_timesteps

if b == number_of_timesteps
        break
else
r = q + timeSkip;

for a=q:timeSkip:r

    if a <= 9
        myfilename1 = sprintf('Matlab_stag00%d.dat', a);
        
    elseif a > 9 && a <= 99
        myfilename1 = sprintf('Matlab_stag0%d.dat', a);
    else
        myfilename1 = sprintf('Matlab_stag%d.dat', a);
        
    end
    filename1 = fullfile(directory, myfilename1); % LOCATION OF THE DATA IS CHANGED HERE


[Data] = ReadDNSData4(filename1);

if a == q
    temp_u1(:,:,:,1) = Data.u;
    temp_v1(:,:,:,1) = Data.v;
    temp_w1(:,:,:,1) = Data.w;
    temp_vof(:,:,:,1) = Data.vof;
else %second timestep
    %if a == r
    temp_u1(:,:,:,2) = Data.u;
    temp_v1(:,:,:,2) = Data.v;
    temp_w1(:,:,:,2) = Data.w;
    temp_vof(:,:,:,2) = Data.vof;
    p1(:,:,:) = Data.p; %second timestep
end
end

u1(:,:,:,1) = temp_u1(:,:,:,1);
u1(:,:,:,2) = temp_u1(:,:,:,2);
u1 = u1*1000; %mm/s

v1(:,:,:,1) = temp_v1(:,:,:,1);
v1(:,:,:,2) = temp_v1(:,:,:,2);
v1 = v1*1000; %mm/s

w1(:,:,:,1) = temp_w1(:,:,:,1);
w1(:,:,:,2) = temp_w1(:,:,:,2);
w1 = w1*1000; %mm/s

vof(:,:,:,1) = temp_vof(:,:,:,1);
vof(:,:,:,2) = temp_vof(:,:,:,2);

%velocity data is given in m/s
%convert from m/s to mm/s

clear vars temp_u1 temp_v1 temp_w1 temp_vof

% Eulerian U_RMS calculation per timestep [EULERIAN]

u_1 = u1(:,:,:,2);
v_1 = v1(:,:,:,2);
w_1 = w1(:,:,:,2);
vof_1 = vof(:,:,:,2);

cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\U_RMS1');
[u_rms, v_rms, w_rms] = U_RMS1_mex(u_1,v_1,w_1,vof_1); 
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions');

u_rms1(1,b) = u_rms;
v_rms1(1,b) = v_rms;
w_rms1(1,b) = w_rms;

clear vars u_rms v_rms w_rms

nth_time = b;

% PARTICLE TRACKING ALGORITHM
[x_track_int, y_track_int, z_track_int, u_track_int, v_track_int, w_track_int,rows_outofDomain_yz, u_last, v_last, w_last] = Particle_Track_Intermediate_PBC(x,y,z,u1,v1,w1,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof,rows_outofDomain_yz,nth_time,number_of_timesteps,u_last,v_last,w_last);

q = q + 1;

x0 = x_track_int;
y0 = y_track_int;
z0 = z_track_int;

x_track(:,b) = x_track_int;
y_track(:,b) = y_track_int;
z_track(:,b) = z_track_int;
u_track(:,b) = u_track_int;
v_track(:,b) = v_track_int;
w_track(:,b) = w_track_int;

% Lagrangian Pressure Computation
% Second timestep
% x0, y0, z0 are updated
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\GetPressure');
[P] = GetPressure_mex(p1,x0,y0,z0,x,y,z); %Get pressure at the next particle position
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions');
P_track(:,b) = P;

% Eulerian Velocity Gradient Invariant conditioned on Lagrangian Ref Frame

% method for differencing is consistent
k = u_1;
[dkdx, dkdy, dkdz] = Finite_diff_3Ddz(k, method, dx, dy, dz);
dudx = dkdx;
dudy = dkdy;
dudz = dkdz;
clearvars dkdx dkdy dkdz k
k = v_1;
[dkdx, dkdy, dkdz] = Finite_diff_3Ddz(k, method, dx, dy, dz);
dvdx = dkdx;
dvdy = dkdy;
dvdz = dkdz;
clearvars dkdx dkdy dkdz k
k = w_1;
[dkdx, dkdy, dkdz] = Finite_diff_3Ddz(k, method, dx, dy, dz);
dwdx = dkdx;
dwdy = dkdy;
dwdz = dkdz;
clearvars dkdx dkdy dkdz k

% dudx filter bubbles
dudx = dudx(:);
dudy = dudy(:);
dudz = dudz(:);
dvdx = dvdx(:);
dvdy = dvdy(:);
dvdz = dvdz(:);
dwdx = dwdx(:);
dwdy = dwdy(:);
dwdz = dwdz(:);
vof = vof_1(:);
for i = 1:numel(vof)
    if vof(i) == 1
        dudx(i) = nan;
        dudy(i) = nan;
        dudz(i) = nan;
        dvdx(i) = nan;
        dvdy(i) = nan;
        dvdz(i) = nan;
        dwdx(i) = nan;
        dwdy(i) = nan;
        dwdz(i) = nan;
    else
    end
end
dudx = reshape(dudx,Ny,Nx,Nz);
dudy = reshape(dudy,Ny,Nx,Nz);
dudz = reshape(dudz,Ny,Nx,Nz);
dvdx = reshape(dvdx,Ny,Nx,Nz);
dvdy = reshape(dvdy,Ny,Nx,Nz);
dvdz = reshape(dvdz,Ny,Nx,Nz);
dwdx = reshape(dwdx,Ny,Nx,Nz);
dwdy = reshape(dwdy,Ny,Nx,Nz);
dwdz = reshape(dwdz,Ny,Nx,Nz);

clearvars u_1 v_1 w_1 vof_1 vof
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\GetVelocityGradientEulcondLag');
[dudx_eL, dudy_eL, dudz_eL, dvdx_eL, dvdy_eL, dvdz_eL, dwdx_eL, dwdy_eL, dwdz_eL] = GetVelocityGradientEulcondLag_mex(x0, y0, z0, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);
clear vof

dudx_eL_second(:,b) = dudx_eL;
dudy_eL_second(:,b) = dudy_eL;
dudz_eL_second(:,b) = dudz_eL;
dvdx_eL_second(:,b) = dvdx_eL;
dvdy_eL_second(:,b) = dvdy_eL;
dvdz_eL_second(:,b) = dvdz_eL;
dwdx_eL_second(:,b) = dwdx_eL;
dwdy_eL_second(:,b) = dwdy_eL;
dwdz_eL_second(:,b) = dwdz_eL;

% Eulerian Energy Dissipation
e = nan([1 Nx*Ny*Nz]);
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\EnergyDissipation')
[e] = EnergyDissipation_mex(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vis_kin,e); %free of bubbles
e_second(b) = e; % average eulerian energy dissipation for the second timestep devoid of bubbles
clearvars dudx dudy_eL dudz_eL dvdx_eL dvdy_eL dvdz_eL dwdx_eL dwdy_eL dwdz_eL dudx_eL dudy_eL dudz_eL dvdx_eL dvdy_eL dvdz_eL dwdx_eL dwdy_eL dwdz_eL
cd('C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions');

%
disp(['Timestep(s) Completed: ' num2str(b)]);

% clear vars u1 v1 w1 vof
end % repitition ends
end

x_trajectory = [x0_original x_track];
y_trajectory = [y0_original y_track];
z_trajectory = [z0_original z_track];

u_trajectory = u_track;
v_trajectory = v_track;
w_trajectory = w_track;

u_trajectory = [u_trajectory u_last];
v_trajectory = [v_trajectory v_last];
w_trajectory = [w_trajectory w_last];

u_rms = [u_rms_first u_rms1];
v_rms = [v_rms_first v_rms1];
w_rms = [w_rms_first w_rms1];

Pressure = [P1 P_track];

dudx_eL = [dudx_eL_first dudx_eL_second];
dudy_eL = [dudy_eL_first dudy_eL_second];
dudz_eL = [dudz_eL_first dudz_eL_second];
dvdx_eL = [dvdx_eL_first dvdx_eL_second];
dvdy_eL = [dvdy_eL_first dvdy_eL_second];
dvdz_eL = [dvdz_eL_first dvdz_eL_second];
dwdx_eL = [dwdx_eL_first dwdx_eL_second];
dwdy_eL = [dwdy_eL_first dwdy_eL_second];
dwdz_eL = [dwdz_eL_first dwdz_eL_second];

clear vars x_track y_track z_track u_track v_track w_track x_track_int y_track_int z_track_int u_track_int v_track_int w_track_int Data dudx dudy dudz dvdx dvdy dvdz dwdx dwdy dwdz dudx_eL_first dudx_eL_second dudy_eL_first dudy_eL_second dudz_eL_first dudz_eL_second dvdx_eL_first dvdx_eL_second dvdy_eL_first dvdy_eL_second dvdz_eL_first dvdz_eL_second dwdx_eL_first dwdx_eL_second dwdy_eL_first dwdy_eL_second dwdz_eL_first dwdz_eL_second
clear vars P_track P1 u_last v_last w_last u_rms1 v_rms1 w_rms1 u_rms_first v_rms_first w_rms_first x_coordinate y_position z_position

dUdX_eL = cat(3,dudx_eL,dudy_eL,dudz_eL,dvdx_eL,dvdy_eL,dvdz_eL,dwdx_eL,dwdy_eL,dwdz_eL);
Trajectory_Lagrangian = cat(3, x_trajectory, y_trajectory, z_trajectory);
Velocity_Lagrangian = cat(3, u_trajectory, v_trajectory, w_trajectory);
EnergyDissipationRate_Eulerian = [e_first e_second]; %average energy dissipation rate
%% Plotting Particle Trajectories in 3D
% 
% 
x_track = x_trajectory;
y_track = y_trajectory;
z_track = z_trajectory;

%trajectory = [x_track; y_track; z_track];

F1 = figure;
for i=1:n
    plot3(x_track(i,:),y_track(i,:),z_track(i,:))
    %hold on
    %plot3(x_track_GP(i,:),y_track_GP(i,:),z_track_GP(i,:))
    %hold on
    %plot3(x_track_RK2(i,:),y_track_RK2(i,:),z_track_RK2(i,:))
    axis square
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    %L=legend('particle track','GetPosition','particle track with RK2','location','northeast')

end
hold on 
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
hold off

toc
