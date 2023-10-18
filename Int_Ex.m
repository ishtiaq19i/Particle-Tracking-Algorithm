% 
% % %INPUT
% % % 
% % L = 1;
% % N = 2;
% % C = 3;
% % M = 4;
% % S = 5;
% % 
% % method_interp = S; %Linear
% 
x = double([0.24242425;0.48484850;0.72727275;0.96969700;1.2121212;1.4545455;1.6969697;1.9393940;2.1818182;2.4242425;2.6666667;2.9090910;3.1515152;3.3939395;3.6363637;3.8787880;4.1212120;4.3636365;4.6060605;4.8484850;5.0909090;5.3333335;5.5757575;5.8181820;6.0606060;6.3030305;6.5454545;6.7878790;7.0303030;7.2727275;7.5151515;7.7575760;8;8.2424240;8.4848480;8.7272730;8.9696970;9.2121210;9.4545450;9.6969700;9.9393940;10.181818;10.424242;10.666667;10.909091;11.151515;11.393939;11.636364;11.878788;12.121212;12.363636;12.606061;12.848485;13.090909;13.333333;13.575758;13.818182;14.060606;14.303030;14.545455;14.787879;15.030303;15.272727;15.515152;15.757576;16;16.242424;16.484848;16.727272;16.969696;17.212122;17.454546;17.696970;17.939394;18.181818;18.424242;18.666666;18.909090;19.151516;19.393940;19.636364;19.878788;20.121212;20.363636;20.606060;20.848484;21.090910;21.333334;21.575758;21.818182;22.060606;22.303030;22.545454;22.787878;23.030304;23.272728;23.515152;23.757576;24;24.242424;24.484848;24.727272;24.969696;25.212122;25.454546;25.696970;25.939394;26.181818;26.424242;26.666666;26.909090;27.151516;27.393940;27.636364;27.878788;28.121212;28.363636;28.606060;28.848484;29.090910;29.333334;29.575758;29.818182;30.060606;30.303030;30.545454;30.787878;31.030304;31.272728;31.515152;31.757576;32;32.242424;32.484848;32.727272;32.969696;33.212120;33.454544;33.696968;33.939392;34.181820;34.424244;34.666668;34.909092;35.151516;35.393940;35.636364;35.878788;36.121212;36.363636;36.606060;36.848484;37.090908;37.333332;37.575756;37.818180;38.060608;38.303032;38.545456;38.787880;39.030304;39.272728;39.515152;39.757576;40;40.242424;40.484848;40.727272;40.969696;41.212120;41.454544;41.696968;41.939392;42.181820;42.424244;42.666668;42.909092;43.151516;43.393940;43.636364;43.878788;44.121212;44.363636;44.606060;44.848484;45.090908;45.333332;45.575756;45.818180;46.060608;46.303032;46.545456;46.787880;47.030304;47.272728;47.515152;47.757576;48;48.242424;48.484848;48.727272;48.969696;49.212120;49.454544;49.696968;49.939392;50.181820;50.424244;50.666668;50.909092;51.151516;51.393940;51.636364;51.878788;52.121212;52.363636;52.606060;52.848484;53.090908;53.333332;53.575756;53.818180;54.060608;54.303032;54.545456;54.787880;55.030304;55.272728;55.515152;55.757576;56;56.242424;56.484848;56.727272;56.969696;57.212120;57.454544;57.696968;57.939392;58.181820;58.424244;58.666668;58.909092;59.151516;59.393940;59.636364;59.878788;60.121212;60.363636;60.606060;60.848484;61.090908;61.333332;61.575756;61.818180;62.060608;62.303032;62.545456;62.787880;63.030304;63.272728;63.515152;63.757576;64;64.242424;64.484848;64.727272;64.969696;65.212120;65.454544;65.696968;65.939392;66.181816;66.424240;66.666664;66.909088;67.151512;67.393936;67.636360;67.878784;68.121216;68.363640;68.606064;68.848488;69.090912;69.333336;69.575760;69.818184;70.060608;70.303032;70.545456;70.787880;71.030304;71.272728;71.515152;71.757576;72;72.242424;72.484848;72.727272;72.969696;73.212120;73.454544;73.696968;73.939392;74.181816;74.424240;74.666664;74.909088;75.151512;75.393936;75.636360;75.878784;76.121216;76.363640;76.606064;76.848488;77.090912;77.333336;77.575760;77.818184;78.060608;78.303032;78.545456;78.787880;79.030304;79.272728;79.515152;79.757576;80]);
y = double([0.24242425;0.48484850;0.72727275;0.96969700;1.2121212;1.4545455;1.6969697;1.9393940;2.1818182;2.4242425;2.6666667;2.9090910;3.1515152;3.3939395;3.6363637;3.8787880;4.1212120;4.3636365;4.6060605;4.8484850;5.0909090;5.3333335;5.5757575;5.8181820;6.0606060;6.3030305;6.5454545;6.7878790;7.0303030;7.2727275;7.5151515;7.7575760;8;8.2424240;8.4848480;8.7272730;8.9696970;9.2121210;9.4545450;9.6969700;9.9393940;10.181818;10.424242;10.666667;10.909091;11.151515;11.393939;11.636364;11.878788;12.121212;12.363636;12.606061;12.848485;13.090909;13.333333;13.575758;13.818182;14.060606;14.303030;14.545455;14.787879;15.030303;15.272727;15.515152;15.757576;16;16.242424;16.484848;16.727272;16.969696;17.212122;17.454546;17.696970;17.939394;18.181818;18.424242;18.666666;18.909090;19.151516;19.393940;19.636364;19.878788;20.121212;20.363636;20.606060;20.848484;21.090910;21.333334;21.575758;21.818182;22.060606;22.303030;22.545454;22.787878;23.030304;23.272728;23.515152;23.757576;24;24.242424;24.484848;24.727272;24.969696;25.212122;25.454546;25.696970;25.939394;26.181818;26.424242;26.666666;26.909090;27.151516;27.393940;27.636364;27.878788;28.121212;28.363636;28.606060;28.848484;29.090910;29.333334;29.575758;29.818182;30.060606;30.303030;30.545454;30.787878;31.030304;31.272728;31.515152;31.757576;32;32.242424;32.484848;32.727272;32.969696;33.212120;33.454544;33.696968;33.939392;34.181820;34.424244;34.666668;34.909092;35.151516;35.393940;35.636364;35.878788;36.121212;36.363636;36.606060;36.848484;37.090908;37.333332;37.575756;37.818180;38.060608;38.303032;38.545456;38.787880;39.030304;39.272728;39.515152;39.757576;40]);
z = y;

nth_time = 50;
u_last = nan(20,1);
v_last = u_last;
w_last = u_last;
number_of_timesteps = 5;

% Change to a specific directory
newDirectory = 'C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE';
cd(newDirectory);

% Check the current working directory
currentDirectory = pwd;
disp(['Current working directory: ' currentDirectory])

u1 = load('u1.mat');
v1 = load('v1.mat');
w1 = load('w1.mat');
vof = load('vof.mat');

% Change to a specific directory
newDirectory = 'C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions';
cd(newDirectory);

% Check the current working directory
currentDirectory = pwd;
disp(['Current working directory: ' currentDirectory])

x0 = double([71.0303030303030;56.9696969696970;13.5757575757576;13.3333333333333;53.3333333333333;0.727272727272727;23.2727272727273;55.5151515151515;31.5151515151515;9.21212121212121;18.9090909090909;63.7575757575758;42.4242424242424;4.12121212121212;5.33333333333333;64.4848484848485;25.6969696969697;43.8787878787879;28.6060606060606;47.0303030303030]);
y0 = double([15.0303030303030;10.9090909090909;6.06060606060606;20.1212121212121;20.3636363636364;31.7575757575758;21.0909090909091;10.4242424242424;16.4848484848485;12.3636363636364;35.3939393939394;8;8;1.21212121212121;18.1818181818182;8.72727272727273;1.93939393939394;36.3636363636364;0.969696969696970;29.5757575757576]);
z0 = double([35.6363636363636;12.3636363636364;39.7575757575758;6.30303030303030;8.24242424242424;37.3333333333333;20.8484848484848;5.09090909090909;14.7878787878788;10.9090909090909;32.4848484848485;23.0303030303030;16.9696969696970;31.2727272727273;32;29.5757575757576;12.8484848484849;16;15.0303030303030;16.7272727272727]);

n = 20;
n = single(n);
delta_t = 0.001;

x_trajectory = nan(n,1); %n = number of particles; row
y_trajectory = nan(n,1);
z_trajectory = nan(n,1);

u_trajectory = nan(n,1);
v_trajectory = nan(n,1);
w_trajectory = nan(n,1);

rows_outofDomain_yz = zeros([n 1]);

% function[x_track_int, y_track_int, z_track_int, u_track_int, v_track_int, w_track_int,rows_outofDomain_yz, u_last, v_last, w_last] = Particle_Track_Intermediate_PBC(x,y,z,u1,v1,w1,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof,rows_outofDomain_yz,nth_time,number_of_timesteps,u_last,v_last,w_last)

u1_t1 = u1.u1(:,:,:,1); 
v1_t1 = v1.v1(:,:,:,1);
w1_t1 = w1.w1(:,:,:,1);

u1_t2 = u1.u1(:,:,:,2);
v1_t2 = v1.v1(:,:,:,2);
w1_t2 = w1.w1(:,:,:,2);

vof_t1 = vof.vof(:,:,:,1);
vof_t2 = vof.vof(:,:,:,2);

% u1_t1 = u1(:,:,:,1); 
% v1_t1 = v1(:,:,:,1);
% w1_t1 = w1(:,:,:,1);
% 
% u1_t2 = u1(:,:,:,2);
% v1_t2 = v1(:,:,:,2);
% w1_t2 = w1(:,:,:,2);
% 
% vof_t2 = vof(:,:,:,2);
% 
% % Change to a specific directory
% newDirectory = 'C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions\codegen\mex\Particle_Track_Core_RK2_PBC';
% cd(newDirectory);
% 
% % Check the current working directory
% currentDirectory = pwd;
% disp(['Current working directory: ' currentDirectory])
        
tic

[x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,rows_outofDomain_yz, u_last, v_last, w_last] = Particle_Track_Core_RK2_PBC(x,y,z,u1_t1,u1_t2,v1_t1,v1_t2,w1_t1,w1_t2,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof_t2,rows_outofDomain_yz,nth_time,number_of_timesteps,u_last,v_last,w_last,vof_t2);
% [x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,rows_outofDomain_yz, u_last, v_last, w_last] = Particle_Track_Core_RK2_PBC_mex(x,y,z,u1_t1,u1_t2,v1_t1,v1_t2,w1_t1,w1_t2,x0,y0,z0,n,delta_t,x_trajectory,y_trajectory,z_trajectory,u_trajectory,v_trajectory,w_trajectory,vof_t2,rows_outofDomain_yz,nth_time,number_of_timesteps,u_last,v_last,w_last);

toc

Lx = max(x);
Ly = max(y);
Lz = max(z);

parfor i=1:numel(x0)
    
    % Periodic Boundary Conditions

    if x_trajectory(i) > max(x)
        x_trajectory(i) = x_trajectory(i) - (Lx-min(x));
%         x_trajectory(i) = nan;
%         u_trajectory(i) = nan;
    elseif x_trajectory(i) < min(x)
        x_trajectory(i) = x_trajectory(i) + (Lx-min(x));
%         x_trajectory(i) = nan;
%         u_trajectory(i) = nan;
    end

    if y_trajectory(i) > max(y)
        y_trajectory(i) = y_trajectory(i) - (Ly-min(y));
        rows_outofDomain_yz(i) = 1;
%         y_trajectory(i) = nan;
%         v_trajectory(i) = nan;
    elseif y_trajectory(i) < min(y)
        y_trajectory(i) = y_trajectory(i) + (Ly-min(y));
        rows_outofDomain_yz(i) = 1;
%         y_trajectory(i) = nan;
%         v_trajectory(i) = nan;
    end

    if z_trajectory(i) > max(z)
        z_trajectory(i) = z_trajectory(i) - (Lz-min(z));
        rows_outofDomain_yz(i) = 1;
%         z_trajectory(i) = nan;
%         w_trajectory(i) = nan;
    elseif z_trajectory(i) < min(z)
        z_trajectory(i) = z_trajectory(i) + (Lz-min(z));
        rows_outofDomain_yz(i) = 1;
%         z_trajectory(i) = nan;
%         w_trajectory(i) = nan;
    end
end

x_trajectory_nodal = round((165/40).*x_trajectory); %in pixel coordinates
y_trajectory_nodal = round((165/40).*y_trajectory);
z_trajectory_nodal = round((165/40).*z_trajectory);

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

parfor i=1:numel(x0)

    if isnan(y_trajectory_nodal(i)) || isnan(x_trajectory_nodal(i)) || isnan(z_trajectory_nodal(i))
        x_trajectory(i) = nan;
        y_trajectory(i) = nan;
        z_trajectory(i) = nan;
    else
        if vof_t2(y_trajectory_nodal(i),x_trajectory_nodal(i),z_trajectory_nodal(i)) == 1 % Bubble Termination

        x_trajectory(i) = nan;
        y_trajectory(i) = nan;
        z_trajectory(i) = nan;

        u_trajectory(i) = nan;
        v_trajectory(i) = nan;
        w_trajectory(i) = nan;

        disp(['Bubble-terminated Particle(s) after RK2: ' num2str(i)]);
        else
        end
    end
end

x_track_int = x_trajectory;
y_track_int = y_trajectory;
z_track_int = z_trajectory;
u_track_int = u_trajectory;
v_track_int = v_trajectory;
w_track_int = w_trajectory;
% 
% % Change to a specific directory
% newDirectory = 'C:\Users\iislam30\Documents\Bubble_DNS_Bruno\CODE\Periodic Boundary Conditions';
% cd(newDirectory);
% 
% %Check the current working directory
% currentDirectory = pwd;
% disp(['Current working directory: ' currentDirectory]);
% 
% end
