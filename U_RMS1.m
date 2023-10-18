%%%%%% This code calculates the U_RMS V_RMS W_RMS values of the nodal
%%%%%% points that are devoid of bubble locations
% 

function [u_rms, v_rms, w_rms] = U_RMS1(u_1,v_1,w_1,vof_1)

u1 = u_1;
v1 = v_1;
w1 = w_1;
vof = vof_1;

nx = numel(u1(1,:,1));
ny = numel(u1(:,1,1));
nz = numel(u1(1,1,:));

u = reshape(u1, [nx*ny*nz, 1]); %mm/s
v = reshape(v1, [nx*ny*nz, 1]); %mm/s
w = reshape(w1, [nx*ny*nz, 1]); %mm/s
vof = reshape(vof, [nx*ny*nz,1]);

% find rows that correspond to bubble locations (logical value of 1)

rows_bubbles = find(vof == 1);
rows_withoutBubbles_idx = true(size(u));
rows_withoutBubbles_idx(rows_bubbles) = false;
u_withoutBubbles = u(rows_withoutBubbles_idx);
v_withoutBubbles = v(rows_withoutBubbles_idx);
w_withoutBubbles = w(rows_withoutBubbles_idx);

% Following calculations are devoid of Bubbles

u_prime = u_withoutBubbles - mean(u_withoutBubbles);
v_prime = v_withoutBubbles - mean(v_withoutBubbles);
w_prime = w_withoutBubbles - mean(w_withoutBubbles);

u_rms = std(u_prime);
v_rms = std(v_prime);
w_rms = std(w_prime);

end