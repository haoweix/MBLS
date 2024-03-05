run ./looping_star_k_traj_fov24.m
% load('traj_slew200_grad1.mat');
load('traj_slew200_grad1_230521.mat');

% scale grad to g/cm
gambar = 4.258; % 4.258 kHz/g
fov = 24; % cm
nspokes = 31; % # of RF pulses and therefore, echos
xres = 64 ;
% grad_scale = grad * xres / fov / 2 / pi;
grad_scale = grad;

% check k-space traj
kx = squeeze(gambar .* cumsum(grad(1,:,:),3)*dt) * fov;
ky = squeeze(gambar .* cumsum(grad(2,:,:),3)*dt) * fov;
kz = squeeze(gambar .* cumsum(grad(3,:,:),3)*dt) * fov;

figure;
plot3(kx',ky',kz','LineWidth',2);

xlabel('k_x')
ylabel('k_y')
zlabel('k_z')
title('max grad 1 k-space trajectory')


%% build ramp up and down for FID module
rfpre = 1;
trf = 3;
tmax = nspokes*tseg; %might only need 64, but for completeness
t = [-rfpre:((nspokes-1)*nseg + trf)]*dt;  %start at RF pulse and end just after last RF pulse

g = cat(3, grad(:,:,end-1:end),grad(:,:,1:((nspokes-1)*nseg + trf)));

gamp_max_plane = max(g(:,:,1),[],'all');

slewramp = 0.5; % g/cm/ms - slew for ramp-up/down to looping.
dt = 0.004; % 4 us
nplanes = size(grad,2);

ramplen = ceil(gamp_max_plane/slewramp/dt); % universal ramp length for every plane
rampup = zeros(3, nplanes, ramplen);
rampdown = zeros(3, nplanes, ramplen);
for ni = 1:3
    for nj = 1:nplanes
        rampup(ni,nj,:) = g(ni,nj,1) .* [0:(ramplen-1)]/ramplen;
        rampdown(ni,nj,:) = g(ni,nj,end) .* [(ramplen-1):-1:0]/ramplen;
    end
end
g1 = cat(3,rampup, g, rampdown);

% check gradient waveform
figure;
subplot(311);
plot(squeeze(g1(1,:,:)).');
xlabel('time');
ylabel('gx with max 1');
subplot(312);
plot(squeeze(g1(2,:,:)).');
xlabel('time');
ylabel('gy with max 1');
subplot(313);
plot(squeeze(g1(3,:,:)).');
xlabel('time');
ylabel('gz with max 1');
%% build RF pulse

rf = zeros(size(g1,3),1);
rflen = 3;
for lp = 1:nspokes
    stloc = ramplen+1+(lp-1)*nseg;
    rf(stloc:(stloc+trf-1)) = repmat(rfamp,rflen,1);
end

%% plot k-traj fid
ctr = ramplen + rfpre + 1;
k1 = gambar.*(cumsum(g1(ctr:end))*dt)*fov; % unit are not cm-1, but 1/FOV or sample number
karray_fid = zeros(3,length(g1),nspokes);
for ispoke = 1:nspokes
    st_pt = ctr + (ispoke - 1) * nseg;
    karray_fid(:,st_pt:end,ispoke) =  gambar.*(cumsum(g1(:,1,st_pt:end))*dt)*fov; % unit are not cm-1, but 1/FOV or sample number
end

figure;
plot3(squeeze(karray_fid(1,:,1)),squeeze(karray_fid(2,:,1)),squeeze(karray_fid(3,:,1)),'.');

xlabel('k_x')
ylabel('k_y')
zlabel('k_z')
title('max grad 1 k-space trajectory')

%% build the GRE module
rfpre = 1; % samples, 8 us - amount before "center of RF"
ctr = ramplen + rfpre + 1;

g = cat(3, grad(:,:,end-nseg+1:end),grad(:,:,1:nspokes*nseg));

rampup = zeros(3, nplanes, ramplen);
rampdown = zeros(3, nplanes, ramplen);
for ni = 1:3
    for nj = 1:nplanes
        rampup(ni,nj,:) = g(ni,nj,1) .* [0:(ramplen-1)]/ramplen;
        rampdown(ni,nj,:) = g(ni,nj,end) .* [(ramplen-1):-1:0]/ramplen;
    end
end

gradarea = sum(g1(:,:,ctr:end),3)*dt;

gradarea2 = -(sum(rampup,3) + sum(g(:,:,1:nseg),3))*dt; % starting spot to get to zero at TE

prewinder_cell = cell([3,nplanes]);
max_len = 0;
for ni = 1:3
    for nj = 1:nplanes
        prewinder_cell{ni,nj} = dotrap(gradarea(ni,nj)-gradarea2(ni,nj),slewramp,5,dt); % blip from end of last one to start of next one
        max_len = max(max_len,length(prewinder_cell{ni,nj}));
    end
end

prewinder = zeros(3, nplanes, max_len);
for ni = 1:3
    for nj = 1:nplanes
        tmp_len = length(prewinder_cell{ni,nj});
        prewinder(ni,nj,1:tmp_len) = flip(prewinder_cell{ni,nj});
        prewinder(ni,nj,:) = flip(prewinder(ni,nj,:));
    end
end

g2 = cat(3, -prewinder, rampup, g, rampdown);

ind_prewinder = 1:size(prewinder,3);
ind_rampup = (1:size(rampup,3))+size(prewinder,3);
ind_grad = (1:size(rampup,3))+size(prewinder,3);

len_prewinder = size(prewinder,3);
len_rampup = size(rampup,3);
len_rampdown = size(rampdown,3);
len_inward = nseg + 1;

gtot = cat(3,g1,g2);

karray_fid = zeros([size(g1),nspokes]);
for ispoke = 1:nspokes
    st_pt = ctr + (ispoke - 1) * nseg;
    karray_fid(:,:,st_pt:end,ispoke) =  gambar.*(cumsum(g1(:,:,st_pt:end),3)*dt)*fov; % unit are not cm-1, but 1/FOV or sample number
end

karray_gre = zeros([size(g2),nspokes]);
ktrim = zeros([size(g2,[1,2]),nseg*nspokes,nspokes]);
% tarray_gre = zeros(length(g2),nspokes);
t = [1:nspokes*nseg]*dt;  % start a inward moving spoke
for lp = 1:nspokes
    sttime = (lp-1)*tseg;
    stsamp = round(sttime/dt)+1 ;
    karray_gre(:,:,:,lp) = karray_fid(:,:,end,lp) + gambar.*(cumsum(g2,3)*dt)*fov; % unit are not cm-1, but 1/FOV or sample number
    ktrim(:,:,:,lp) = karray_gre(:,:,len_prewinder  + len_rampup + len_inward:end-len_rampdown,lp); % units are not cm-1, but 1/FOV or sample number
%     tarray_gre(:,lp) = dt*length(g1(ctr:end)) + dt*(1:length(g2)) - nseg * dt * (lp - 1);
end

figure;
ktmp1 = squeeze(ktrim(1,1,:,:));
ktmp2 = squeeze(ktrim(2,1,:,:));
ktmp3 = squeeze(ktrim(3,1,:,:));

plot3(ktmp1,ktmp2,ktmp3,'.');

figure;
subplot(311);
plot(squeeze(gtot(1,:,:)).');
xlabel('time');
ylabel('gx');
title('Gradient in TOPPE (gauss/cm)')
subplot(312);
plot(squeeze(gtot(2,:,:)).');
xlabel('time');
ylabel('gy');
subplot(313);
plot(squeeze(gtot(3,:,:)).');
xlabel('time');
ylabel('gz');

ktmp = gambar .* cumsum(gtot(:,:,(ctr):end),3) * dt;
figure;
subplot(311);
plot(squeeze(ktmp(1,:,:)).');
xlabel('time');
ylabel('kx');
title('k traj in TOPPE (1/fov)')
subplot(312);
plot(squeeze(ktmp(2,:,:)).');
xlabel('time');
ylabel('ky');
subplot(313);
plot(squeeze(ktmp(3,:,:)).');
xlabel('time');
ylabel('kz');

figure;
subplot(211);
plot(ctr:size(gtot,3),squeeze(gtot(1,:,ctr:end)).');
xlabel('time');
ylabel('gx');
title('Gradient in TOPPE (gauss/cm)')
subplot(212);
plot(ctr:size(gtot,3),squeeze(ktmp(1,:,:)).');
xlabel('time');
ylabel('kx');
title('k traj in TOPPE (1/fov)')

close all