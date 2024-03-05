
% units: cm, cm-1, ms, khz, Gauss
% constants
gambar = 4.258; % 4.258 kHz/g
gam = gambar*2*pi;

% pulse seq parameters
fov = 22; % cm
dt = 0.004; % 4 us
nspokes = 31; % # of RF pulses and therefore, echos
nseg = round(210); % samples
tseg = nseg*dt; % 840 us per Dionisio paper
trf = 3; % samples, 12 us per Haowei's simulations
rfpre = 1; % samples, 8 us - amount before "center of RF"
slewramp = 0.5; % g/cm/ms - slew for ramp-up/down to looping.
maxrf = 0.2; % g (needs to be refined)
flipangle = 2*pi/180; % radians
rfshape = [1 1 1];
rfamp = flipangle/(sum(rfshape)*dt*gam);
% setting sinusoidal grads (for simplicity) - this may not be exactly what they did
xres = 64;
res = fov/xres;
kmax1 = 1/res/2; % max k-space for each segment, cm-1
gamp = kmax1/gambar/tseg; % g/cm

% Segment 1
tmax = nspokes*tseg; %might only need 64, but for completeness
t = [-rfpre:((nspokes-1)*nseg + trf)]*dt;  %start at RF pulse and end just after last RF pulse
% t = [-rfpre:((nspokes)*nseg + trf)]*dt;  %start at RF pulse and end just after last RF pulse

% start with a Fourier polygon in k and convert to g
% reference
% https://demonstrations.wolfram.com/FourierConstructionOfRegularPolygonsAndStarPolygons/#more
% norder = 1;
% nar = [-norder:norder];
% c = (nspokes./(nar*nspokes+1).*sin(pi/nspokes)).^2; c = c./max(c);
% % c(1) = c(1)/2; c(3)= c(3)/2;
% k = sum((c'*ones(size(t))).*exp(1i*om1*((nspokes*nar'+1)*t)),1);
% % rotate shift and scale
% % scale k so that it is at the edge k-space at the end of 1 seg
% k = -1i*(k - 1);
% k = k*kmax1/abs(k(round(tseg/dt)+1)); 
% 
% g = diff(k)/dt/gambar;
% g = [g g(end)]; %make the same length

om1 = 2*pi/(tmax);
g = gamp.*exp(1i*om1*t); % g/cm
ramplen = ceil(gamp/slewramp/dt);
rampup = [0:(ramplen-1)]*g(1)/ramplen;
rampdown = [(ramplen-1):-1:0]*g(end)/ramplen;
g1 = [rampup g rampdown];
tu = t(1) + [-ramplen:-1]*dt;
td = t(end) + [1:ramplen]*dt;
t = [ tu t td]; 
ctr = ramplen + rfpre + 1;
k = gambar.*cumsum(g1(ctr:end))*dt*fov; % unit are not cm-1, but 1/FOV or sample number
s1 = [0 diff(g1)]/dt; % g/cm/ms
gradarea = sum(g1(ctr:end))*dt;
te1 = length(g1(ctr:end))*dt;
for lp = 1:nspokes
    gradareaar(lp) = sum(g1((ctr + (lp-1)*nseg):end))*dt;
    te1ar(lp) = length(g1((ctr + (lp-1)*nseg):end))*dt;
end
k1 = gambar.*(cumsum(g1(ctr:end))*dt)*fov; % unit are not cm-1, but 1/FOV or sample number
karray_fid = zeros(length(g1),nspokes);
for ispoke = 1:nspokes
    st_pt = ctr + (ispoke - 1) * nseg;
    karray_fid(st_pt:end,ispoke) =  gambar.*(cumsum(g1(st_pt:end))*dt)*fov; % unit are not cm-1, but 1/FOV or sample number
end

figure;
plot(squeeze(real(karray_fid(:,:))),squeeze(imag(karray_fid(:,:))),'.');

xlabel('k_x')
ylabel('k_y')
title('max grad 1 k-space trajectory')

% build RF pulse
rf = zeros(size(g1));
for lp = 1:nspokes
    stloc = ramplen+1+(lp-1)*nseg;
    rf(stloc:(stloc+trf-1)) = rfshape*rfamp;
end

slewx = gradient(real(g));
slewy = gradient(imag(g));

figure;
subplot(211);
plot(slewx*1e3/4);
xlabel('time');
ylabel('slew rate in x');
title('Slew Rate gauss/cm/ms')
subplot(212);
plot(slewy*1e3/4);
xlabel('time');
ylabel('slew rate in y');

%% Create gre module
t = [-nseg:nspokes*nseg]*dt;  % start a inward moving spoke
g = gamp.*exp(1i*om1*t); % g/cm
% k = sum((c'*ones(size(t))).*exp(1i*om1*((nspokes*nar'+1)*t)),1);
% % rotate shift and scale
% % scale k so that it is at the edge k-space at the end of 1 seg
% k = -1i*(k - 1);
% k = k*kmax1/abs(k(round(tseg/dt)+1)); 
% g = diff(k)/dt/gambar;
% g = [g g(end)]; %make the same length

ramplen = ceil(gamp/slewramp/dt);
rampup = [0:(ramplen-1)]*g(1)/ramplen;
rampdown = [(ramplen-1):-1:0]*g(end)/ramplen;
gradarea2 = -(sum(rampup) + sum(g(1:nseg)))*dt; % starting spot to get to zero at TE
prewinder = dotrap(gradarea-gradarea2,slewramp,5,dt); % blip from end of last one to start of next one
lprew = length(prewinder);
g2 = [-prewinder rampup g rampdown];
t = [1:length(g2)]*dt;
ctr2 = ramplen + lprew + nseg
ssitime = 0.110;  % 110 us ???
te = te1 + ssitime + ctr2*dt
k2 = gambar.*(cumsum(g2)*dt)*fov; % unit are not cm-1, but 1/FOV or sample number
k = k1(end) + k2; % unit are not cm-1, but 1/FOV or sample number
s2 = [0 diff(g2)]/dt; % g/cm/ms

te = nspokes*tseg; % ms - this sets omega1
l1 = te/dt;
% now specify k-space for 32 RF pulses (assume overlapping)
lk = length(k);

len_prewinder = length(prewinder); % per looping_star.m
len_rampup = length(rampup); % per looping_star.m
len_rampdown = length(rampdown); % per looping_star.m
len_inward = nseg+1;
karray = k1(end) + gambar.*cumsum(g2(1:end))*dt*fov;
ktrim = zeros(nspokes,nseg*nspokes);
karray_gre = zeros(length(g2),nspokes);
tarray_gre = zeros(length(g2),nspokes);
t = [1:nspokes*nseg]*dt;  % start a inward moving spoke
g = gamp.*exp(1i*om1*t); % g/cm
for lp = 1:nspokes
    sttime = (lp-1)*tseg;
    stsamp = round(sttime/dt)+1 ;
    karray_gre(:,lp) = karray_fid(end,lp) + gambar.*(cumsum(g2)*dt)*fov; % unit are not cm-1, but 1/FOV or sample number
    ktrim(lp,:) = karray_gre(len_prewinder  + len_rampup + len_inward+1:end-len_rampdown,lp); % units are not cm-1, but 1/FOV or sample number
    tarray_gre(:,lp) = dt*length(g1(ctr:end)) + dt*(1:length(g2)) - nseg * dt * (lp - 1);
end


if 0 
    gtot = cat(2,g1,g2);
    ktmp = gambar * real(cumsum(gtot(ctr:end),2)) * dt;
    figure;
    subplot(211);
    plot(ctr:size(gtot,2),squeeze(real(gtot(ctr:end))).');
    xlabel('time');
    ylabel('gx');
    title('Gradient in TOPPE (gauss/cm)')
    subplot(212);
    plot(ctr:size(gtot,2),squeeze(ktmp).');
    xlabel('time');
    ylabel('kx');
    title('k traj in TOPPE (1/fov)')
end
%% Make simulation of k-space trajectory with gradient ramp up/down
% t2star = 40;
% sim_sig = fshepp2(real(karray_gre),imag(karray_gre)) .* exp(-tarray_gre/t2star);
% sig = sum(sim_sig,2);
