% recon_looping_star.m


% reconstruct looping star images and using sensitivity
% maps from BART

% Haowei Xiang
% Mar 2024

%% Read path
import toppe.*
import toppe.utils.*
%% Load data
phan = '01-24-24-human';
lspath = ['./',phan,'/'];
sapath_gre = ['./',phan,'/'];

switch phan
    case '01-24-24-human'
        lsname = 'P46592.7'; 
end
% reconstruct (magnitude) images
echo = 1;
[dat, rdb_hdr] = toppe.utils.loadpfile(lsname,echo);
[nk, ncoil, nsli, ntp, nview] = size(dat);
%% Extract sensitivity maps
%% coil compress raw data
sfname = 'nc6_te0_fov24_im80_1a_smaps_gre.mat';
load(sfname,'sens_3d_gre_te1','vr_3d');
sensmaps = sens_3d_gre_te1;
sensmaps = flip(flip(flip(sensmaps,1),2),3);

ncc = size(sensmaps,4);
dat_2d = permute(dat,[1 3 4 5 2]); % permute coil index to the last
dat_2d = reshape(dat_2d, [nk*ntp*nsli*nview,ncoil]); % reshape for matrix multiplication 
dat_2d_cc = dat_2d * vr_3d; % compress coil
dat_cc_rs = reshape(dat_2d_cc, [nk,nsli,ntp,nview,ncc]); % reshape 2d to 5d dat
dat_cc = permute(dat_cc_rs, [1 5 2 3 4]); % permute the coil index to 2nd index
dat = dat_cc;
clear dat_cc dat_cc_rs dat_2d_cc dat_2d
% extract the 3rd 3D acquisition to reach steady state
nrun_idle = 2; % number of runs to reach steady-state
ind_rep = (nrun_idle+1):size(dat,5); % use first 2 to reach steady-state
dat = squeeze(dat(:,:,:,:,ind_rep)); % 
dat = flip(dat,1); % flip dat spoke because unknown reason
% plot k-data from the largest sinular value
if 0 
    dat_tmp = squeeze(dat(:,1,:,1));
    figure; plot(abs(dat_tmp(:,:)));
    figure; plot(angle(dat_tmp(:,:)));
end

%% Extract corresponding k-space locations
run ./looping_star_cv521_grad.m

k_tot = zeros(3,nplanes,nseg,nspokes,nspokes);
for lp = 1:nspokes % looping over number of RF pulses
    ind = (lp -1) * nseg + (1:nseg); % time point index

    tmp = circshift(ktrim(:,:,ind,:),[0 0 0 nspokes-lp+1]); % circshift to shift spoke index
    k_tot(:,:,:,lp,:) = tmp; % k-space locations for all spokes and all sections
end
k1 = squeeze(k_tot(:,:,:,:,1));
k2 = squeeze(k_tot(:,:,:,:,2));

k1x_tot = permute(squeeze(k1(1,:,:,:)),[2 3 1]);
k1y_tot = permute(squeeze(k1(2,:,:,:)),[2 3 1]);
k1z_tot = permute(squeeze(k1(3,:,:,:)),[2 3 1]);
k2x_tot = permute(squeeze(k2(1,:,:,:)),[2 3 1]);
k2y_tot = permute(squeeze(k2(2,:,:,:)),[2 3 1]);
k2z_tot = permute(squeeze(k2(3,:,:,:)),[2 3 1]);

%% arguments for Gmri object
fov = 24; % cm
xres = 120;
ig = image_geom('nx', xres, 'nz', xres, ...
    'offsets', 'dsp', ... % [-n/2:n/2-1] for mri
    'fov', [fov fov fov]); % 20 cm transaxial FOV
N = ig.dim;
ig.mask = true(N);
nufft_args = {N, 6*ones(size(N)), 2*N, N/2, 'table', 2^10, 'minmax:kb'};
x0 = zeros(N);
niter = 30;

%% Extract valid looping star data
len_prewinder = length(prewinder);
len_rampup = length(rampup);
len_inward = nseg + 1;

st_pts = len_prewinder + len_rampup + len_inward ; % approximately len_prewinder  + len_rampup + len_inward, but need to shift due to hardware delay
st_pts = st_pts ; % manually corrected by looking at the phase
% ktrim = karray(:,len_prewinder  + len_rampup + len_inward+1:end - len_rampdown);
dat3mc = dat(st_pts + 1: st_pts + nseg*nspokes,:,:,:);
%     clear dat
dat3mc = permute(dat3mc,[1 3 4 2]); % [nk, nspoke, ncoil];
ntime = size(dat3mc,4);

kpercent = 0.5; % how much data to use with IFT
klen = nseg;
% ind_klen = floor(klen * kpercent/2):klen * kpercent;
ind_klen = 1:klen * kpercent;
ind_kspk = 1:nspokes-1;
ind_kp = 1:nplanes;
ind_kview = 1;

%% Gridding using sparse matrix
GOF = 2; % gridding over-sampling factor
Ni=GOF*ig.nx; 

%% grid with echo-in spoke
ind_klen = 1:klen * kpercent;

KXr=round(GOF*k1y_tot(ind_klen,ind_kspk,ind_kp)); 
KYr=round(GOF*k1x_tot(ind_klen,ind_kspk,ind_kp)); 
KZr=round(GOF*k1z_tot(ind_klen,ind_kspk,ind_kp));

ind_klen = (klen * kpercent+1): klen;

% tmpD=sum(abs(dat3mc(ind_klen,ind_kspk,ind_kp,ind_kview,:)),4);
KXr2=round(GOF*k2y_tot(ind_klen,ind_kspk,ind_kp)); 
KYr2=round(GOF*k2x_tot(ind_klen,ind_kspk,ind_kp)); 
KZr2=round(GOF*k2z_tot(ind_klen,ind_kspk,ind_kp));

ix=[KXr(:);KXr2(:)]+Ni/2+1; iy=[KYr(:);KYr2(:)]+Ni/2+1; iz=[KZr(:);KZr2(:)]+Ni/2+1;
indi=(iz-1)*Ni^2+(iy-1)*Ni+ix;

% 
iNonZero = (1:numel([KXr(:);KXr2(:)]));
nNonZero = (numel([KXr(:);KXr2(:)]));

% Sparse Gridding Matrix (GM) and Density Compensation (DC)
GM=sparse(indi,iNonZero,ones(nNonZero,1),Ni^3,nNonZero);
DC=reshape(full(GM*sparse(ones(nNonZero,1))),[Ni,Ni,Ni]); 
DC=smooth3(DC,'gauss',2*ceil(GOF/2)+1); MSK=DC>1e-6;

ii=Ni/2+(-Ni/2+1:Ni/2); Nii=length(ii);
BILD=zeros(Nii,Nii,Nii,ncc,'single');
ims_ini = Ni/2 + (-Ni/4+1:Ni/4);

BILD_sens=zeros([ig.dim,ncc],'single');

for icoil=1:ncc
    dat_tmp = reshape(dat3mc(:,:,:,icoil),[nseg,nspokes,nplanes,nview-nrun_idle]);
    dtmp1 = dat_tmp(1:klen * kpercent,ind_kspk,ind_kp,ind_kview);
    dtmp2 = dat_tmp((klen * kpercent+1): klen,ind_kspk,ind_kp,ind_kview);
    
%     dtmp1 = rand(size(dtmp1)) + 1i * rand(size(dtmp1));
%     dtmp2 = rand(size(dtmp2)) + 1i * rand(size(dtmp2));
    DAT=sparse([dtmp1(:);dtmp2(:)]);
    CK=reshape(full(GM*DAT),[Ni,Ni,Ni]);
    CK(MSK)=CK(MSK)./DC(MSK);
    TMP=single(fftshift(fftn(fftshift(CK))));
    TMP = permute(TMP,[2 1 3]);
    sens_map_i = imresize3(sensmaps(:,:,:,icoil),ig.dim); % conjugate phase ??? why need it to work
    BILD(:,:,:,icoil)= sqrt(TMP(ii,ii,ii).*conj(TMP(ii,ii,ii)));
    BILD_sens(:,:,:,icoil)= conj(sens_map_i) .*TMP(ims_ini,ims_ini,ims_ini);
%     icoil
end

grid_echo = abs(sum(BILD,4));
grid_echo_sens = abs(sum(abs(BILD_sens(:,:,:,:)),4));

grid_echo = grid_echo(ims_ini,ims_ini,ims_ini);

figure;im((grid_echo_sens(:,:,:)),'gridding recon with sense map');colormap(gray);cbar

%% Model-based reconstruction
st_pts = len_prewinder + len_rampup + len_inward; % approximately len_prewinder  + len_rampup + len_inward, but need to shift due to hardware delay
st_pts = st_pts ; % manually corrected by looking at the phase
dat3mc = dat(st_pts + 1: st_pts + nseg*nspokes,:,:,:);
%     clear dat
dat3mc = permute(dat3mc,[1 3 4 2]); % [nk, nspoke, ncoil];
ntime = size(dat3mc,3);

if 0
    dat3mc = zeros(size(dat3mc));
end
%% Construct system matrix
kpercent = 1; % how much data to use with IFT
klen = nseg;
ind_klen = 1:klen * kpercent;
ind_klen2 = (klen * kpercent+1):klen;
ind_kspk = 1:nspokes-1; % choose the last spoke since its with double TE
ind_kp = 1:nplanes;
% ind_kview = 1:nview-nrun_idle;
ind_kview = 1;

%% 1st system matrix
kxtmp = -k1x_tot(ind_klen,ind_kspk,ind_kp,ind_kview);
kytmp = -k1y_tot(ind_klen,ind_kspk,ind_kp,ind_kview);
kztmp = -k1z_tot(ind_klen,ind_kspk,ind_kp,ind_kview);

kspace = [kxtmp(:),kytmp(:),kztmp(:)]/fov; % from 1/fov to 1/cm
Gm1 = Gmri(kspace, ig.mask, ...
'fov', ig.fovs, 'basis', {'rect'}, 'nufft', nufft_args);
%% 2nd system matrix
kxtmp = -k2x_tot(ind_klen,ind_kspk,ind_kp,ind_kview);
kytmp = -k2y_tot(ind_klen,ind_kspk,ind_kp,ind_kview);
kztmp = -k2z_tot(ind_klen,ind_kspk,ind_kp,ind_kview);

kspace = [kxtmp(:),kytmp(:),kztmp(:)]/fov; % from 1/fov to 1/cm
Gm2 = Gmri(kspace, ig.mask, ...
'fov', ig.fovs, 'basis', {'rect'}, 'nufft', nufft_args);
%% System matrix for multiple coil
Gm_mc_cell = cell(ncc,1);
Gm_mc_cell2 = cell(ncc,1);
for icoil = 1:ncc
    sens_map_i = imresize3(sensmaps(:,:,:,icoil),[xres,xres,xres]);
    tmp = Gdiag(reshape(sens_map_i(ig.mask),[],1),'mask',ig.mask);
    Gm_mc_cell{icoil} = Gm1 * tmp;
    Gm_mc_cell2{icoil} = Gm2 * tmp;
end
Gm_tot = cat(1, Gm_mc_cell{:});
Gm_tot2 = cat(1, Gm_mc_cell2{:});

dat_mc = cell(ncc,1);
dat_mc2 = cell(ncc,1);
for ic = 1:ncc
    dat_tmp = reshape(dat3mc(:,:,:,ic),[nseg,nspokes,nplanes,nview-nrun_idle]);
    dat_tmp2 = mean(dat_tmp(ind_klen,ind_kspk,ind_kp,1:6),4); % data for echo-out spokes
    dat_tmp2 = dat_tmp2(:);
    
    dat_mc{ic} = dat_tmp2;
end
dat_tot = cat(1, dat_mc{:});
%% With Gaussian window
alpha = 3;
wtmp = gausswin(nseg,alpha);
w = repmat(wtmp, [1, nspokes, nplanes, nview, ncc]);
w = w(ind_klen,ind_kspk,ind_kp,ind_kview,:); 
w = reshape(w,[],ncc);
W = Gdiag(w(:));

beta = 2^-12 * numel(dat_tot); % good for quadratics

Rq = Reg1(ig.mask, 'beta', beta); % beta is divided by 2 because seperate recon has half signal length compared to joint recon
tmp2 = qpwls_psf(Gm_tot + Gm_tot2, Rq, 1, ig.mask, W, 'fwhmtype', 'profile');
%         xpcg_mc2 = qpwls_pcg(x0(ig.mask), Gm_tot + Gm_tot2, W, dat_tot(:), 0, Rq.C, 1, niter, ig.mask);
xpcg_mc2 = qpwls_pcg(x0(ig.mask), Gm_tot + Gm_tot2, W, dat_tot(:), 0, Rq.C, 1, niter, ig.mask);
mb2_ims = ig.embed(xpcg_mc2(:,end));
  

figure;
im(mb2_ims(:,:,:),'Model-based reconstruction'); colormap(gray)