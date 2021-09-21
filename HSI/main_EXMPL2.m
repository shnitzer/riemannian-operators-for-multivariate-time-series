clear; clc; close all;
load('casi_2013.mat');
load('lidar_2013.mat');

ldr_im = x2013_IEEE_GRSS_DF_Contest_LiDAR;
hps_im = x2013_IEEE_GRSS_DF_Contest_CASI;

fig1  = figure(100); set(gcf,'Position',[83,400,1337,235]);
ax1 = axes('Parent',fig1);
imagesc(ldr_im); axis equal;

fig2  = figure(200); set(gcf,'Position',[83,400,1337,235]);
ax2 = axes('Parent',fig2);
imagesc(mean(hps_im,3)); axis equal;
linkaxes([ax1,ax2],'xy');

rng_x = [151,210]; rng_y = [321,410];

%%
use_spat_nn =0;

neig = 20;

if use_spat_nn
spat_neig     = 20;         
[Xm,Ym]       = meshgrid(1:(rng_y(2)-rng_y(1)+1),1:(rng_x(2)-rng_x(1)+1));
spat_dist     = squareform(pdist([Xm(:),Ym(:)]));
spat_dist_bin = spat_dist<spat_neig;
spat_dist     = spat_dist.*(spat_dist>spat_neig);
spat_dist     = 0.5*(spat_dist + spat_dist.');
end

c_ldr = double(ldr_im(rng_x(1):rng_x(2),rng_y(1):rng_y(2)));
c_hps = double(hps_im(rng_x(1):rng_x(2),rng_y(1):rng_y(2),:));

figure(50)
subplot(231);imagesc(c_ldr); title('ldr');                  colorbar;
subplot(232);imagesc(c_hps(:,:,1)); title('hps - 1');       colorbar;
subplot(233);imagesc(mean(c_hps,3)); title('hps - avg');    colorbar;

prc_ldr = [0,99];
prc_hps = [5,95];
nrm_c_ldr        = c_ldr/std(reshape(c_ldr,1,[]));
[~,TF_ldr,~,~,~] = filloutliers(nrm_c_ldr(:),'nearest','percentile',prc_ldr);
c_ldr            = regionfill(nrm_c_ldr,reshape(TF_ldr,size(c_ldr)));

nrm_c_hps        = bsxfun(@rdivide,reshape(c_hps,numel(c_ldr),[]),std(reshape(c_hps,numel(c_ldr),[])));
[~,TF_hps,~,~,~] = filloutliers(nrm_c_hps,'nearest','percentile',prc_hps);
for ii = 1:size(c_hps,3)
    c_hps(:,:,ii) = regionfill(reshape(nrm_c_hps(:,ii),size(c_ldr)),reshape(TF_hps(:,ii),size(c_ldr)));
end

figure(50)
subplot(234);imagesc(c_ldr); title('ldr');                  colorbar;
subplot(235);imagesc(c_hps(:,:,1)); title('hps - 1');       colorbar;
subplot(236);imagesc(mean(c_hps,3)); title('hps - avg');    colorbar;

%%%% for paper: %%%%
figure(1000)
imagesc(c_ldr);
title('LiDAR');
set(gcf,'Position',[1000,100,400,240]);
axis equal; axis tight;
colorbar;

figure(2000)
imagesc(mean(c_hps,3));
title('avg(HSI)');
set(gcf,'Position',[1000,100,400,240]);
axis equal; axis tight;
colorbar;
%%%%%%%%%%%%%%%%%%%%


DIS_ldr = squareform(pdist(c_ldr(:)));
DIS_hps = squareform(pdist(reshape(c_hps,numel(c_ldr),[])));

if use_spat_nn
    DIS_hps_l = (DIS_hps/std(DIS_hps(:))).*(1 + spat_dist/max(spat_dist(:)));
    DIS_ldr_l = (DIS_ldr/std(DIS_ldr(:))).*(1 + spat_dist/max(spat_dist(:)));
    
    nnbrs   = 500;
    alpha   = 2;
    
    [min_vec_ldr,min_ind_ldr] = mink(DIS_ldr_l,nnbrs,2);
    bin_mat_neig_ldr          = sparse(reshape(repmat((1:size(DIS_ldr_l,1)).',1,nnbrs),[],1),min_ind_ldr(:),1,size(DIS_ldr_l,1),size(DIS_ldr_l,2));
    bin_mat_neig_ldr          = logical(bin_mat_neig_ldr + bin_mat_neig_ldr.');
    nbrs_ldr_mat              = ones(size(bin_mat_neig_ldr)); % spat_dist_bin; % bin_mat_neig_ldr; % 
    
    [min_vec_hps,min_ind_hps] = mink(DIS_hps_l,nnbrs,2);
    bin_mat_neig_hps          = sparse(reshape(repmat((1:size(DIS_hps_l,1)).',1,nnbrs),[],1),min_ind_hps(:),1,size(DIS_hps_l,1),size(DIS_hps_l,2));
    bin_mat_neig_hps          = logical(bin_mat_neig_hps + bin_mat_neig_hps.');
    nbrs_hps_mat              = ones(size(bin_mat_neig_hps)); % spat_dist_bin; % bin_mat_neig_hps; % 

else
    DIS_hps_l = (DIS_hps/std(DIS_hps(:)));
    DIS_ldr_l = (DIS_ldr/std(DIS_ldr(:)));
end

if use_spat_nn
    K_ldr   = exp(-DIS_ldr_l.^2/(0.05*median(DIS_ldr_l(nbrs_ldr_mat(:)~=0).^2)));
    K_ldr   = K_ldr.*nbrs_ldr_mat;
else
    K_ldr   = exp(-DIS_ldr_l.^2/(1*median(DIS_ldr_l(:).^2)));
end
if 1
K_ldr_dm   = diag(1./(sum(K_ldr,1))) * K_ldr * diag(1./(sum(K_ldr,1)));
K_ldr_dm   = diag(1./sqrt(sum(K_ldr_dm,1))) * K_ldr_dm * diag(1./sqrt(sum(K_ldr_dm,1)));
[V_dm_ldr,D_dm_ldr] = eigs(K_ldr_dm,neig);
DMdist_ldr          = squareform(pdist(V_dm_ldr * D_dm_ldr));
end
K_ldr   = diag(1./(sum(K_ldr,1))) * K_ldr * diag(1./(sum(K_ldr,1)));
K_ldr   = diag(1./sqrt(sum(K_ldr,1))) * K_ldr * diag(1./sqrt(sum(K_ldr,1)));
K_ldr   = 0.5* (K_ldr + K_ldr.');

if use_spat_nn
    K_hps   = exp(-DIS_hps_l.^2/(1*median(DIS_hps_l(nbrs_hps_mat(:)~=0).^2)));
    K_hps   = K_hps.*nbrs_hps_mat;
else
    K_hps   = exp(-DIS_hps_l.^2/(0.5*median(DIS_hps_l(:).^2)));
end
if 1
K_hps_dm   = diag(1./(sum(K_hps,1))) * K_hps * diag(1./(sum(K_hps,1)));
K_hps_dm   = diag(1./sqrt(sum(K_hps_dm,1))) * K_hps_dm * diag(1./sqrt(sum(K_hps_dm,1)));
[V_dm_hps,D_dm_hps] = eigs(K_hps_dm,neig);
DMdist_hps          = squareform(pdist(V_dm_hps * D_dm_hps));
end
K_hps   = diag(1./(sum(K_hps,1))) * K_hps * diag(1./(sum(K_hps,1)));
K_hps   = diag(1./sqrt(sum(K_hps,1))) * K_hps * diag(1./sqrt(sum(K_hps,1)));
K_hps   = 0.5* (K_hps + K_hps.');

max(K_ldr(:))
max(K_hps(:))

tmp = eigs(K_ldr,50);
tmp2 = eigs(K_hps,50);


%% AD

% dim_ad = 11;
% % K_ldr_ad = diag(1./sum(K_ldr,1)) * K_ldr;
% % K_hps_ad = diag(1./sum(K_hps,1)) * K_hps;
% 
% K_ldr_ad = K_ldr;
% K_hps_ad = K_hps;
% 
% P = K_ldr_ad * K_hps_ad;
% 
% DIS_ad = squareform(pdist(P));
% K_ad   = exp(-DIS_ad.^2/median(DIS_ad(:)));
% K_ad   = diag(1./sum(K_ad,1)) * K_ad;
% 
% [Vad,Dad] = eigs(K_ad,dim_ad);
% 
% dad = diag(Dad);
% 
% gindx_ad  = sum(abs(Vad)*diag(1./(max(abs(Vad),[],1)*sqrt(dad))),2);
% 
% figure
% imagesc(reshape(gindx_ad,size(c_ldr)))
%% S & A
dim = neig;
t_S = 0.5;

S = FixedGeodes_eff( K_hps, K_ldr, t_S, dim );
A = FixedGeodes_eff_proj( S, K_ldr, dim );  % A neg refers to hps in this case
[VS,DS] = eigs(S,dim);
[VA,DA] = eigs(A,dim);
tmpVA = VA; tmpDA = DA;

dS      = diag(DS);
dA      = diag(DA);
dA_pos  = dA(dA > 0);
VA_pos  = VA(:,dA > 0);

A = FixedGeodes_eff_proj( S, K_hps, dim );

[VA,DA] = eigs(A,dim);
dA      = diag(DA);
dA_neg  = abs(dA(dA > 0));
VA_neg  = VA(:,dA > 0);

%%
% Calculating the function from - "The geometry of nodal sets and outlier
% detection"

dim_g    = min(size(VA_pos,2),size(VA_neg,2));
dim_g_dm = neig;
g_rng    = 1:dim_g;
g_rng_dm = 1:dim_g_dm;
gindx_S  = sum(abs(VS(:,1:dim_g))   *diag(1./(max(abs(VS(:,g_rng)),[],1)    *sqrt(dS(g_rng)))),2);
gindx_A1 = sum(abs(VA_pos(:,g_rng)) *diag(1./(max(abs(VA_pos(:,g_rng)),[],1)*sqrt(dA_pos(g_rng)))),2);
gindx_A2 = sum(abs(VA_neg(:,g_rng)) *diag(1./(max(abs(VA_neg(:,g_rng)),[],1)*sqrt(dA_neg(g_rng)))),2);
% gindx_A  = sum(abs(VA)*diag(1./(max(abs(VA),[],1)*sqrt(abs(dA)))),2);


vind = 0;
figure(300)
subplot(1,3,1); imagesc(reshape(VS(:,vind+1)  ,size(c_ldr))); colorbar; title('$$V_S^{(1)}$$','Interpreter','latex');
subplot(1,3,2); imagesc(reshape(VS(:,vind+2)  ,size(c_ldr))); colorbar; title('$$V_S^{(2)}$$','Interpreter','latex');
subplot(1,3,3); imagesc(reshape(VS(:,vind+3)  ,size(c_ldr))); colorbar; title('$$V_S^{(3)}$$','Interpreter','latex');

figure(301)
subplot(1,3,1); imagesc(reshape(VA_pos(:,vind+1)  ,size(c_ldr))); colorbar; title('$$V_{A+}^{(1)}$$','Interpreter','latex');
subplot(1,3,2); imagesc(reshape(VA_pos(:,vind+2)  ,size(c_ldr))); colorbar; title('$$V_{A+}^{(2)}$$','Interpreter','latex');
subplot(1,3,3); imagesc(reshape(VA_pos(:,vind+3)  ,size(c_ldr))); colorbar; title('$$V_{A+}^{(3)}$$','Interpreter','latex');

figure(302)
subplot(1,3,1); imagesc(reshape(VA_neg(:,vind+1)  ,size(c_ldr))); colorbar; title('$$V_{A-}^{(1)}$$','Interpreter','latex');
subplot(1,3,2); imagesc(reshape(VA_neg(:,vind+2)  ,size(c_ldr))); colorbar; title('$$V_{A-}^{(2)}$$','Interpreter','latex');
subplot(1,3,3); imagesc(reshape(VA_neg(:,vind+3)  ,size(c_ldr))); colorbar; title('$$V_{A-}^{(3)}$$','Interpreter','latex');

figure
subplot(1,3,1); imagesc(reshape(gindx_S,size(c_ldr)));  colorbar; title('$$f_N(S)$$','Interpreter','latex');
subplot(1,3,2); imagesc(reshape(gindx_A1,size(c_ldr))); colorbar; title('$$f_N(A+)$$','Interpreter','latex');
subplot(1,3,3); imagesc(reshape(gindx_A2,size(c_ldr))); colorbar; title('$$f_N(A-)$$','Interpreter','latex');

d_dm_hps = diag(D_dm_hps);
d_dm_ldr = diag(D_dm_ldr);
gindx_ldr = sum(abs(V_dm_ldr(:,g_rng_dm)) *diag(1./(max(abs(V_dm_ldr(:,g_rng_dm)),[],1)*sqrt(d_dm_ldr(g_rng_dm)))),2);
gindx_hps = sum(abs(V_dm_hps(:,g_rng_dm)) *diag(1./(max(abs(V_dm_hps(:,g_rng_dm)),[],1)*sqrt(d_dm_hps(g_rng_dm)))),2);

figure(400)
subplot(1,2,1); imagesc(reshape(gindx_ldr,size(c_ldr))); colorbar;
subplot(1,2,2); imagesc(reshape(gindx_hps,size(c_ldr))); colorbar;


figure(303)
subplot(1,3,1); imagesc(reshape(V_dm_ldr(:,vind+1)  ,size(c_ldr))); colorbar; title('$$V_{LDR}^{(1)}$$','Interpreter','latex');
subplot(1,3,2); imagesc(reshape(V_dm_ldr(:,vind+2)  ,size(c_ldr))); colorbar; title('$$V_{LDR}^{(2)}$$','Interpreter','latex');
subplot(1,3,3); imagesc(reshape(V_dm_ldr(:,vind+3)  ,size(c_ldr))); colorbar; title('$$V_{LDR}^{(3)}$$','Interpreter','latex');

figure(304)
subplot(1,3,1); imagesc(reshape(V_dm_hps(:,vind+1)  ,size(c_ldr))); colorbar; title('$$V_{HPS}^{(1)}$$','Interpreter','latex');
subplot(1,3,2); imagesc(reshape(V_dm_hps(:,vind+2)  ,size(c_ldr))); colorbar; title('$$V_{HPS}^{(2)}$$','Interpreter','latex');
subplot(1,3,3); imagesc(reshape(V_dm_hps(:,vind+3)  ,size(c_ldr))); colorbar; title('$$V_{HPS}^{(3)}$$','Interpreter','latex');


%%
%%%%% for paper %%%%%
figure
imagesc(reshape(abs(VS(:,1))  ,size(c_ldr))); colorbar; title('$$|\psi_1^{(\mathbf{S})}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(VS(:,2))  ,size(c_ldr))); colorbar; title('$$|\psi_2^{(\mathbf{S})}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(VS(:,3))  ,size(c_ldr))); colorbar; title('$$|\psi_3^{(\mathbf{S})}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(VS(:,4))  ,size(c_ldr))); colorbar; title('$$|\psi_4^{(\mathbf{S})}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);

figure
imagesc(reshape(abs(VA_pos(:,1))  ,size(c_ldr))); colorbar; title('$$|\psi_1^{(\mathbf{A}^+)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(VA_pos(:,2))  ,size(c_ldr))); colorbar; title('$$|\psi_2^{(\mathbf{A}^+)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(VA_pos(:,3))  ,size(c_ldr))); colorbar; title('$$|\psi_3^{(\mathbf{A}^+)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);

figure
imagesc(reshape(abs(VA_neg(:,1))  ,size(c_ldr))); colorbar; title('$$|\psi_1^{(\mathbf{A}^-)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(VA_neg(:,2))  ,size(c_ldr))); colorbar; title('$$|\psi_2^{(\mathbf{A}^-)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(VA_neg(:,3))  ,size(c_ldr))); colorbar; title('$$|\psi_3^{(\mathbf{A}^-)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);

figure
imagesc(reshape(abs(V_dm_hps(:,1))  ,size(c_ldr))); colorbar; title('$$|\psi_1^{(2)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(V_dm_hps(:,2))  ,size(c_ldr))); colorbar; title('$$|\psi_2^{(2)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(V_dm_hps(:,3))  ,size(c_ldr))); colorbar; title('$$|\psi_3^{(2)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(V_dm_hps(:,4))  ,size(c_ldr))); colorbar; title('$$|\psi_4^{(2)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);

figure
imagesc(reshape(abs(V_dm_ldr(:,1))  ,size(c_ldr))); colorbar; title('$$|\psi_1^{(1)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(V_dm_ldr(:,2))  ,size(c_ldr))); colorbar; title('$$|\psi_2^{(1)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(V_dm_ldr(:,3))  ,size(c_ldr))); colorbar; title('$$|\psi_3^{(1)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);
figure
imagesc(reshape(abs(V_dm_ldr(:,4))  ,size(c_ldr))); colorbar; title('$$|\psi_4^{(1)}|$$','Interpreter','latex','FontSize',14);
axis equal; axis tight; set(gcf,'Position',[900,250,340,220]);


%%%%%%%%%%%%%
% saveas(gcf,'figures_for_paper\data2\Apos_psi_3_arrow.png')
% saveas(gcf,'figures_for_paper\data2\Apos_psi_3_arrow.fig')
% % saveas(gcf,'figures_for_paper\data2\Apos_psi_3_arrow','epsc')
%%%%%%%%%%%%%%%%%%%%%


%% Correlations:
figure
subplot(1,3,1); imagesc(abs((V_dm_hps.')*VS));      colorbar; colormap(gray); caxis([0,1]); title('corr(hsi,S)')
subplot(1,3,2); imagesc(abs((V_dm_hps.')*VA_pos));  colorbar; colormap(gray); caxis([0,1]); title('corr(hsi,A+)')
subplot(1,3,3); imagesc(abs((V_dm_hps.')*VA_neg));  colorbar; colormap(gray); caxis([0,1]); title('corr(hsi,A-)')
figure
subplot(1,3,1); imagesc(abs((V_dm_ldr.')*VS));      colorbar; colormap(gray); caxis([0,1]); title('corr(ldr,S)')
subplot(1,3,2); imagesc(abs((V_dm_ldr.')*VA_pos));  colorbar; colormap(gray); caxis([0,1]); title('corr(ldr,A+)')
subplot(1,3,3); imagesc(abs((V_dm_ldr.')*VA_neg));  colorbar; colormap(gray); caxis([0,1]); title('corr(ldr,A-)')
% figure
% subplot(1,3,1); imagesc(abs((V_dm_fus.')*VS));      colorbar; colormap(gray); caxis([0,1]); title('corr(fus,S)')
% subplot(1,3,2); imagesc(abs((V_dm_fus.')*VA_pos));  colorbar; colormap(gray); caxis([0,1]); title('corr(fus,A+)')
% subplot(1,3,3); imagesc(abs((V_dm_fus.')*VA_neg));  colorbar; colormap(gray); caxis([0,1]); title('corr(fus,A-)')


%% K-means clustering
% k = 7;
% 
% [X,Y] = meshgrid(1:100,1:100);
% 
% xy_coord = ones(numel(c_ldr),1);%1e-3*[X(:),Y(:)];%1e-2*c_ldr(:);%,5*];
% 
% idx_S  = kmeans([VS*diag(dS),mean(mean(VS*diag(dS)))*xy_coord],k);
% idx_A1 = kmeans([VA_pos*diag(dA_pos),mean(mean(VA_pos*diag(dA_pos)))*xy_coord],k);
% idx_A2 = kmeans([VA_neg*diag(dA_neg),mean(mean(VA_neg*diag(dA_neg)))*xy_coord],k);
% % idx_ad = kmeans([gindx_ad,xy_coord],k);
% 
% % idx_S  = kmeans(gindx_S,k);
% % idx_A1 = kmeans(gindx_A1,k);
% % idx_A2 = kmeans(gindx_A2,k);
% 
% % idx_S  = kmeans([VS*diag(dS),VA_pos*diag(dA_pos),VA_neg*diag(dA_neg)],k);
% 
% figure
% subplot(131); imagesc(reshape(idx_S ,size(c_ldr)))
% subplot(132); imagesc(reshape(idx_A1,size(c_ldr)))
% subplot(133); imagesc(reshape(idx_A2,size(c_ldr)))
% 
% 
% % Yo Dog K-means on K-means:
% idx_kmns = kmeans([idx_S,idx_A1,idx_A2],k);
% figure
% imagesc(reshape(idx_kmns,size(c_ldr)));

%%

params.k        = 20;
params.sigma2   = 0.5;
params.k_modes  = 7; % also 5 gives interesting results for A
params.siz      = size(c_ldr);
params.r_s      = 5;

params_fus = params;
params_fus.k_modes = 7;

k_dm = 100;

lambda    = norm(c_hps(:),'fro')/norm(c_ldr,'fro');
DIS_fus   = squareform(pdist([reshape(c_hps,numel(c_ldr),[]),lambda*c_ldr(:)]));
DIS_fus_l = (DIS_fus/std(DIS_fus(:)));

K_fus     = exp(-DIS_fus_l.^2/median(DIS_fus_l(:).^2));
K_fus     = diag(1./sum(K_fus,1)) * K_fus;

[V_dm_fus,D_dm_fus] = eigs(K_fus,neig);
DM_dis_fus          = squareform(pdist(V_dm_fus * (D_dm_fus.^1)));

DM_dis_al = squareform(pdist([VS * DS ,VA_pos * diag(dA_pos) ,VA_neg * diag(dA_neg)]));
DM_dis_S  = squareform(pdist(VS * DS));
DM_dis_Ap = squareform(pdist(VA_pos * diag(dA_pos)));
DM_dis_An = squareform(pdist(VA_neg * diag(dA_neg)));

[modes_ind_fus,y_lbld_fus]  = diffusion_clustering(DIS_fus_l,  DM_dis_fus,    DM_dis_fus, params_fus);
[modes_ind_al,  y_lbld_al]  = diffusion_clustering(DIS_fus_l,  DM_dis_al,     DM_dis_al,  params);
[modes_ind_S, y_lbld_S]     = diffusion_clustering(DIS_fus_l,  DM_dis_S,      DM_dis_S,   params);
[modes_ind_Ap, y_lbld_Ap]   = diffusion_clustering(DIS_fus_l,  DM_dis_Ap,     DM_dis_Ap,  params);
[modes_ind_An, y_lbld_An]   = diffusion_clustering(DIS_fus_l,  DM_dis_An,     DM_dis_An,  params);
[modes_ind_mx, y_lbld_mx]   = diffusion_clustering(DIS_fus_l,  DM_dis_S,      DM_dis_Ap,  params);

figure
subplot(4,3,[4,7]); imagesc(reshape(y_lbld_fus,params.siz));  title('Reference','Interpreter','latex'); 
subplot(4,3,[2,5]); imagesc(reshape(y_lbld_al,params.siz));   title('$$[V_{S},V_{A+},V_{A-}]$$','Interpreter','latex'); colorbar;
subplot(4,3,[3,6]); imagesc(reshape(y_lbld_S,params.siz));    title('$$V_{S}$$','Interpreter','latex'); colorbar;
subplot(4,3,[8,11]); imagesc(reshape(y_lbld_Ap,params.siz));   title('$$V_{A+}$$','Interpreter','latex'); colorbar;
subplot(4,3,[9,12]); imagesc(reshape(y_lbld_An,params.siz));   title('$$V_{A-}$$','Interpreter','latex'); colorbar;

figure
subplot(131); imagesc(reshape(y_lbld_fus,params.siz));  title('Reference','Interpreter','latex');
subplot(132); imagesc(reshape(y_lbld_Ap,params.siz));   title('$$V_{A+}$$','Interpreter','latex');
subplot(133); imagesc(reshape(y_lbld_mx,params.siz));   title('$$V_{A+}\ and\ V_{S}$$','Interpreter','latex');
