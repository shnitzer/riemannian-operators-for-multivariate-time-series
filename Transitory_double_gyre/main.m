% MAIN - Transitory Double Gyre Simulation
close all; clear; clc;

gendata  = 0; % 0 - load previously generated data,             1 - generate new data
procdata = 1; % 0 - use saved (previously analyzed) operators,  1 - process data (don't use saved results)

% Data parameters:
N     = 2500;      % Number of data points in each time frame (No. of trajectories)
dt    = 0.0039;    % Time difference between frames
ep    = 0.5;%1;0.1;  % Kernel scale
fname = ['Saved_eigs\correct_norm_N_',num2str(N),'_dt_',num2str(dt),'_ep_',num2str(ep)];


if procdata == 1

    if gendata == 0
        
        load(['gyre_transitory_data_N_',num2str(N),'_dt_',num2str(dt),'.mat']);
        
    elseif gendata == 1
        
        ptn_x    = sqrt(N);
        ptn_y    = sqrt(N);
        p.dt     = dt;
        p.t_fin  = 1;
        epsi     = 1e-2;
        
        [init_x, init_y] = meshgrid(linspace(epsi,1-epsi,ptn_x),linspace(epsi,1-epsi,ptn_y));
        init_x   = init_x(:);
        init_y   = init_y(:);
        t_vec    = 0:p.dt:p.t_fin;
        p.init_x = init_x;
        p.init_y = init_y;
        p.t_vec  = t_vec;
        xilen    = length(init_x);
        
        p.f1 = 2*pi;
        p.f2 = pi;
        p.A1 = 2;
        p.A2 = 10;
        psip_dx = @(t,x) p.A1*p.f1*cos(p.f1*x(1:xilen)).*sin(p.f2*x((xilen+1):end));
        psip_dy = @(t,x) p.A1*p.f2*sin(p.f1*x(1:xilen)).*cos(p.f2*x((xilen+1):end));
        psif_dx = @(t,x) p.A2*p.f2*cos(p.f2*x(1:xilen)).*sin(p.f1*x((xilen+1):end));
        psif_dy = @(t,x) p.A2*p.f1*sin(p.f2*x(1:xilen)).*cos(p.f1*x((xilen+1):end));
        s_t     = @(t) t.^2 .* (3-2*t);
        
        xd          = @(t,x) [-(1-s_t(t)).*psip_dy(t,x) - s_t(t).*psif_dy(t,x); (1-s_t(t)).*psip_dx(t,x) + s_t(t).*psif_dx(t,x)];
        [~,sig1]    = ode45(xd,t_vec,[init_x; init_y]);
        
        x   = sig1(:,1:length(init_x)).';
        y   = sig1(:,(length(init_x)+1):end).';
        
        save(['gyre_transitory_data_N_',num2str(size(x,1)),'_dt_',num2str(p.dt),'.mat'],'x','y','p');
        
    end
    
    %% Construct SPD kernels
    
    T       = length(p.t_vec);
    K       = cell(1,T);
    h       = waitbar(0,'Please wait while the SPD kernels are calculated');
    epparam = ep;
    
    for t = 1:T
        waitbar(t/T,h);
        
        data = [x(:,t),y(:,t)];
        
        K{t} = construct_spd_kern( data, epparam );
        
    end
    close(h)
    
    %% Construct and save wavelet kernels
    
    saveyn    = 1;              % If the kernel size and number of kernels is large, it is better to save the kernels for each level separately.
    
    sav.yn    = saveyn;
    sav.fname = fname;
    r         = 100;            % Chosen empirically
    nE        = 10;             % Chosen empirically
    nL        = floor(log2(T));
    
    riemann_wavelets( K, sav, r, nE, nL );
    
else
    load(['gyre_transitory_data_N_',num2str(N),'_dt_',num2str(dt),'.mat']);
    T = size(x,2);
end

%% Display results
%% Eigenvectors

nL     = floor(log2(T));
figpos = [100,200,1360,420];
axang  = [-8,50];
veind  = 3;%[1,2,2,3]; % 1
v.type = 'eig';

savgyre.yn = 0;
if ep == 0.5
    figttl     = @(cop,cveind,clev)['eps05\psi_',cop,cveind,'_lev',clev,'_t'];
elseif ep == 1
    figttl     = @(cop,cveind,clev)['eps1\psi_',cop,cveind,'_lev',clev,'_t'];
else
    figttl     = @(cop,cveind,clev)['psi_',cop,cveind,'_lev',clev,'_t'];
end
savptdisp  = 0;

for lev = 8
    load([fname,'_wavelet_eigs_level_',num2str(lev),'.mat'])
    
%     v.data     = VS_l;
%     v.vind     = veind+1;
%     thresh     = 0.1*N;%0.03;
%     ftitle     = ['$$\psi^{(\mathbf{S}_{',num2str(lev),'})}_{',num2str(veind+1),'}$$'];
%     savgyre.nm = figttl('S',num2str(veind+1),num2str(lev));
% %     display_gyre_results_analysis( x, y, v, lev, ftitle, figpos, axang, savgyre );
%     display_gyre_results( x, y, v, lev, ftitle, figpos, axang, savgyre );
%     point_displacement_analysis_gyre( x, y, v, thresh, lev, ftitle )
%     switch lev
%         case {7,8};     set(gcf,'Position',[100,200,300+(8-lev)*100,310]);
%         case {5,6};     set(gcf,'Position',[100,200,300+(8-lev)*150,310]);
%         case {1,2,3,4}; set(gcf,'Position',[100,200,300+(8-lev)*200,310]);
%     end
%     if savptdisp
%         if ep == 0.5
%             saveas(gcf,['figures\eps05\lev',num2str(lev),'_eign',num2str(veind+1),'_ptdist_S_var_LineWidth.fig'])
%             saveas(gcf,['figures\eps05\lev',num2str(lev),'_eign',num2str(veind+1),'_ptdist_S_var_LineWidth.png'])
%             saveas(gcf,['figures\eps05\lev',num2str(lev),'_eign',num2str(veind+1),'_ptdist_S_var_LineWidth'],'epsc')
%         elseif ep == 1
%             saveas(gcf,['figures\eps1\lev',num2str(lev),'_eign',num2str(veind+1),'_ptdist_S_var_LineWidth.fig'])
%             saveas(gcf,['figures\eps1\lev',num2str(lev),'_eign',num2str(veind+1),'_ptdist_S_var_LineWidth.png'])
%         else
%             saveas(gcf,['figures\lev',num2str(lev),'_ptdist_S_var_LineWidth.fig'])
%             saveas(gcf,['figures\lev',num2str(lev),'_ptdist_S_var_LineWidth.png'])
%         end
%     end
    
    
    v.data     = VA_l;
    v.vind     = veind;
    thresh     = 0.1*N; %0.03;
    ftitle     = ['$$\psi^{(\mathbf{A}_{',num2str(lev),'})}_{',num2str(veind),'}$$'];
    ftitle     = ['$$\psi^{(\mathbf{A}_{',num2str(lev),'})}_{',num2str(-1),'}$$'];
    savgyre.nm = figttl('A',num2str(veind),num2str(lev));
    display_gyre_results( x, y, v, lev, ftitle, figpos, axang, savgyre );
%     display_gyre_results_analysis( x, y, v, lev, ftitle, figpos, axang, savgyre );
    point_displacement_analysis_gyre( x, y, v, thresh, lev, ftitle )
    switch lev
        case {7,8};     set(gcf,'Position',[100,200,300+(8-lev)*100,310]);
        case {5,6};     set(gcf,'Position',[100,200,300+(8-lev)*150,310]);
        case {1,2,3,4}; set(gcf,'Position',[100,200,300+(8-lev)*200,310]);
    end
%     if savptdisp
%         if ep == 0.5
%             saveas(gcf,['figures\eps05\lev',num2str(lev),'_eign',num2str(veind),'_ptdist_A_var_LineWidth.fig'])
%             saveas(gcf,['figures\eps05\lev',num2str(lev),'_eign',num2str(veind),'_ptdist_A_var_LineWidth.png'])
%             saveas(gcf,['figures\eps05\lev',num2str(lev),'_eign',num2str(veind),'_ptdist_A_var_LineWidth'],'epsc')
%         elseif ep == 1
%             saveas(gcf,['figures\eps1\lev',num2str(lev),'_eign',num2str(veind),'_ptdist_A_var_LineWidth.fig'])
%             saveas(gcf,['figures\eps1\lev',num2str(lev),'_eign',num2str(veind),'_ptdist_A_var_LineWidth.png'])
%         else
%             saveas(gcf,['figures\lev',num2str(lev),'_eign',num2str(veind),'_ptdist_A_var_LineWidth.fig'])
%             saveas(gcf,['figures\lev',num2str(lev),'_eign',num2str(veind),'_ptdist_A_var_LineWidth.png'])
% %             saveas(gcf,['figures\lev',num2str(lev),'_ptdist_A_var_LineWidth.fig'])
% %             saveas(gcf,['figures\lev',num2str(lev),'_ptdist_A_var_LineWidth.png'])
%         end
%     end
end

%% Clustering
if 0
nKs     = 5; % Number of clusters for the eigenvectors of operator S
nKa     = 4; % Number of Clusters for the eigenvectors of operator A
serng   = 1:3; % Eigenvector range to use from S
aerng   = 1:3; % Eigenvector range to use from A
idx_s   = cell(1,nL); % Cluster index for all points based on S
idx_a   = cell(1,nL); % Cluster index for all points based on A

idx_s_all = cell(1,nL);
idx_a_all = cell(1,nL);

for lev = 1:nL
    
    load([fname,'_wavelet_eigs_level_',num2str(lev),'.mat'])
    
    len         = length(VS_l);
    idx_s{lev}  = zeros(N,len);
    idx_a{lev}  = zeros(N,len);
    
    for jj = 1:len
        if jj == 1
            [idx_s{lev}(:,jj),C_s,~,cd_S] = kmeans(VS_l{jj}(:,serng),nKs);
            [idx_a{lev}(:,jj),C_a,~,cd_A] = kmeans(VA_l{jj}(:,aerng),nKa);
            [~,cS_i] = min(cd_S,[],1);
            [~,cA_i] = min(cd_A,[],1);
        else
            idx_s{lev}(:,jj) = kmeans(VS_l{jj}(:,serng),nKs,'Start',VS_l{jj}(cS_i,serng));
            idx_a{lev}(:,jj) = kmeans(VA_l{jj}(:,aerng),nKa,'Start',VA_l{jj}(cA_i,aerng));
        end
    end
   
    idx_s_all{lev} = zeros(N,1);
    idx_a_all{lev} = zeros(N,1);
    tmp_S  = [];
    tmp_A  = [];
    for jj = 1:len
        tmp_S = [tmp_S,VS_l{jj}(:,serng)];
        tmp_A = [tmp_A,VA_l{jj}(:,serng)];
    end
    idx_s_all{lev} = kmeans(tmp_S,nKs);
    idx_a_all{lev} = kmeans(tmp_A,nKa);
end

%% Display clustering
ccolor = [0 0 1;...
          0 1 0;...
          1 0 0;...
          1 1 0;...
          0 1 1];

v.map  = ccolor;
v.type = 'clust';

figpos = [100,200,1360,420];
axang  = [-8,50];

for lev = 5
    v.data = idx_s{lev};
    ftitle = ['S clustering - level ',num2str(lev)];
    display_gyre_results( x, y, v, lev, ftitle, figpos, axang );
    
    v.data = idx_a{lev};
    ftitle = ['A clustering - level ',num2str(lev)];
    display_gyre_results( x, y, v, lev, ftitle, figpos, axang );
end

%% Display unified clustering for all operators in the same level

ccolor = [0 0 1;...
          0 1 0;...
          1 0 0;...
          1 1 0;...
          0 1 1];

v.map  = ccolor;
v.type = 'uclust';

figpos = [100,200,1360,420];
axang  = [-8,50];

for lev = 7
    v.data = idx_s_all{lev};
    ftitle = ['S unified clustering - level ',num2str(lev)];
    display_gyre_results( x, y, v, lev, ftitle, figpos, axang );
    
    v.data = idx_a_all{lev};
    ftitle = ['A unified clustering - level ',num2str(lev)];
    display_gyre_results( x, y, v, lev, ftitle, figpos, axang );
end
end
%%
% %% Display clustering - with consecutive kernels - not correct!
% % since the consecutive frame is already taken into account in the original
% % display function
% 
% ccolor = [0 0 1;...
%           0 1 0;...
%           1 0 0;...
%           1 1 0;...
%           0 1 1];
% 
% v.map  = ccolor;
% v.type = 'clust';
% 
% figpos = [100,200,1360,420];
% axang  = [-8,50];
% 
% for lev = 5
%     v.data = idx_s{lev};
%     ftitle = ['S clustering - level ',num2str(lev)];
%     display_gyre_results_consec( x, y, v, lev, ftitle, figpos, axang );
%     
%     v.data = idx_a{lev};
%     ftitle = ['A clustering - level ',num2str(lev)];
%     display_gyre_results_consec( x, y, v, lev, ftitle, figpos, axang );
% end
% 
% %% Display eigenvectors - with consecutive kernels - not correct!
% % since the consecutive frame is already taken into account in the original
% % display function
% 
% nL     = floor(log2(T));
% figpos = [100,200,1360,420];
% axang  = [-8,50];
% veind  = 1;
% v.type = 'eig';
% 
% for lev = 6
%     load([fname,'_wavelet_eigs_level_',num2str(lev),'.mat'])
%     
%     v.data = VS_l;
%     v.vind = veind;
%     ftitle = ['S eigenvector no. ',num2str(veind),', level ',num2str(lev)];
%     display_gyre_results_consec( x, y, v, lev, ftitle, figpos, axang );
%     
%     v.data = VA_l;
%     ftitle = ['A eigenvector no. ',num2str(veind),', level ',num2str(lev)];
%     display_gyre_results_consec( x, y, v, lev, ftitle, figpos, axang );
% end