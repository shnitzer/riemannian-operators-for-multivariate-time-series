function display_gyre_results_analysis( x, y, v, lev, ftitle, figpos, axang, savfig )

if nargin<8
   if nargin<6
       figpos = [370,445,1360,420];
       axang  = [-8,50];
   end
   savfig.yn  = 0;
end

fnum     = 4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4
T        = size(x,2);
len      = floor(T/(2^lev));
axang    = axang - [8,0];
tinds    = [5,6]; % [15,16]; % [13,14]; % [3,4];%

if strcmp(v.type,'clust') || strcmp(v.type,'eig')
    jump     = max([floor(2^lev/fnum),1]);
    t        = 1:(len*2^lev);
    trng_all = reshape(t,[],len).';
    trng_all = trng_all(:,1:jump:end);
    trng_all = trng_all(tinds,:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% grayout
elseif strcmp(v.type,'uclust')
    jump     = 8;    
    t        = 1:(len*2^lev);
    trng_all = reshape(t,[],4).';
    trng_all = trng_all(:,1:jump:end);
end



for ii = 1:size(trng_all,1)
    
    trng     = trng_all(ii,:);
    if strcmp(v.type,'clust') || strcmp(v.type,'eig')
        opidx = ceil(trng(end)/2^lev);
    elseif strcmp(v.type,'uclust')
        opidx = 1;
    end
    
    if strcmp(v.type,'clust') || strcmp(v.type,'uclust')
        mcolor = zeros(size(x,1),size(v.map,2));
        for ci = 1:size(v.map,1)
            mcolor(v.data(:,opidx)==ci,:) = repmat(v.map(ci,:),[sum(v.data(:,opidx)==ci),1]);
        end
%         mcolor = repmat(mcolor,[length(trng),1]);
    elseif strcmp(v.type,'eig')
        mcolor = v.data{opidx}(:,v.vind);%repmat(v.eig{opidx}(:,v.vind),[1,length(trng)]);
    end
    
    ftitle_cut = [ftitle(1:23),'^{(',num2str(tinds(ii)),')}',ftitle(24:end)];
    
    x_dat = x(:,trng);
    y_dat = y(:,trng);
    t_dat = repmat(trng,[size(x,1),1]);
    figure; set(gcf,'Position',figpos);
    if ii<size(trng_all,1)
        sp1 = subplot(1,2,1);
        scatter3(t_dat(:),x_dat(:),y_dat(:),10,repmat(mcolor,[length(trng),1]),'fill');
        sp1.Position = sp1.Position + [-0.05, 0, 0.06, -0.02];
    elseif ii>1
        sp1 = subplot(1,2,2);
        scatter3(t_dat(:),x_dat(:),y_dat(:),10,repmat(mcolor,[length(trng),1]),'fill');
        colorbar;
        sp1.Position = sp1.Position + [-0.05, 0, 0.075, -0.02];
    end
    xlabel('t','FontSize',14); ylabel('x','FontSize',14); zlabel('y','FontSize',14); xlim([trng(1)-1,trng(end)+1]);
    title(ftitle_cut,'Interpreter','Latex','FontSize',16);
    caxis([-0.04,0.04]);
    set(gca,'view',axang);
        
    if ii<size(trng_all,1)
        trng_nxt = trng_all(ii+1,:);
        x_dat = x(:,trng_nxt);
        y_dat = y(:,trng_nxt);
        t_dat = repmat(trng_nxt,[size(x,1),1]);
    elseif ii>1
        trng_nxt = trng_all(ii-1,:);
        x_dat = x(:,trng_nxt);
        y_dat = y(:,trng_nxt);
        t_dat = repmat(trng_nxt,[size(x,1),1]);
    end
    
    if ii<size(trng_all,1)
        sp2 = subplot(1,2,2);
        scatter3(t_dat(:),x_dat(:),y_dat(:),10,repmat(mcolor,[length(trng_nxt),1]),'fill');
        colorbar;
        sp2.Position = sp2.Position + [-0.05, 0, 0.075, -0.02];
    elseif ii>1
        sp2 = subplot(1,2,1);
        scatter3(t_dat(:),x_dat(:),y_dat(:),10,repmat(mcolor,[length(trng_nxt),1]),'fill');
        sp2.Position = sp2.Position + [-0.05, 0, 0.06, -0.02];
    end
    xlabel('t','FontSize',14); ylabel('x','FontSize',14); zlabel('y','FontSize',14); xlim([trng_nxt(1)-1,trng_nxt(end)+1]);
    title(ftitle_cut,'Interpreter','Latex','FontSize',16);
    caxis([-0.04,0.04]);
    set(gca,'view',axang);
        
end

