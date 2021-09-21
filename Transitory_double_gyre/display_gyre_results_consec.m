function display_gyre_results_consec( x, y, v, lev, ftitle, figpos, axang )

if nargin<6
    figpos = [370,445,1360,420];
    axang  = [-8,50];
end

fnum     = 8;
T        = size(x,2);
len      = floor(T/(2^lev));

if strcmp(v.type,'clust') || strcmp(v.type,'eig')
    jump     = max([floor(2^lev/fnum),1]);
    t        = 1:(len*2^lev);
    trng_all = reshape(t,[],len).';
    trng_all = trng_all(:,1:jump:end);
elseif strcmp(v.type,'uclust')
    jump     = 8;    
    t        = 1:(len*2^lev);
    trng_all = reshape(t,[],4).';
    trng_all = trng_all(:,1:jump:end);
end



for ii = 1:(size(trng_all,1)-1)
    
    trng     = trng_all(ii,:);
    trng_nxt = trng_all(ii+1,:);
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
    x_dat     = x(:,trng);
    y_dat     = y(:,trng);
    t_dat     = repmat(trng,[size(x,1),1]);
    x_dat_nxt = x(:,trng_nxt);
    y_dat_nxt = y(:,trng_nxt);
    t_dat_nxt = repmat(trng_nxt,[size(x,1),1]);    
    
    figure; set(gcf,'Position',figpos);
    subplot(2,1,1); scatter3(t_dat(:),x_dat(:),y_dat(:),10,repmat(mcolor,[length(trng),1]),'fill');
    set(gca,'view',axang);
    xlabel('t'); ylabel('x'); zlabel('y'); xlim([trng(1)-1,trng(end)+1]);
    colorbar;
    subplot(2,1,2); scatter3(t_dat_nxt(:),x_dat_nxt(:),y_dat_nxt(:),10,repmat(mcolor,[length(trng_nxt),1]),'fill');
    set(gca,'view',axang);
    xlabel('t'); ylabel('x'); zlabel('y'); xlim([trng_nxt(1)-1,trng_nxt(end)+1]);
    colorbar;
    suptitle(ftitle);
end

