function display_gyre_results( x, y, v, lev, ftitle, figpos, axang, savfig )

if nargin<8
   if nargin<6
       figpos = [370,445,1360,420];
       axang  = [-8,50];
   end
   savfig.yn  = 0;
end

fnum     = 8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4
T        = size(x,2);
len      = floor(T/(2^lev));


if strcmp(v.type,'clust') || strcmp(v.type,'eig')
    jump     = max([floor(2^lev/fnum),1]);
    t        = 1:(len*2^lev);
    trng_all = reshape(t,[],len).';
    trng_all = trng_all(:,1:jump:end);
%     trng_all = trng_all(15:16,:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% grayout
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
        if length(v.vind) == size(trng_all,1)
            mcolor = v.data{opidx}(:,v.vind(ii));%repmat(v.eig{opidx}(:,v.vind),[1,length(trng)]);
        else
            mcolor = v.data{opidx}(:,v.vind);
        end
    end
    x_dat = x(:,trng);
    y_dat = y(:,trng);
    t_dat = repmat(trng,[size(x,1),1]);
    figure; set(gcf,'Position',figpos);
    scatter3(t_dat(:),x_dat(:),y_dat(:),10,repmat(mcolor,[length(trng),1]),'fill');
    set(gca,'view',axang);
    xlabel('t','FontSize',14); ylabel('x','FontSize',14); zlabel('y','FontSize',14); xlim([trng(1)-1,trng(end)+1]);
    colorbar;
    title(ftitle,'Interpreter','Latex','FontSize',16);
    caxis([-0.04,0.04]);
    if savfig.yn
        saveas(gcf,['figures\eigs\',savfig.nm,num2str(ii),'.fig']);
        saveas(gcf,['figures\eigs\',savfig.nm,num2str(ii),'.png']);
        saveas(gcf,['figures\eigs\',savfig.nm,num2str(ii)],'epsc');
    end
end

