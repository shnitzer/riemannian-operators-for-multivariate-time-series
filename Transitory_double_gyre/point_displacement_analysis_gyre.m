function point_displacement_analysis_gyre(x,y,v,thresh,lev, ftitle)

sbplt    = 1;

fnum     = 8;
T        = size(x,2);
len      = floor(T/(2^lev));

jump     = max([floor(2^lev/fnum),1]);
t        = 1:(len*2^lev);
trng_all = reshape(t,[],len).';
trng_all = trng_all(:,1:jump:end);

cfll = [0,0,1;   0,1,0;   0.5,0.5,0.5; 1,0,1    ];
cplt = [0,0,0.5; 0,0.5,0; 0,0,0;       0.5,0,0.5];
lplt = {'-','--','-.',':'};

if sbplt; figure('Name',ftitle); end
for ii = 1:size(trng_all,1)
    
    trng     = trng_all(ii,:);
    opidx    = ceil(trng(end)/2^lev);
    
%     pt_idx   = abs(v.data{opidx}(:,v.vind))>thresh;
    if length(v.vind) == size(trng_all,1)
        [~,pt_idx] = maxk(abs(v.data{opidx}(:,v.vind(ii))),thresh); 
    else
        [~,pt_idx] = maxk(abs(v.data{opidx}(:,v.vind)),thresh); 
    end
    
    x_dat    = x(pt_idx,trng);
    y_dat    = y(pt_idx,trng);
    dist_dat = sqrt((x_dat(:,1:end)-x_dat(:,1)).^2 + (y_dat(:,1:end)-y_dat(:,1)).^2);
    
    dist_dat_max = sqrt((x(:,trng)-x(:,trng(1))).^2 + (y(:,trng)-y(:,trng(1))).^2);
    dist_dat_max = max(dist_dat_max,[],1);
        
    kclust   = 3;
    clustidx = kmeans(dist_dat,kclust);
    
    if sbplt; %subplot(1,size(trng_all,1),ii);
    else; figure('Name',ftitle); end
    pltx = [trng, fliplr(trng)];
        
    dist_avg = zeros(kclust,length(trng));
    dist_std = zeros(kclust,length(trng));
    for ci   = 1:kclust
        dist_avg(ci,:) = mean(dist_dat(clustidx==ci,:),1);
        dist_std(ci,:) = std(dist_dat(clustidx==ci,:),[],1);
        
        pltfill = [dist_avg(ci,:)+dist_std(ci,:),fliplr(dist_avg(ci,:)-dist_std(ci,:))];
        fill(pltx, pltfill, cfll(ci,:),'FaceAlpha',0.1,'EdgeAlpha',0.15); hold on;
        plot(trng, dist_avg(ci,:),'Color',cplt(ci,:),'LineWidth',2*kclust*sum(clustidx==ci)/thresh,'LineStyle',lplt{ci});
    end
%     plot(trng,dist_dat_max,'Color',[0.1,0.1,0.1],'LineStyle',':')
%     xlim([trng(1),trng(end)]); 
    xlim([trng_all(1),trng_all(end)]); 
    ylim([0,sqrt(2)]);
    xlabel('t');
    if sbplt; if ii == 1; ylabel('Distance from origin point'); end
    else; ylabel('Distance from origin point'); end

%     text(trng(1),sqrt(2)-0.1,[num2str(sum(pt_idx)),' total pts']);
end


% y = rand(1,10); % your mean vector;
% x = 1:numel(y);
% std_dev = 1;
% curve1 = y + std_dev;
% curve2 = y - std_dev;
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% fill(x2, inBetween, 'g');
% hold on;
% plot(x, y, 'r', 'LineWidth', 2);