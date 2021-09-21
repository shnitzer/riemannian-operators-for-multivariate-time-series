clear; clc; close all;

%% Full rank

v1   = [1;1;1;1]; v2   = [1;-1;1;-1]; v3   = [-1;1;1;-1]; v4   = [1;1;-1;-1];
e(1) = 0.5;       e(2) = 0.2;         e(3) = 0.01;        e(4) = 1;
d(1) = 0.01;      d(2) = 0.2;         d(3) = 0.5;         d(4) = 1;

V  = [v1/norm(v1),v2/norm(v2),v3/norm(v3),v4/norm(v4)];

Vplt = V;
i_sr = 1:4;

colormat = [0,0,1;0,1,1;0,1,0.5;0,1,0];
figpos   = [675,220,560,530];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ma = V*diag(e)*(V.');
Mb = V*diag(d)*(V.');

[Va,Da] = eig(Ma);
[Vb,Db] = eig(Mb);

d1_diag = diag(Da);
d2_diag = diag(Db);

S = sqrtm(Ma) * sqrtm( inv(sqrtm(Ma)) * Mb * inv(sqrtm(Ma))) * sqrtm(Ma);
A = sqrtm(S)  * logm(  inv(sqrtm(S))  * Ma * inv(sqrtm(S)) ) * sqrtm(S);

[VS,DS] = eig(S);
[VA,DA] = eig(A);

ds_diag  = diag(DS);
da_diag  = diag(DA);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% sorting the eigenvectors %%%%%%%%%%%

[~,si_a] = sort(sum(abs(diff(Va))),'ascend');
[~,si_b] = sort(sum(abs(diff(Vb))),'ascend');
[~,si_S] = sort(sum(abs(diff(VS))),'ascend');

si_A     = zeros(1,4);
for ii = 1:4
    [~,si_A(ii)] = max(abs(Va(:,si_a(ii)).'*VA));
end

% [~,si_A] = sort(sum(abs(diff(VA))),'ascend');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% plot eigenvectors & eigenvalues %%%%%%%%

% 
figure('name','full rank all');
for ii = 1:4
    all_evals = [d1_diag(si_a(ii)),d2_diag(si_b(ii)),ds_diag(si_S(ii)),da_diag(si_A(ii))];
    all_color = [1,0.7,1;1,0,1;0,0,1;0,1,1];
    subplot(4,4,((ii-1)*4+1):((ii-1)*4+2)); br = bar(1:4,all_evals,'FaceColor','flat','EdgeColor','k');
    for jj = 1:4
        br.CData(jj,:) = all_color(jj,:);
        tx = text(jj,max(all_evals(jj),0),num2str(all_evals(jj),3),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom','Color','k'); 
    end
    ylabel(['$$\lambda^{(\ell)}_',num2str(ii),'$$'],'Interpreter','latex');
    xlim([0.5,4.5]);
    ybds = ylim;
    ylim([ybds(1)-0.1,ybds(2)+0.6]);
    if ii == 4
        set(gca,'FontSize',12); set(gca,'xtick',1:4,'xticklabels',{'\lambda^{(1)}_n','\lambda^{(2)}_n','\lambda^{(S)}_n','\lambda^{(A)}_n'}); %      xlim([0.5,4.5]);    ylim([-0.01,2.01]);
    else
        set(gca,'FontSize',12); set(gca,'xticklabels',[]);
    end
    if ii == 1
        title('Eigenvalues');
    end
        
    subplot(4,4,((ii-1)*4+3):(ii*4)); plot(1:4,sign(Va(1,si_a(ii)))*Va(:,si_a(ii)).','Color',all_color(1,:),'LineWidth',3);                     hold on;
    subplot(4,4,((ii-1)*4+3):(ii*4)); plot(1:4,sign(Vb(1,si_b(ii)))*Vb(:,si_b(ii)).','Color',all_color(2,:),'LineWidth',3); 
    subplot(4,4,((ii-1)*4+3):(ii*4)); plot(1:4,sign(VS(1,si_S(ii)))*VS(:,si_S(ii)).','Color',all_color(3,:),'LineWidth',2,'LineStyle',':'); 
    sb_c = subplot(4,4,((ii-1)*4+3):(ii*4)); plot(1:4,sign(VA(1,si_A(ii)))*VA(:,si_A(ii)).','Color',all_color(4,:),'LineWidth',2,'LineStyle','--');    hold off;
    sb_pos = get(sb_c,'position'); 
    xlim([0.8,4.2]);
    if ii == 4
        set(gca,'xtick',1:4,'xticklabels',{'\psi_n[1]','\psi_n[2]','\psi_n[3]','\psi_n[4]'},'FontSize',12);
        ylim([-0.8,0.8])
    else
        set(gca,'xticklabels',[],'FontSize',12); 
        ylim([-0.7,0.7]);
    end
    annotation('textbox', [sb_pos(1)+0.985*sb_pos(3), sb_pos(2)+0.8*sb_pos(4), 0, 0], 'string', ['\psi_',num2str(ii)],'FontSize',12)
    if ii == 1
        ylim([0,1]);
        title('Eigenvectors')
        h_lg = legend('M_1','M_2','S','A');
        set(h_lg,'Position',[0.932,0.78,0.06,0.2]);
    end
end

set(gcf,'Position',[170,240,920,420]);

% saveas(gcf,'figures for paper\full_rank.fig');
% saveas(gcf,'figures for paper\full_rank.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
figure('name','full rank expected eigenvalues');
for ii = 1:4
    l1 = abs(d1_diag(si_a(ii))); l2 = abs(d2_diag(si_b(ii)));
    all_evals = [sqrt(l1*l2),0.5*sqrt(l1*l2)*log(l1/l2)];
    all_color = [0,0,1;0,1,1];
    subplot(4,1,ii); br = bar(1:2,all_evals,'FaceColor','flat','EdgeColor','k');
    for jj = 1:2
        br.CData(jj,:) = all_color(jj,:);
        tx = text(jj,max(all_evals(jj),0),num2str(all_evals(jj),3),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom','Color','k'); 
    end
    ylabel(['$$\lambda^{(\ell)}_',num2str(ii),'$$'],'Interpreter','latex');
    xlim([0.5,2.5]);
    ybds = ylim;
    ylim([ybds(1)-0.1,ybds(2)+0.6]);
    if ii == 4
        set(gca,'FontSize',12); set(gca,'xtick',1:4,'xticklabels',{'\lambda^{(S)}_n','\lambda^{(A)}_n'}); %      xlim([0.5,4.5]);    ylim([-0.01,2.01]);
    else
        set(gca,'FontSize',12); set(gca,'xticklabels',[]);
    end
    if ii == 1
        title('Expected Eigenvalues');
    end
    
end
set(gcf,'Position',[700,40,260,420]);

% saveas(gcf,'figures for paper\full_rank_expected_eval.fig');
% saveas(gcf,'figures for paper\full_rank_expected_eval.png');


%% Not full rank

v1   = [1;1;1;1]; v2   = [1;-1;1;-1]; v3   = [-1;1;1;-1]; v4   = [1;1;-1;-1];
e(1) = 0.5;       e(2) = 0;           e(3) = 0.01;        e(4) = 1;
d(1) = 0.01;      d(2) = 0;           d(3) = 0.5;         d(4) = 1;

V  = [v1/norm(v1),v2/norm(v2),v3/norm(v3),v4/norm(v4)];

Vplt = V;
i_sr = 1:4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ma = V*diag(e)*(V.');
Mb = V*diag(d)*(V.');

[Va,Da] = eig(Ma);
[Vb,Db] = eig(Mb);

d1_diag = diag(Da);
d2_diag = diag(Db);

Slr = FixedGeodes_eff( Mb,Ma,0.5,3 );
Alr = FixedGeodes_eff_proj( Slr,Ma,3 );

[VS1,DS1] = eig(Slr);
[VA1,DA1] = eig(Alr);

ds_diag  = diag(DS1);
da_diag  = diag(DA1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% sorting the eigenvectors %%%%%%%%%%%

[~,si_a] = sort(sum(abs(diff(Va))),'ascend');
[~,si_b] = sort(sum(abs(diff(Vb))),'ascend');
[~,si_S] = sort(sum(abs(diff(VS1))),'ascend');

si_A     = zeros(1,4);
for ii = 1:4
    [~,si_A(ii)] = max(abs(Va(:,si_a(ii)).'*VA1));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% plot eigenvectors & eigenvalues %%%%%%%%

figure('name','rank 3 all');
for ii = 1:4
    all_evals = [d1_diag(si_a(ii)),d2_diag(si_b(ii)),ds_diag(si_S(ii)),da_diag(si_A(ii))];
    all_color = [1,0.7,1;1,0,1;0,0,1;0,1,1];
    subplot(4,4,((ii-1)*4+1):((ii-1)*4+2)); br = bar(1:4,all_evals,'FaceColor','flat','EdgeColor','k');
    for jj = 1:4
        br.CData(jj,:) = all_color(jj,:);
        tx = text(jj,max(all_evals(jj),0),num2str(all_evals(jj),3),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom','Color','k'); 
    end
    ylabel(['$$\lambda^{(\ell)}_',num2str(ii),'$$'],'Interpreter','latex');
    xlim([0.5,4.5]);
    ybds = ylim;
    ylim([ybds(1)-0.1,ybds(2)+0.6]);
    if ii == 4
        set(gca,'FontSize',12); set(gca,'xtick',1:4,'xticklabels',{'\lambda^{(1)}_n','\lambda^{(2)}_n','\lambda^{(S)}_n','\lambda^{(A)}_n'}); %      xlim([0.5,4.5]);    ylim([-0.01,2.01]);
    else
        set(gca,'FontSize',12); set(gca,'xticklabels',[]);
    end
    if ii == 1
        title('Eigenvalues');
    end
        
    subplot(4,4,((ii-1)*4+3):(ii*4)); plot(1:4,sign(Va(1,si_a(ii)))*Va(:,si_a(ii)).','Color',all_color(1,:),'LineWidth',3);                     hold on;
    subplot(4,4,((ii-1)*4+3):(ii*4)); plot(1:4,sign(Vb(1,si_b(ii)))*Vb(:,si_b(ii)).','Color',all_color(2,:),'LineWidth',3); 
    subplot(4,4,((ii-1)*4+3):(ii*4)); plot(1:4,sign(VS1(1,si_S(ii)))*VS1(:,si_S(ii)).','Color',all_color(3,:),'LineWidth',2,'LineStyle',':'); 
    sb_c = subplot(4,4,((ii-1)*4+3):(ii*4)); plot(1:4,sign(VA1(1,si_A(ii)))*VA1(:,si_A(ii)).','Color',all_color(4,:),'LineWidth',2,'LineStyle','--');    hold off;
    sb_pos = get(sb_c,'position'); 
    xlim([0.8,4.2]);
    if ii == 4
        set(gca,'xtick',1:4,'xticklabels',{'\psi_n[1]','\psi_n[2]','\psi_n[3]','\psi_n[4]'},'FontSize',12);
        ylim([-0.8,0.8])
    else
        set(gca,'xticklabels',[],'FontSize',12); 
        ylim([-0.7,0.7]);
    end
    annotation('textbox', [sb_pos(1)+0.985*sb_pos(3), sb_pos(2)+0.8*sb_pos(4), 0, 0], 'string', ['\psi_',num2str(ii)],'FontSize',12)
    if ii == 1
        ylim([0,1]);
        title('Eigenvectors')
        h_lg = legend('M_1','M_2','S','A');
        set(h_lg,'Position',[0.932,0.78,0.06,0.2]);
    end
end

set(gcf,'Position',[170,240,920,420]);

% saveas(gcf,'figures for paper\rank_3.fig');
% saveas(gcf,'figures for paper\rank_3.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Expected eigenvalues

figure('name','rank 3 expected eigenvalues');
for ii = 1:4
    l1 = abs(d1_diag(si_a(ii))); l2 = abs(d2_diag(si_b(ii)));
    all_evals = [sqrt(l1*l2),0.5*sqrt(l1*l2)*log(l1/l2)];
    all_color = [0,0,1;0,1,1];
    subplot(4,1,ii); br = bar(1:2,all_evals,'FaceColor','flat','EdgeColor','k');
    for jj = 1:2
        br.CData(jj,:) = all_color(jj,:);
        tx = text(jj,max(all_evals(jj),0),num2str(all_evals(jj),3),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom','Color','k'); 
    end
    ylabel(['$$\lambda^{(\ell)}_',num2str(ii),'$$'],'Interpreter','latex');
    xlim([0.5,2.5]);
    ybds = ylim;
    ylim([ybds(1)-0.1,ybds(2)+0.6]);
    if ii == 4
        set(gca,'FontSize',12); set(gca,'xtick',1:4,'xticklabels',{'\lambda^{(S)}_n','\lambda^{(A)}_n'}); %      xlim([0.5,4.5]);    ylim([-0.01,2.01]);
    else
        set(gca,'FontSize',12); set(gca,'xticklabels',[]);
    end
    if ii == 1
        title('Expected Eigenvalues');
    end
    
end
set(gcf,'Position',[700,40,260,420]);

% saveas(gcf,'figures for paper\rank_3_expected_eval.fig');
% saveas(gcf,'figures for paper\rank_3_expected_eval.png');