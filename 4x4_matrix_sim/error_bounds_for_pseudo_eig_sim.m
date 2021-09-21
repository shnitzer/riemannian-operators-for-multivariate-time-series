clear; clc; %close all;

%% Full rank

v1   = [1;1;1;1]; v2   = [1;-1;1;-1]; v3   = [-1;1;1;-1]; v4   = [1;1;-1;-1];
e(1) = 0.5;       e(2) = 0.2;         e(3) = 1e-2;        e(4) = 1;
d(1) = 1e-2;      d(2) = 0.2;         d(3) = 0.5;         d(4) = 1;

V  = [v1/norm(v1),v2/norm(v2),v3/norm(v3),v4/norm(v4)];

epsi = 1e-4:1e-4:8e-3;
Nr  = 100;

err_S = zeros(length(epsi),Nr,4);
err_A = zeros(length(epsi),Nr,4);
err_S_ref = zeros(length(epsi),Nr,4);
err_A_ref = zeros(length(epsi),Nr,4);


for ei = 1:length(epsi)
    epsi_c = epsi(ei);
    
    for jj = 1:Nr
    
        V_tmp = zeros(4);
        
%         loopin = 1;
%         while loopin
        for ii = 1:4
            err_vec = randn(4,1);
            err_vec = (epsi_c*sqrt(d(ii)/e(ii))/max(abs(d-d(ii))))*err_vec./norm(err_vec);
%             err_vec = (epsi_c)*err_vec./norm(err_vec);
            V_tmp(:,ii) = V(:,ii) + err_vec;
            V_tmp(:,ii) = V_tmp(:,ii)/norm(V_tmp(:,ii));
        end
        
%         Mb = V_tmp*diag(d)*(V_tmp.');
%         db = eig(Mb);
%         
%         testvec = zeros(1,4);
%         for ii = 1:4
%             testvec(ii) = (epsi_c*sqrt(db(ii)/e(ii))/max(db(db~=db(ii)))) <= (epsi_c*sqrt(d(ii)/e(ii))/max(d(d~=d(ii))));
%         end
%         if all(testvec)
%             loopin = 0;
%         end
%         end
% V_tmp = V + epsi_c*randn(4,4);%[1;2;3;4]/norm([1;2;3;4]);
% V_tmp = V;
% V_tmp(:,4) = V(:,4) + epsi_c*randn(4,1);%[1;2;3;4]/norm([1;2;3;4]);
% V_tmp(:,4) = V_tmp(:,4)/norm(V_tmp(:,4));
%%%%%%%%%%%%%%% plot eigs:

% Vplt = V;
% i_sr = 1:4;

% figure('name','full rank');
% ax(1) = subplot(5,4,1:4);   sc = scatter(1:4,e,100,[0,0,1;0,1,1;0,1,0.5;0,1,0],'fill');     xlim([0.5,4.5]);    ylim([-0.01,2.01]);
% set(gca,'FontSize',12); ylabel('\lambda_a');  set(gca,'xticklabels',[]); for ii = 1:4; datatip(sc,ii,e(ii)); end
% ax(2) = subplot(5,4,5:8);   sc = scatter(1:4,d,100,[0,0,1;0,1,1;0,1,0.5;0,1,0],'fill');     xlim([0.5,4.5]);    ylim([-0.01,2.01]);
% set(gca,'FontSize',12); ylabel('\lambda_b');  set(gca,'xticklabels',[]); for ii = 1:4; datatip(sc,ii,d(ii)); end
% ax(3) = subplot(5,4,[9,13,17]);  scatter(Vplt(:,i_sr(1)),1:4,100,[0,0,1],'fill');           xlim([-0.75,0.75]); ylim([0,5]);
% set(gca,'FontSize',12); xlabel('v_1');      set(gca,'yticklabels',[]);
% ax(4) = subplot(5,4,[10,14,18]); scatter(Vplt(:,i_sr(2)),1:4,100,[0,1,1],'fill');           xlim([-0.75,0.75]); ylim([0,5]);
% set(gca,'FontSize',12); xlabel('v_2');      set(gca,'yticklabels',[]);
% ax(5) = subplot(5,4,[11,15,19]); scatter(Vplt(:,i_sr(3)),1:4,100,[0,1,0.5],'fill');         xlim([-0.75,0.75]); ylim([0,5]);
% set(gca,'FontSize',12); xlabel('v_3');      set(gca,'yticklabels',[]);
% ax(6) = subplot(5,4,[12,16,20]); scatter(Vplt(:,i_sr(4)),1:4,100,[0,1,0],'fill');           xlim([-0.75,0.75]); ylim([0,5]);
% set(gca,'FontSize',12); xlabel('v_4');      set(gca,'yticklabels',[]);
% 
% set(gcf,'Position',[675,420,560,530])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ma = V*diag(e)*(V.');
Mb = V_tmp*diag(d)*(V_tmp.');

[Vab,Dab] = eig(Ma,Mb*Ma);

[Va,Da] = eig(Ma);
[Vb,Db] = eig(Mb);

db_diag = diag(Db);
db_diag = sort(db_diag,'descend');
[~,srI] = sort(d,'descend');
db_diag = db_diag(srI);

S = sqrtm(Ma) * sqrtm( inv(sqrtm(Ma)) * Mb * inv(sqrtm(Ma))) * sqrtm(Ma);
A = sqrtm(S)  * logm(  inv(sqrtm(S))  * Ma * inv(sqrtm(S)) ) * sqrtm(S);

[VS,DS] = eig(S);
[VA,DA] = eig(A);

for ii = 1:4
    err_S(ei,jj,ii) = norm((S-sqrt(e(ii)*db_diag(ii))*eye(4))*V(:,ii));
    err_A(ei,jj,ii) = norm((A-0.5*sqrt(e(ii)*db_diag(ii))*log(e(ii)/db_diag(ii))*eye(4))*V(:,ii));
%     err_S_ref(ei,jj,ii) = norm(S*V(:,ii));
%     err_A_ref(ei,jj,ii) = norm(A*V(:,ii));
end

    end
end

figure
for ii = 1:4
% subplot(2,4,ii); errorbar(epsi,mean(err_S(:,:,ii),2),std(err_S(:,:,ii),[],2));
% subplot(2,4,ii); errorbar(epsi,mean(err_S(:,:,ii),2),max(err_S(:,:,ii),[],2));
subplot(2,4,ii); plot(epsi,mean(err_S(:,:,ii),2),'b');
hold on;
lambda2t_max = max(abs(d-d(ii)));
eps_S_bound = 1;%/(sqrt(db_diag(ii)/e(ii))/max(abs(db_diag-db_diag(ii))));%sqrt(e(ii)/d(ii))*lambda2t_max;
% eps_S_bound = max(1/(sqrt(db_diag(ii)/e(ii))/max(abs(db_diag-db_diag(ii)))),1/(sqrt(e(ii)/db_diag(ii))/max(abs(e-e(ii)))));
subplot(2,4,ii); plot(epsi,epsi*eps_S_bound,'r');
title(['Bound for S for \lambda_1=',num2str(e(ii)),', \lambda_2=',num2str(d(ii))]); xlabel('\epsilon');

subplot(2,4,ii+4); plot(epsi,mean(err_A(:,:,ii),2),'b');
hold on;
% eps_A_bound = 2/(d(ii)/(min(d)*sqrt(min(e)*e(ii)))+max([abs(log(min(e))),abs(log(min(d)))]));
eps_A_bound = 0.5*(d(ii)/(min(d)*sqrt(min(e)*e(ii)))+max([abs(log(min(e))),abs(log(min(d)))]));
% eps_A_bound = 1;%0.5*(db_diag(ii)/(sqrt(min(db_diag))*sqrt(min(e)*e(ii))));
subplot(2,4,ii+4); plot(epsi,epsi*eps_A_bound*eps_S_bound,'r');
title(['Bound for A for \lambda_1=',num2str(e(ii)),', \lambda_2=',num2str(db_diag(ii))]); xlabel('\epsilon');

end

figure
for ii = 1:4
    plot(epsi,mean(err_S(:,:,ii),2)); hold on;
end
ylabel('$$\left\Vert (\mathbf{S}_R-\sqrt{\lambda^{(1)}_i\lambda^{(2)}_i}\mathrm{I})\psi^{(1)}_i\right\Vert_2$$','Interpreter','Latex');
xlabel('$$\epsilon$$','Interpreter','Latex');
plot(epsi,epsi*eps_S_bound,'r'); hold off;
legend({'i=1','i=2','i=3','i=4','$$\epsilon$$'},'Interpreter','Latex');
set(gca,'FontSize',11);
