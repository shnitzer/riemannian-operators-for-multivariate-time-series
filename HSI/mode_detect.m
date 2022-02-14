function [p,modes_ind] = mode_detect(DM_dist,Euc_dist,siz,k,sigma2,k_modes)

% Implementation of "Diffusion geometric methods for fusion of remotely
% sensed data", Murphy & Maggioni.
%
% DM_dist  - Diffusion distance between points (NxN)
% Euc_dist - Euclidean distance between points (NxN)
% siz      - Image size (m x n = N)
% k        - Number of neigbors to use for density calculations
% sigma2   - Kernel scale
% k_modes  - Number of interest points (to become classes)

len = siz(1)*siz(2);

knn_dis = sort(Euc_dist,2,'ascend'); knn_dis = knn_dis(:,2:(k+1));
p0      = sum(exp(-knn_dis.^2/(sigma2*mean(knn_dis(:).^2))),2);
p       = p0/sum(p0);

[~,m_ind_p]      = max(p);
rho              = zeros(1,len);
DM_dist_inf_diag = DM_dist + diag(Inf*ones(1,len));
for ii  = 1:len
    if ii ~= m_ind_p
        rho(ii) = min(DM_dist_inf_diag(ii,(p>=p(ii))));
    else
        rho(ii) = max(DM_dist(ii,:));
    end
end
rho = rho/max(rho);
 
[~,modes_ind] = maxk(rho(:).*p(:),k_modes);

end

