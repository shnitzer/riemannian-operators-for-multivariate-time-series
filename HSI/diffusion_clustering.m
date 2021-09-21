function [modes_ind,y_lbld] = diffusion_clustering(Euc_dis,DM_dis,DM_dis2, params)
% Implementation of "Diffusion geometric methods for fusion of remotely
% sensed data", Murphy & Maggioni.

[p,modes_ind]       = mode_detect(DM_dis,Euc_dis,params.siz,params.k,params.sigma2,params.k_modes);
DMdist_inf_diag     = DM_dis2 + diag(Inf*ones(1,params.siz(1)*params.siz(2)));
y_lbld              = spec_spat_lbl(DMdist_inf_diag,p,modes_ind,params.r_s,params.siz);

end


