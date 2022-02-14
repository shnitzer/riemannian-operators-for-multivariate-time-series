function y_lbld = spec_spat_lbl(DM_dist_inf,p,modes_ind,r_s,siz)

% Implementation of "Diffusion geometric methods for fusion of remotely
% sensed data", Murphy & Maggioni.
%
% DM_dist_inf - Diffusion distance between points with Inf on its diagonal(NxN)
% p           - 
% modes_ind   - Indices of chosen modes
% r_s         -
% siz         - Image size (m x n = N)

y_lbld              = zeros(siz);
y_lbld(modes_ind)   = (1:length(modes_ind)).';

[~,srtd_ind_p]      = sort(p,'descend');
srtd_ind_p_unlbld   = setdiff(srtd_ind_p,modes_ind,'stable');

% 1st iteration:
for ii = srtd_ind_p_unlbld.'
    % Compute spatial consensus
    n_ind   = get_neighbors(ii,siz,r_s);
    lbl_frq = histcounts(y_lbld(n_ind),(0:length(modes_ind)) + 0.5);
    if any(lbl_frq > length(n_ind)/2)
        [~,spat_cons_lbl] = max(lbl_frq);
    else
        spat_cons_lbl = [];
    end
    
    % Compute the spectral label
    tmp_inds     = find((p>=p(ii)) & (y_lbld(:)~=0));
    [~,spec_ind] = min(DM_dist_inf(ii,tmp_inds));
    spec_lbl     = y_lbld(tmp_inds(spec_ind));
    
    if isempty(spat_cons_lbl) || (spec_lbl==spat_cons_lbl)
        y_lbld(ii) = spec_lbl;
    end
end

srtd_ind_p_unlbld   = setdiff(srtd_ind_p,[modes_ind;find(y_lbld(:)~=0)],'stable');
% 2nd iteration:
for ii = srtd_ind_p_unlbld.'
    % Compute spatial consensus
    n_ind   = get_neighbors(ii,siz,r_s);
    lbl_frq = histcounts(y_lbld(n_ind),(0:length(modes_ind)) + 0.5);
    if any(lbl_frq > length(n_ind)/2)
        [~,spat_cons_lbl] = max(lbl_frq);
        y_lbld(ii)        = spat_cons_lbl;
    else
        % Compute the spectral label
        tmp_inds     = find((p>=p(ii)) & (y_lbld(:)~=0));
        [~,spec_ind] = min(DM_dist_inf(ii,tmp_inds));
        spec_lbl     = y_lbld(tmp_inds(spec_ind));
        y_lbld(ii)   = spec_lbl;
    end
end

end

