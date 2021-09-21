function n_ind = get_neighbors(ind,siz,r_s)
%GET_NEIGHBORS returns the indices of all the neigbors in a patch sized
%(2*'r_s'+1)^2 surrounding the pixel with index 'ind' in an image of size
%'siz'.


r_s_vec     = -r_s:1:r_s;
[row,col]   = ind2sub(siz,ind);

n_row                   = row + r_s_vec;
n_row(n_row < 1)        = [];
n_row(n_row > siz(1))   = [];

n_col                   = col + r_s_vec;
n_col(n_col < 1)        = [];
n_col(n_col > siz(2))   = [];

[n_col_all,n_row_all] = meshgrid(n_col,n_row);

n_ind = sub2ind(siz,n_row_all(:),n_col_all(:));
n_ind = setdiff(n_ind,ind,'stable');

end

