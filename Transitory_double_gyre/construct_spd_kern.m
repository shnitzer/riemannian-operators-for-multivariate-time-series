function [ K ] = construct_spd_kern( data, epspar )
%CONSTRUCT_SPD_KERN constructs a symmetric and (semi-) positive-definite
%kernel, 'K' (NxN), for 'data' (N x feat), with the normalization factor 
%'epspar'.

nandata = 0;

if any(isnan(data))
    data_nan = data;
    data(isnan(data_nan)) = 1e5;
    nandata = 1;
end

dist = squareform(pdist(data));
if nandata
    dist_nan = squareform(pdist(data_nan));
    epsmed   = nanmedian(dist_nan(:).^2);
else
    epsmed   = median(dist(:).^2);
end
W = exp(-dist.^2/(epspar*epsmed));
P = diag(1./sum(W,1)) * W * diag(1./sum(W,1));
P = diag(1./sqrt(sum(P,1))) * P * diag(1./sqrt(sum(P,1)));
K = 0.5*(P + P.');

end

