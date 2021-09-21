function [ VS, VA, DS, DA ] = riemann_wavelets( K, sav, r, nE, nL )
%RIEMANN_WAVELETS 
% Input:    K   -   cell array containing the SPD kernels
%           sav -   struct including a binary variable that indicates 
%                   whether to save the wavelet eigenvectors (1) or not (0)
%                   and a string with the file name: sav.yn, sav.fname.
%                   If sav.yn = 1 the function returns VS = [], VA = [].
%           r   -   kernel rank
%           nE  -   number of eigenvectors to keep
%           nL  -   number of levels to calculate

T = length(K);

if nargin<5
    nL = floor(log2(T));
end

if ~sav.yn
    VS_all = cell(1,nL);
    VA_all = cell(1,nL);
    
    DS_all = cell(1,nL);
    DA_all = cell(1,nL);
end

Sprev = K;

for lev = 1:nL
    
    Tc    = floor(length(Sprev)/2);
    Scurr = cell(1,Tc);
    VS_l  = cell(1,Tc);
    VA_l  = cell(1,Tc);
    DS_l  = cell(1,Tc);
    DA_l  = cell(1,Tc);
    
    ind   = 1;
    h     = waitbar(0,['Calculating kernels for level ',num2str(lev),' out of ',num2str(nL)]);
    
    for t = 1:2:(length(Sprev)-1)
        waitbar(t/length(Sprev),h);
        
        S          = FixedGeodes_eff( Sprev{t}, Sprev{t+1}, 0.5, r );
        Scurr{ind} = 0.5*(S + S.');
        A          = FixedGeodes_eff_proj( Scurr{ind}, Sprev{t+1}, r );
        A          = 0.5*(A + transpose(A));
        
        [VS_l{ind},DS] = eigs(Scurr{ind},nE); DS_l{ind} = diag(DS);
        [VA_l{ind},DA] = eigs(A,nE);          DA_l{ind} = diag(DA);
                
        ind = ind + 1;
    end
    Sprev = Scurr;
    close(h);
    
    if ~sav.yn
        VS_all{lev} = VS_l;
        VA_all{lev} = VA_l;
        
        DS_all{lev} = DS_l;
        DA_all{lev} = DA_l;
    else
        save([sav.fname,'_wavelet_eigs_level_',num2str(lev),'.mat'],'VS_l','VA_l','DS_l','DA_l');
    end

end

if ~sav.yn
    VS = VS_all{lev};
    VA = VA_all{lev};
    
    DS = DS_all{lev};
    DA = DA_all{lev};
else
    VS = [];
    VA = [];
    
    DS = [];
    DA = [];
end