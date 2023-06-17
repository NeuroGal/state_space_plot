function [dim_red_patterns, pca_info] = do_dim_red(segs, method, time_inds, varargin)
% Wrapper for Matlab dimensionality reduction functions. Performs
% dimensionality reduction over dim2 (supposed to be electrodes\other
% recording channel). Plot the output using plot_mv_patterns.m.
%
% Written while working on this paper (please cite):
%   Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
%
% Input details:
%   segs   - dims: n_stim(\conditions\trials) x n_elec(\other recording channels) x n_time (can also be 1)
%   method - for dimensionality reduction, options: 'tsne', 'pca', 'cmds'
%         (classical mds), 'mds' (non-classical mds)
%       * 'tsne', 'mds', 'cmds' are performed using euclidean distance
%         unless provided with a 'metric' optional argument.
%   time_inds - time-indices to run over (you may not want to plot
%         everything afterwards)
% 
% Optional: (all in the format of 'string', setting)
%   'metric' - default: euclidean, used only for 'tsne', 'mds', 'cmds'. 
%       Goes into Matlab's 'pdist' so you can use anything that recieves, 
%       e.g.: 'euclidean', 'correlation', 'cosine', 'mahalanobis'
%   'n_dim' - determines number of dimensions in the output (just 2\3, default 2).
%   'link_plots' - true\false, default: false.
%       method = 'pca': perform pca on all tp together (default is separately for each point).
%           * currently uses all segment time-points (not just time_inds) starting from 'pca_stat_idx'.
%       method = 'mds' or 'tsne' the previous time-point solution is fed as an initial
%           estimate to the next point. This eases comparison between time-points, 
%           but isn't recommended in case of large temporal gaps or fast changes
%           since the algorithm may not converge.
%       * TODO: doing tsne\mds on all tp together (just needs to add calculation of the distances between different timepoints)
%       * Recommended to match this to the settings in plot_mv_patterns.
%   'n_pc_tsne' - number of pcs to use for 'tsne'. Default is min(50, n_features)
%   'pca_start_idx' - what index (out of n_time) to start the pca on (goes together with 'link_plots')
%   'verbose' - print user feedback? true\false, default: true.
%
% Output:
%   dim_red_patterns: n_stim x n_dim x n_time_inds
%   pca_info: currently only has the field 'explained' with the % explained
%       variance by each pc (pc# x time). More fields will be added
%       according to the needs in plot_mv_patterns.
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2019-2022
% Send bug reports and requests to gal.vishne@gmail.com

[n_stim, n_elec, n_time] = size(segs); 
n_time_inds = length(time_inds); distance_metric = 'euclidean'; link_plots = false;
n_dim = 2; pca_start_idx = 1; n_pc_tsne = min(n_elec,50); verbose = true;
pca_info = [];

arg  = 1; 
while arg <= size(varargin,2)
    switch varargin{arg}
        case 'metric'
            distance_metric = varargin{arg+1};
        case 'n_dim'
            n_dim = varargin{arg+1};
            if ~(n_dim == 2 || n_dim == 3); warning('You tried n_dim = %d, using 2 pcs instead', n_dim); n_dim = 2; end
        case 'link_plots'
            link_plots = varargin{arg+1};
        case 'n_pc_tsne'
            n_pc_tsne = varargin{arg+1};
            if n_elec < n_pc_tsne; warning('Using %d pc for tsne (n_pc needs to be smaller than the number of features)', n_elec); n_pc_tsne = n_elec; end
        case 'pca_start_idx'
            pca_start_idx = varargin{arg+1};
            if pca_start_idx > n_time; warning('Out of range pca_start, using the start of the segments'); end
        case 'vebose'
            verbose = varargin{arg+1};
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
    arg = arg + 2;
end

dim_red_patterns = nan(n_stim, n_dim, n_time_inds);

if strcmp(method, 'tsne') && verbose; fprintf('Starting tsne calculations:\n'); end
% calculations
if strcmp(method, 'pca') && link_plots
    pca_dat = reshape(permute(segs(:,:,pca_start_idx:end),[2,1,3]), n_elec, n_stim*(n_time-pca_start_idx+1))';
    [pca_info.coeff,~,~,~,pca_info.explained,pca_mu] = pca(pca_dat, 'NumComponents', n_dim);
    pca_all_dat = reshape(permute(segs,[2,1,3]), n_elec, n_stim*n_time)';
    % not using the function's output because we want to project all data, not just the one used for PCA
    pca_scores = permute(reshape( ((pca_all_dat-pca_mu) *pca_info.coeff)',n_dim, n_stim, n_time),[2 1 3]); 
    dim_red_patterns = pca_scores(:,:,time_inds);
else
    pca_info.explained = nan(min(n_elec,n_stim-1), n_time_inds);
    for t_i = time_inds
        curr_patterns = segs(:,:,t_i); t_num = find(time_inds==t_i);
        if strcmp(method, 'mds') || strcmp(method, 'cmds')
            curr_patterns = squareform(pdist(curr_patterns, distance_metric));
        end
        switch method
            case 'tsne'
                if verbose; fprintf('done with %d/%d\n', t_num, n_time_inds); end
                inputs = {curr_patterns,'Distance',distance_metric, 'NumPCAComponents', n_pc_tsne, 'NumDimensions', n_dim};
                if t_num>1 && link_plots;inputs = [inputs, {'InitialY', dim_red_pat}];end
                dim_red_pat = tsne(inputs{:});
            case 'mds'
                inputs = {curr_patterns, n_dim};
                if t_num>1 && link_plots;inputs = [inputs, {'Start', dim_red_pat}];end
                dim_red_pat = mdscale(inputs{:});
            case 'cmds'
                dim_red_pat = cmdscale(curr_patterns, n_dim);
            case 'pca' % this is only the case where ~link_plots
                [~,dim_red_pat,~,~,pca_info.explained(:,t_num)] = pca(curr_patterns, 'NumComponents', n_dim);
            otherwise
                error(['Unknown dim red method: ' method '.']);
        end
        dim_red_patterns(:, :, t_num) = dim_red_pat;
    end
end
end