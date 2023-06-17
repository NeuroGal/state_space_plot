function fig = plot_mv_patterns(dim_red_patterns, time_vec, varargin)
% Visualization of multivariate data (plots the output of do_dim_red.m)
%
% The format was originally intended for time-resolved data (e.g. M\EEG) so
%   it expects as input a time dimension, but you can also just plot one
%   time-point.
% Output: figure with 2D\3D representation of the given segments at specific
%   time-points. Default: each time-point in a different subplot, and
%   each segment is represented as a point, but also includes a mode for
%   plotting trajectories.
%
% Written while working on this paper (please cite):
%   Vishne et al., Cell Reports 2023 (biorxiv DOI, to be updated
%   when formally published): https://doi.org/10.1101/2022.08.02.502469
%   'Distinct Ventral Stream and Prefrontal Cortex Representational
%   Dynamics during Sustained Conscious Visual Perception'
% -> see Figure 2a and Supplementary Figure 11a examples of what it can do,
%    but there are a lot more options than this.
%
% If you are interested in producing these types of plots please do get in
% touch! Gal Vishne, gal.vishne@gmail.com
%
% * Uses cbrewer2 for default colors: https://github.com/scottclowe/cbrewer2
% * and tight_subplot.m (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
% -> in both cases the version I used is in the folder but don't forget to credit them too.
%
% Input details:
%   dim_red_patterns - dims: n_stim x n_dim x n_time (n_time is actually the number of time_indices inserted to do_dim_red.m)
%   time_vec         - in ms, corresponding to n_time. 
%
% Optional:
%   'trajectories' - followed by struct with optional fields:
%           * to simply move to trajectory plotting mode you can just insert [] here
%           > add_time_text - T\F (default: false) add more text next to time-points indicating the plotted time
%           > TODO: add back the 'insets' and 'add_main_ax' options
%   'categories' - followed by struct with fields: (grouping of the stimuli to categories\conditions)
%           > cat_inds (n_stim x 1) - use numbers from 1 to n_categ
%           > cat_names (cellstr to be used in legend)
%   'plot_settings' - followed by struct with fields: (all can be just true\false)
%           > mean - plot the mean of all stimuli\ all stimuli per group\category if 'categories' used too (see below). Default: false.
%           > std - show circles\spherese\cylinders representing std in each group\ across all points. default: false.
%           > single - plot single images\conditions (dim1 of the dim_red_patterns). default: true.
%           > use_sem - switches the std to sem (if relevant). default: true.
%   'color_code' - followed by a numeric array size n_stim x 1\2 or n_categ x 1\2:
%           Column 1 - color of the points, column 2 (if supplied) - color
%           for surrounding circles. (both are pooled together to the same colormap, so if you use 1:n in 
%           column 1 and 1:m in column 2 they will share colors, so start from n+1 in column 2 if you don't want that)
%           Default: each stimulus has it's own color.
%           If categories were inserted this takes precendece over the first color_code input, unless it is in 
%           length n_categ which enables coloring some categories in the same color (but still averaging separately).
%   'colormap' - followed by array n_colors x 3(rgb) x 1\n_time.
%           n_colors = n_stim\number of unique color_code indices\n_categories etc.
%           if the last dimension is 1 and you are using trajecoties mode
%           the entire trajectory will have the same color.
%   'viz_settings' - followed by a structure with the following fields:
%           > link_plots: T\F sync the axis limits of all plots (default: false, 
%           changed to true for trajectories mode) - match this to the settings in do_dim_red.m
%           > new_fig: T\F make new figure (default: true)
%           > show_ticks: T\F show ticks of axes (default: true)
%           > show_legends: T\F print legends (for category etc, default: true)
%           > fig_title: add figure title as annotation (default: [])
%           > panel_titles: default: cellstr(string(time_vec)+" ms");
%   'pca_info' - changes the mode to is_pca which currently adds an inset
%           with % explained variance by each PC, and adds PC# titles to the axes.
%           followed by a struct with the following fields: (created by
%           do_dim_red when you use 'pca' as your method)
%           > explained: % explained variance by each PC (#PC x time)
%           > add back the 'add_main_ax' option and plot angles etc, which
%           also require more exports from do_dim_red.m
%
% Cool bonuses:
%   'video' - turn output is a video with the trajectories progressing through time, instead of a static image.
%           followed by struct with fields:
%           > FrameRate
%           > filename (where to save the video)
%   'images' - add a depiction of the actual stimuli used in the study (*currently only grayscale)
%           followed by struct with fields:
%           > folder - where to find the images
%           > filenames - file names, cell array in the order of the first dimension of dim_red_patterns
%           Optional:
%           > frame_size for the image (colored according to color_code) - in pixels (insert a number that makes sense relative to the size of the image) [default: 45]
%           > side_crop  - remove the sides, also in pixels [default: 30]
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2019-2022
% Send bug reports and requests to gal.vishne@gmail.com

% parse inputs
[n_stim, n_dim, n_time] = size(dim_red_patterns);
viz_settings = []; plot_settings = [];
color_code  = (1:n_stim)'; cat_inds = ones(n_stim,1);
gen_inset = @(pos) axes('Units','Normalized','Position',pos,'FontSize',12.5);
arg  = 1; 
while arg <= size(varargin,2)
    switch varargin{arg}
        case 'categories'
            cat_inputs = varargin{arg+1};
        case 'color_code'
            color_code = varargin{arg+1};
        case 'video'
            video_settings = varargin{arg+1};
        case 'pca_info'
            pca_info = varargin{arg+1};
        case 'images'
            image_settings = varargin{arg+1};
        case 'colormap' 
            cmap = varargin{arg+1};
        case 'plot_settings'
            plot_settings = varargin{arg+1};
        case 'trajectories'
            traj_settings = varargin{arg+1};
        case 'viz_settings'
            viz_settings = varargin{arg+1};
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
    arg = arg + 2;
end

% continue parsing & setting defaults
vars_tf = {'categories','trajectories','is_pca','video','images'};
struct_names = {'cat_inputs','traj_settings','pca_info','video_settings','image_settings'};
for val = 1:length(vars_tf); eval(sprintf("%s = exist('%s','var');",vars_tf{val},struct_names{val})); end
if ~trajectories; traj_settings = []; end
[plot_settings,viz_settings,traj_settings] = input_parser(plot_settings,viz_settings,traj_settings,time_vec);
if trajectories; viz_settings.link_plots = true; end
mult_plot_mode = ~(video || trajectories || ~viz_settings.new_fig);
if any(strcmp(varargin,'color_code')) && ~plot_settings.single
    warning("Not using your colorcode - you didn't ask for single patterns")
end 
if categories
    cat_inds = cat_inputs.cat_inds; 
    if any(strcmp(varargin,'color_code'))
        if size(color_code,1)==n_stim % color code per stim
            if size(color_code,2)==2
            warning('You inserted a color-code per stimulus - using only the 1st color code (basic color-code is by categories)');
            color_code(:,2) = color_code(:,1);
            end
        elseif size(color_code,1)==numel(unique(cat_inds)) % color code per category (enable coloring multiple categories in the same way)
            color_code = color_code(cat_inds,:);
        end
    else
        color_code(:,1) = cat_inds;
    end
end
n_cat = length(unique(cat_inds));

% compute mean and SEM per category (plotted later if asked for it)
mean_all   = nan(n_cat, n_dim, n_time);
std_all    = nan(n_cat, n_dim, n_time); % will be sem if that's requested
for categ = 1:n_cat
    cat_dat = dim_red_patterns(cat_inds==categ,:,:);
    mean_all(categ,:,:) = mean(cat_dat,1); std_all(categ,:,:) = std(cat_dat,[],1);
    if plot_settings.use_sem; std_all(categ,:,:) = std_all(categ,:,:)/sqrt(sum(cat_inds==categ)); end
end

% check how many distinct colors we need (in the stim\category axis, not considering time yet):
n_needed_col = plot_settings.single*numel(unique(color_code)); % 0 if no single patterns plotted, # unique color code if yes [already incorporates the category case (see above)
if plot_settings.mean||plot_settings.std
    if ~categories
        n_needed_col = n_needed_col+1; % add one more color for the mean (i.e. just one if no single patterns)
    else
        n_needed_col = numel(unique(color_code(:,1))); % just the number we need for the categories
    end
end
[~, unique_cat_idx] = unique(cat_inds); unique_cat_idx = color_code(unique_cat_idx,1); % needed later (pointer to one stimulus of each category)
if ~categories; unique_cat_idx = n_needed_col; end % go to the last one in this case

% define colormap [n_stim\n_catetc x rgb x n_time]
if exist('cmap','var') % if user input given test it's fine [currently not adjusting it, but see suggestion below]
    % dim 1 (stim\category)
    if size(cmap,1)~=n_needed_col; error('Wrong number of colors in the first dimension'); end    % in the future: can add black\gray for the mean in the no categories case
    % dim 2 (rgb)
    if size(cmap,2)~=3; error('The colormap must have rgb as the second dimension'); end
    % dim 3 (time)
    if ndims(cmap)<3 % no time dimension, so we create it
        cmap = repmat(cmap,1,1,n_time);
        if trajectories; warning(['Note all time-points are going to be colored in the same way, ',...
                'to change that add a 3rd dimension for time (or let the defaults do their magic)']); end
            % in the future we can use adjust_cmap function in this file to account for this
    end
else % we need to define the colormap
    if ~trajectories % we don't really need a time axis
        if n_needed_col < 10; cbrewer_inputs = {'qual','Set1'}; else; cbrewer_inputs = {'div','Spectral'}; end
        cmap = cbrewer(cbrewer_inputs{:}, max(n_needed_col,3), 'spline'); cmap = cmap(1:n_needed_col,:);
        cmap = repmat(cmap,1,1,n_time);
    else
        if n_needed_col == 1
            cmap = permute(cbrewer('div', 'Spectral', n_time, 'spline'),[3 2 1]);
        elseif n_needed_col < 9
            cmap = nan(n_needed_col, 3, n_time); 
            cbrew_name = {'Purples','Reds','Greens','Blues','Oranges','PuRd','YlGnBu','Greys'}; % changed easily
            for col = 1:n_needed_col
                CT = cbrewer('seq', cbrew_name{col}, round(n_time*1.15), 'spline');
                cmap(col,:,:) = CT((round(n_time*1.15)-n_time+1):end,:)';
            end
        else
            error('Please insert a colormap - there are not enough distinct trajectory colors in the current defaults')
        end
    end
end
cmap = max(min(cmap,1),0);

% set up figure\ axes - the dimensions here are probably very specific to my screen 
if mult_plot_mode
    nrow = ceil(n_time/5); ncol = round(n_time/nrow);
    fig_size = [0.05 0.1 0.1+0.09*ncol 0.09+0.14*nrow];
    marg_w = [0.05 0.02]; gap = [0.08 0.02]; marg_h = [0.08 0.07];
    if viz_settings.show_ticks; gap(1) = gap(1)+0.025;end
    if categories; marg_w(1) = 0.12; fig_size(3) = fig_size(3)+0.06; end
    if images; fig_size(3:4) = fig_size(3:4)*1.7; end
    fig = figure('Units','Normalized','Position',fig_size);
    ha = tight_subplot(nrow, ncol, gap, marg_h, marg_w); % Nh, Nw, gap [gap_h gap_w], marg_h [lower upper], marg_w [left right]
    set(ha,'TickDir','out'); set(ha((n_time + 1) : end),'visible','off');
elseif viz_settings.new_fig
    fig = figure('Units','Normalized','Position',[0.1 0.1 0.5 0.5],'color',ones(1,3)*0.98); 
    if ~video; small_ax = axes('Units','Normalized','Position',[0.75 0.1 0.2 0.2*(0.8/0.62)],'TickDir','out'); grid on; end
    ax_handle = axes('Units','Normalized','Position',[0.23 0.1 0.62 0.8],'TickDir','out');
    if ~video
    if n_dim == 2
        set(small_ax,'XAxisLocation','origin','YAxisLocation','origin');
        linkaxes([ax_handle small_ax])
    else
        small_ax.XRuler.FirstCrossoverValue = 0; small_ax.XRuler.SecondCrossoverValue = 0; % X crossover with Y,Z axis
        small_ax.YRuler.FirstCrossoverValue = 0; small_ax.YRuler.SecondCrossoverValue = 0; % Y crossover with X,Z axis
        small_ax.ZRuler.FirstCrossoverValue = 0; small_ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with X,Y axis
        Link = linkprop([ax_handle small_ax],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);
    end
    end
else
    fig = gcf;
    ax_handle = gca;
end

% if you asked for a video - open the file here
if video
    set(fig,'units','centimeters','Visible', 'off');
    vid = VideoWriter(video_settings.filename);
    vid.FrameRate = video_settings.FrameRate;
    open(vid);
end
annotation('textbox',[0.005,0.9,0.1,0.1],'String',viz_settings.fig_title,'FontSize',18,'LineStyle','none');

% if you asked to add images of the stimuli - some setup here
if images
    % setup
    im_fields = {'frame_size','side_crop'}; im_defaults = {45,30};
    for f = 1:length(im_fields); if ~isfield(image_settings,im_fields{f}); image_settings.(im_fields{f}) = im_defaults{f}; end; end
    im_cmap = [nan(1,3);repmat(linspace(0,1,256)',1,3)];
    % preload images (assumes the files are all the same size)
    im_tmp = imread([image_settings.folder,image_settings.filenames{1}]); im_size = size(im_tmp); % get the size of one image to prepare   
    im_resid_size = im_size-2*image_settings.side_crop; im_new_size = im_resid_size+2*image_settings.frame_size;
    im_matrix = zeros([im_new_size, n_stim]);
    for stim = 1:n_stim
        im_tmp = imread([image_settings.folder,image_settings.filenames{stim}]) + 1; % add one because the first value of cmap will be for the frame
        im_matrix(image_settings.frame_size+(1:im_resid_size(1)),image_settings.frame_size+(1:im_resid_size(2)),stim) = ...
            im_tmp(image_settings.side_crop+(1:im_resid_size(1)),image_settings.side_crop+(1:im_resid_size(2)));
    end
    set(fig,'units','centimeters'); fig_pos = get(fig,'position');
    im_size = 0.08.*[1, fig_pos(3)/fig_pos(4)]; % size of images (normalized units)
end
    
% last settings before plotting
primary_plot_args = {30,[],'filled','Marker','o','MarkerFaceAlpha',0.8}; if trajectories; primary_plot_args{1} = 25; end
secondary_plot_args = {42,[],'Linewidth',1.25};
mean_plot_args = {80,[],'filled','Marker','p','MarkerEdgeColor','k'};
% settings for plotting category legends\ keeping consistent axis limits if this is a video \ trajectory with images
if plot_settings.single
    if categories; idx = unique_cat_idx; else; idx = 1; end
    cat_locs = dim_red_patterns(idx,:,:); cat_plot_args = primary_plot_args; lims_array = dim_red_patterns;
else
    cat_locs = mean_all; cat_plot_args = mean_plot_args; lims_array = mean_all;
end
if plot_settings.std; lims_array = [lims_array; mean_all+std_all; mean_all-std_all]; end
lims = cellfun_wrap(@(x) limer(x), num2cell(lims_array,[1 3])', true);
% add category legend
if categories && ~video && viz_settings.show_legends
    rel_tp = 1*(~trajectories) + round(n_time/2)*(trajectories); 
    if mult_plot_mode; ax_handle = ha(rel_tp); end
    leg = add_cat_legend(ax_handle, cat_inputs.cat_names, cmap(unique_cat_idx,:,rel_tp), cat_locs(:,:,rel_tp), cat_plot_args);
end

% finally getting to plot the main data!
for t_idx=1:n_time 
    if mult_plot_mode; ax_handle = ha(t_idx); axes(ax_handle); end
    cur_patterns = dim_red_patterns(:, :, t_idx);
    cur_colormap = cmap(:,:,t_idx);
    if plot_settings.single
        primary_plot_args{2} = cur_colormap(color_code(:,1),:);
        scatter_add(cur_patterns, primary_plot_args);
        if size(color_code,2)==2 % add secondary plot
            secondary_plot_args{2} = cur_colormap(color_code(:,2),:);
            scatter_add(cur_patterns, secondary_plot_args);
        end
    end
    if plot_settings.mean
        mean_plot_args{2} = cur_colormap(unique_cat_idx,:);
        scatter_add(mean_all(:,:,t_idx), mean_plot_args);
    end
    if plot_settings.std
        std_to_use = std_all(:,:,t_idx); % plot the average dispersion in both directions, probably change this if you opt for something else in the 2pt variance below
        if trajectories;std_to_use = mean(std_to_use,2);end
        for categ = 1:n_cat
            add_var_shade_1tp(mean_all(categ,:,t_idx), std_to_use(categ), cur_colormap(unique_cat_idx(categ),:));
        end
    end
    if trajectories && t_idx > 1 % add connecting lines between points
        sect_plot_args = {'LineWidth',0.5+4/n_time*t_idx,'Color',[]};
        if plot_settings.single
            for stim = 1:n_stim
                sect_plot_args{end} = cur_colormap(color_code(stim,1),:);
                curr_pts = squeeze(dim_red_patterns(stim,:,(t_idx-1):t_idx))';
                plot_add(curr_pts, sect_plot_args);
            end
        end
        if plot_settings.mean
            for categ = 1:n_cat
                sect_plot_args{end} = cur_colormap(unique_cat_idx(categ),:);
                curr_pts = squeeze(mean_all(categ,:,(t_idx-1):t_idx))';
                plot_add(curr_pts, sect_plot_args);
            end
        end
        if plot_settings.std
            for categ = 1:n_cat
                curr_pts = squeeze(mean_all(categ,:,(t_idx-1):t_idx))';
                curr_std = squeeze(std_all(categ,:,(t_idx-1):t_idx))';
                var_total = [mean(curr_std(1,:)), mean(curr_std(2,:))];
              %  var_total = nan(2,1);   % if we want to plot only the variability in the orthogonal direction
              %  var_total(1) = calc_var_dir(all_patterns(cat_inds==categ, :, t_idx-1), curr_pts(2,:)-curr_pts(1,:));
              %  var_total(2) = calc_var_dir(all_patterns(cat_inds==categ, :, t_idx), curr_pts(2,:)-curr_pts(1,:));
                add_var_shade_2tp(curr_pts(1,:), curr_pts(2,:), var_total, cur_colormap(unique_cat_idx(categ),:));
            end
        end
    end
    
    % some visual settings
    if ~trajectories || video
        title(viz_settings.panel_titles(t_idx))
    end
    set(gca,'FontSize',14);
    if ~viz_settings.show_ticks; set(gca,'XTick',[],'YTick',[],'ZTick',[]); end
    
    % add information about the % of explained variance by the PCs
    if is_pca && (~viz_settings.link_plots || t_idx == 1 || (video && ~trajectories)) % if link_plots add for t_idx=1, otherwise for every t_idx. if it's video but not trajectories replot this every frame.
        if trajectories
            pc_inset = gen_inset([0.825 0.75 0.1 0.175]); n_pc_to_show = 5;
        else
            set(gca,'Units','Normalized'); pos = get(gca,'Position'); pos = [pos(1:2)+pos(3:4)*0.8,pos(3:4)*0.16];
            pc_inset = gen_inset(pos); n_pc_to_show = n_dim;
        end
        bar(1:n_pc_to_show, pca_info.explained(1:n_pc_to_show,t_idx),0.7, 'EdgeColor','none', 'FaceColor',0.65*ones(1,3),'FaceAlpha',0.8); box off
        xlabel('PC'); xlim([0.5 n_pc_to_show+0.5])
        ylb = ylabel('%','Rotation',0,'Units','normalized'); ylb.Position = ylb.Position - [0.5/n_pc_to_show, 0.1, 0];
        set(pc_inset, 'FontSize', 12.5)
        if trajectories; title('Explained variance','FontWeight','normal'); end
        set(fig,'CurrentAxes',ax_handle)
    end
    
    if video || (images && trajectories); xlim(lims(1,:)); ylim(lims(2,:));if n_dim == 3; zlim(lims(3,:));end; end
    % add images
    if images && (video || ~trajectories || t_idx == round(n_time/2)) 
        if video % need to verify this is enough for video mode
            getframe(fig); % not sure why, but this fixes the fig size for the first frame
            delete(findobj(fig,'Tag','im_ax')); % delete previous images
        end
        for stim = 1:n_stim
            new_ax = open_new_ax(ax_handle, cur_patterns(stim,:), im_size);
            im_cmap(1,:) = cur_colormap(color_code(stim,1),:);
            imshow(im_matrix(:,:,stim),im_cmap,'Parent',new_ax);
            set(new_ax,'Tag','im_ax'); % relevant for the video case (to delete the previous time-point images)
        end
        set(fig,'CurrentAxes',ax_handle)
    end
    if video
        set(fig,'CurrentAxes',ax_handle)
        if categories && viz_settings.show_legends 
            leg = add_cat_legend(ax_handle, cat_inputs.cat_names, cur_colormap(unique_cat_idx,:), cat_locs(:,:,t_idx), cat_plot_args);
        end
        uistack(ax_handle,'bottom')
        if ~trajectories; hold off; end % so the next plot overrides the current things
        frame = getframe(fig); writeVideo(vid,frame);
        fprintf('Done with frame %d/%d\n',t_idx,n_time)
    end
end

% add legend about the distance between tp, and add some more markings of time
if trajectories && ~video
    time_spaces = unique(diff(time_vec)); % if this is just 1 value it means spacing is equal
    if numel(time_spaces)==1 && viz_settings.show_legends 
        points_inset = gen_inset([0.15 0.8 0.08 0.08]); cat_plot_args{2} = 0.65*ones(1,3);
        scatter(points_inset,[0 1], [0 0], cat_plot_args{:}); hold on
        plot(points_inset,[0 1], [0 0], 'LineWidth', 2, 'Color', [0.65*ones(1,3) 0.8])
        text(0.5,0.1,sprintf('%dms',time_spaces),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12.5);
        set(points_inset,'Visible','off','FontSize',12.5)
        set(fig,'CurrentAxes',ax_handle)
    end
    if ~traj_settings.add_time_text % minimal text
        txt_i = [find(time_vec==0) n_time];
    else % a lot of text
        if n_time > 10
            txt_i = linspecer(1,n_time,9); txt_i = txt_i(2:end);
        else
            txt_i = 1:n_time;
        end
    end
    txt = arrayfun(@(x) sprintf(' %dms',x), time_vec(txt_i),'UniformOutput',false);
    pts = reshape(permute(mean_all(1, :, txt_i),[2 1 3]),n_dim,[])';
    if n_dim == 2
        text(pts(:,1),pts(:,2),txt,'FontSize',12)
    else
        text(pts(:,1),pts(:,2),pts(:,3),txt,'FontSize',12)
    end
end
% final touches
if video; close(vid); end
if viz_settings.link_plots && ~(trajectories || video)
    if n_dim == 2
        linkaxes(ha,'xy');
    else
        all_lims = arrayfun(@(x) [x.XLim;x.YLim;x.ZLim] , ha , 'UniformOutput' , false);
        all_lims = cat(3,all_lims{:});
        new_lims = [squeeze(min(all_lims(:,1,:),[], 3)), squeeze(max(all_lims(:,2,:),[], 3))];
        Link = linkprop(ha,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);
        xlim(new_lims(1,:)); ylim(new_lims(2,:)); zlim(new_lims(3,:));
    end
end
if exist('small_ax','var')
    set(ax_handle,'visible','off','Tag','main_ax')
    set(small_ax,'GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[],'FontSize',12.5)
    if is_pca
        xlb = xlabel(small_ax,'PC1'); xlb.Position = [ax_handle.XAxis.Limits(2) 0 0];
        ylb = ylabel(small_ax,'PC2'); ylb.Position = [0 ax_handle.YAxis.Limits(2) 0]; 
        if n_dim == 3; zlb = zlabel(small_ax,'PC3'); zlb.Position = [0 0 ax_handle.ZAxis.Limits(2)]; end
    end
end
if exist('leg','var'); set(leg,'Position', [0.02 leg.Position(2:4)]); end
end

%% functions

function p = scatter_add(pts, plot_args)
dim = size(pts, 2);
if dim == 2
    p = scatter(pts(:,1),pts(:,2),plot_args{:});hold on
elseif dim == 3
    p = scatter3(pts(:,1),pts(:,2),pts(:,3),plot_args{:});hold on
else
    error('points should be 2/3d')
end
end

function plot_add(pts, plot_args)
dim = size(pts, 2);
if size(pts,1) ~= 2
    error('expected only 2 points')
end     
if dim == 2
    plot(pts(:,1),pts(:,2),plot_args{:});hold on
elseif dim == 3
    plot3(pts(:,1),pts(:,2),pts(:,3),plot_args{:});hold on
else
    error('points should be 2/3d')
end
end

function leg = add_cat_legend(ax_handle, cat_names, colors_all, dot_locs, plot_args)
tmp_ax = gca; set(gcf,'CurrentAxes',ax_handle); p = []; n_cat = length(cat_names);
for categ = 1:n_cat
    col_code = colors_all(categ,:);plot_args{2} = col_code;
    p(categ) = scatter_add(dot_locs(categ,:), plot_args);
end
leg_inputs = {p, cat_names,'AutoUpdate','off','Box','off','NumColumns',floor(n_cat/3),'Location','northwest', 'FontSize',12.5}; % 'Position',[0.02,0.8,0.05,0.1],
if plot_args{1}>30
    [leg,icons] = legend(leg_inputs{:});
    set(findobj(icons,'Type','patch'),'MarkerSize',6+(plot_args{1}-30)/12);
else
    leg = legend(leg_inputs{:}); % oddly numcolumns doesn't work when you ask for icons... so I am asking for them only if it will change something
end
set(gcf,'CurrentAxes',tmp_ax)
end

function [plot_settings,viz_settings,traj_settings] = input_parser(plot_settings,viz_settings,traj_settings,time_vec)
viz_fields = {'link_plots','new_fig','show_ticks','show_legends','fig_title','panel_titles'};
viz_defaults = {false,true,true,true,[],cellstr(string(time_vec)+" ms")};
for f = 1:length(viz_fields); if ~isfield(viz_settings,viz_fields{f}); viz_settings.(viz_fields{f}) = viz_defaults{f}; end; end
plot_fields = {'mean','std','single','use_sem'};
plot_defaults = {false,false,true,true};
for f = 1:length(plot_fields); if ~isfield(plot_settings,plot_fields{f}); plot_settings.(plot_fields{f}) = plot_defaults{f}; end; end
if ~(plot_settings.mean || plot_settings.single)
     plot_settings.mean = true;
end
if ~isfield(traj_settings,'add_time_text'); traj_settings.add_time_text = false; end
end

function cmap = adj_cmap(cmap, bright_const)
cmap_hsv = rgb2hsv(cmap);
cmap_hsv(:,3) = bright_const;
cmap = hsv2rgb(cmap_hsv);
end