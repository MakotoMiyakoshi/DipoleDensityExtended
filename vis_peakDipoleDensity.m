function h = vis_peakDipoleDensity(MRI, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm, varargin)
% vis_peakDipoleDensity() - Visualize dipole density over an MNI MRI template with correct mm axes.
%
% DESCRIPTION
%   Renders a slice of the dipole-density volume (typically from
%   calc_peakDipoleDensity) over a high-resolution MNI MRI (1 mm), with
%   proper world-axis labels in millimeters. The function supports drawing
%   directly into a user-specified axes (for tiled layouts/GUI) or into a
%   new figure. MRI brightness can be darkened independently from density,
%   and the density overlay uses a dead-zone + gamma ramp so weak near-zero
%   values do not haze the background.
%
% USAGE
%   h = vis_peakDipoleDensity(MRI, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, ...
%                             x_mm, y_mm, z_mm, 'key', value, ...);
%
% REQUIRED INPUTS
%   MRI          : EEGLAB MRI struct with fields .anatomy (nx×ny×nz) and
%                  .transform (4×4), or a numeric 3-D volume. If numeric,
%                  you must provide 'MRITransform' (4×4).
%   peakXYZ_mm   : [1×3] Peak location in mm (metadata for title).
%   peakCount    : Scalar peak value (metadata for title).
%   ijk_hist     : [nx×ny×nz] Dipole density volume (double).
%   ijk_peak     : [1×3] Peak voxel indices (1-based) in ijk_hist.
%   x_mm,y_mm,z_mm : 1-D coordinate vectors (mm) matching ijk_hist dims.
%
% NAME–VALUE PAIRS (optional)
%   Axes            : Target axes handle to draw into (default: create new figure).
%   View            : 'axial'|'coronal'|'sagittal' (default 'axial').
%   Slice           : Slice position of the selected View (mm by default; see SliceUnits). Default = through peak.
%   SliceUnits      : 'mm' (default) | 'index' (voxel index in ijk_hist).
%
%   % Density display
%   Colormap        : 3-col RGB colormap for density (default hot(256)).
%   CLim            : [cmin cmax] fixed density range (default auto via CLimQuantile).
%   CLimQuantile    : High quantile for cmax when CLim not given (default 0.995).
%
%   % Alpha shaping (overlay opacity)
%   Alpha           : Max opacity of density in [0..1] (default 1).
%   AlphaCutoff     : Dead-zone (normalized, 0..1). Values ≤ cutoff are fully
%                     transparent (default 0.001).
%   AlphaGamma      : Ramp nonlinearity (>0). <1 boosts mid-range; >1
%                     emphasizes peaks (default 0.4).
%   AlphaSmoothSigma: Gaussian smoothing sigma (pixels) for the alpha field
%                     (default 0; 0 = off).
%   MinPatchArea    : Remove small islands before alpha shaping (pixels; default 25; 0 = off).
%
%   % MRI underlay
%   MRIScale        : Brightness scale after normalization & gamma (default 0.4).
%   MRIGamma        : Gamma for MRI intensities (>1 darkens highlights; default 1.1).
%   MRIFixCLim      : If true, freeze MRI CLim to [0 1] so MRIScale truly dims (default true).
%
%   % Visuals
%   ShowColorbar    : true|false (default true; attached to density).
%   ShowContours    : true|false (default true; draws clean outlines over density).
%   ContourLevels   : Values OR normalized [0..1] levels. Default = [0.6 0.85]*cmax.
%   ContourColor    : [r g b] (default [0 0 0]).
%   ContourLineWidth: Scalar (default 1.2).
%   Neurological    : true|false. If true, flip X axes so right=right (default true).
%   MRITransform    : 4×4 affine (required only if MRI is numeric).
%   TitleStyle      : 'none'|'currentSlice'|'peakXYZ' (default 'peakXYZ').
%   omitAxesLabels  : true|false (default false). If true, hides x/y ticks, labels,
%                     and axis lines on the base axes to minimize visual clutter in
%                     dense tiled layouts; the title (and colorbar, if shown) remain.
%
% OUTPUT
%   h : struct with fields
%       .figure          – Figure handle (if created)
%       .axMRI           – Axes of MRI layer
%       .axDensity       – Axes of density overlay
%       .imgDensity      – Image handle of density
%       .view            – View string
%       .slice_idx_density – Slice index in ijk_hist
%       .slice_mm        – Slice position in mm
%       .CLim            – Applied density CLim [cmin cmax]
%
% EXAMPLES
%   % Basic (axial @ Z=0 mm)
%   vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, ...
%       x_mm, y_mm, z_mm, 'View','axial','Slice',0,'Neurological',true);
%
%   % Draw into existing tile and emphasize density
%   ax = nexttile;
%   vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, ...
%       x_mm, y_mm, z_mm, 'Axes',ax, ...
%       'MRIScale',0.35,'MRIGamma',1.3,'MRIFixCLim',true, ...
%       'Alpha',0.9,'AlphaCutoff',0.25,'AlphaGamma',0.6,'AlphaSmoothSigma',1, ...
%       'MinPatchArea',25,'ShowContours',true);
%
% HISTORY
%   10/23/2025 – Initial version (Makoto & ChatGPT-5).
%   10/24/2025 – Axes support, resize listeners, MRI dimming, robust CLim,
%                simplified alpha (cutoff+gamma), optional contours, doc update.
% -------------------------------------------------------------------------

% -------------------- Parse options --------------------
p = inputParser;
addParameter(p,'Axes',[],@(h) isempty(h) || isgraphics(h,'axes'));
addParameter(p,'View','axial',@(s)ischar(s)||isstring(s));
addParameter(p,'Slice',[],@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'SliceUnits','mm',@(s) any(validatestring(s,{'mm','index'})));

addParameter(p,'Colormap',jet(256),@(x)isnumeric(x)&&size(x,2)==3);
addParameter(p,'CLim',[],@(x)isnumeric(x)&&(isempty(x)||(isvector(x)&&numel(x)==2)));
addParameter(p,'CLimQuantile',0.995,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);

addParameter(p,'Alpha', 0.85,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'AlphaCutoff',0.001,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
addParameter(p,'AlphaGamma',0.4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'AlphaSmoothSigma',0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MinPatchArea', 0,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MRIScale',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MRIGamma',1.1,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MRIFixCLim',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'ShowColorbar',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'ShowContours',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'ContourLevels',[],@(x)isnumeric(x)||isempty(x));
addParameter(p,'ContourColor',[0 0 0],@(x)isnumeric(x)&&numel(x)==3);
addParameter(p,'ContourLineWidth',1.2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Neurological',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'MRITransform',[],@(x)isnumeric(x)&&isequal(size(x),[4 4]) || isempty(x));
addParameter(p,'TitleStyle','peakXYZ',@(s)ischar(s)||isstring(s));
addParameter(p,'omitAxesLabels',false,@(x)islogical(x)||ismember(x,[0 1]));

parse(p,varargin{:});
opt = p.Results;

% -------------------- Basic checks ---------------------
assert(ndims(ijk_hist)==3,'ijk_hist must be 3-D.');
[nx,ny,nz] = size(ijk_hist);
assert(numel(x_mm)==nx && numel(y_mm)==ny && numel(z_mm)==nz, ...
    'x_mm,y_mm,z_mm must match ijk_hist size.');

% -------------------- MRI input ------------------------
if isstruct(MRI)
    V_mri = MRI.anatomy;
    T_mri = MRI.transform;
else
    V_mri = MRI;
    if isempty(opt.MRITransform)
        warning('MRI numeric without MRITransform; inferring from density axes.');
        T_mri = eye(4);
        T_mri(1,1) = mean(diff(x_mm)); T_mri(1,4) = x_mm(1) - T_mri(1,1);
        T_mri(2,2) = mean(diff(y_mm)); T_mri(2,4) = y_mm(1) - T_mri(2,2);
        T_mri(3,3) = mean(diff(z_mm)); T_mri(3,4) = z_mm(1) - T_mri(3,3);
    else
        T_mri = opt.MRITransform;
    end
end
assert(ndims(V_mri)==3,'MRI must be 3-D.');

% MRI axis vectors (assume diag+offset affine as in EEGLAB templates)
xmri_all = T_mri(1,1).*(1:size(V_mri,1)) + T_mri(1,4);
ymri_all = T_mri(2,2).*(1:size(V_mri,2)) + T_mri(2,4);
zmri_all = T_mri(3,3).*(1:size(V_mri,3)) + T_mri(3,4);

% -------------------- Pick slice -----------------------
view = lower(string(opt.View));
switch view
    case "axial"    % Z-const → (X,Y)
        idxD = resolveSlice(opt.Slice, z_mm, ijk_peak(3), opt.SliceUnits);
        slice_mm = z_mm(idxD);
        Dslice   = ijk_hist(:,:,idxD);
        idxM     = nearestIndex(zmri_all, slice_mm);
        Mslice   = double(V_mri(:,:,idxM));
        x_den_mm = x_mm;  y_den_mm = y_mm;
        x_mri_mm = xmri_all; y_mri_mm = ymri_all;

    case "coronal"  % Y-const → (X,Z)
        idxD = resolveSlice(opt.Slice, y_mm, ijk_peak(2), opt.SliceUnits);
        slice_mm = y_mm(idxD);
        Dslice   = squeeze(ijk_hist(:,idxD,:));
        idxM     = nearestIndex(ymri_all, slice_mm);
        Mslice   = double(squeeze(V_mri(:,idxM,:)));
        x_den_mm = x_mm;  y_den_mm = z_mm;
        x_mri_mm = xmri_all; y_mri_mm = zmri_all;

    case "sagittal" % X-const → (Y,Z)
        idxD = resolveSlice(opt.Slice, x_mm, ijk_peak(1), opt.SliceUnits);
        slice_mm = x_mm(idxD);
        Dslice   = squeeze(ijk_hist(idxD,:,:));
        idxM     = nearestIndex(xmri_all, slice_mm);
        Mslice   = double(squeeze(V_mri(idxM,:,:)));
        x_den_mm = y_mm;  y_den_mm = z_mm;
        x_mri_mm = ymri_all; y_mri_mm = zmri_all;

    otherwise
        error('Unknown View: %s', opt.View);
end

% Orient for imagesc(..., axis xy)
Dslice = Dslice.';  % columns → X, rows → Y
Mslice = Mslice.';

% Neurological vs Radiological (flip X axis vectors only)
if ~opt.Neurological
    x_den_mm = fliplr(x_den_mm);
    x_mri_mm = fliplr(x_mri_mm);
end

% -------------------- Axes setup -----------------------
if ~isempty(opt.Axes)
    ax1 = opt.Axes;
    fig = ancestor(ax1,'figure');
    ax2 = axes('Parent',fig,'Position',get(ax1,'Position'), ...
               'Color','none','XAxisLocation','top','YAxisLocation','right', ...
               'HitTest','off','HandleVisibility','off');
else
    fig = figure('Color','w');
    ax1 = axes('Parent',fig); axis(ax1,'xy'); axis(ax1,'equal'); hold(ax1,'on');
    ax2 = axes('Parent',fig,'Position',get(ax1,'Position'), ...
               'Color','none','XAxisLocation','top','YAxisLocation','right', ...
               'HitTest','off','HandleVisibility','off');
end
set([ax1 ax2],'Units','normalized');

% Keep ax2 glued to ax1 on resize/redraw
lA = addlistener(fig, 'SizeChanged', @(~,~) set(ax2,'Position',get(ax1,'Position')));
lB = addlistener(ax1, 'MarkedClean', @(~,~) set(ax2,'Position',get(ax1,'Position')));
setappdata(ax1,'OverlaySyncListeners',[lA lB]);

% -------------------- MRI render -----------------------
axis(ax1,'xy'); axis(ax1,'equal'); hold(ax1,'on');

lo = prctile(Mslice(:),1);
hi = prctile(Mslice(:),99);
if hi > lo
    Mslice = (Mslice - lo) / (hi - lo);
    Mslice = min(max(Mslice,0),1);
else
    Mslice = zeros(size(Mslice));
end
Mslice = Mslice .^ opt.MRIGamma;
Mslice = Mslice * opt.MRIScale;

imagesc(ax1, x_mri_mm, y_mri_mm, Mslice);
colormap(ax1, gray(256));
if opt.MRIFixCLim
    set(ax1,'CLim',[0 1],'CLimMode','manual');
end

% -------------------- Density render -------------------
axis(ax2,'xy'); axis(ax2,'equal'); hold(ax2,'on');
D = Dslice;  % keep raw values
img2 = imagesc(ax2, x_den_mm, y_den_mm, D);
colormap(ax2, opt.Colormap);

% Density CLim (robust unless provided)
if isempty(opt.CLim)
    Dflat = D(:); Dflat = Dflat(isfinite(Dflat));
    if isempty(Dflat), cMax = 1;
    else
        cMax = quantile(Dflat, opt.CLimQuantile);
        if cMax <= 0
            cMax = max(Dflat);
            if cMax <= 0, cMax = 1; end
        end
    end
    cLimD = [0 cMax];
else
    cLimD = opt.CLim;
end
set(ax2,'CLim',cLimD);

% ---- Alpha shaping: dead-zone + gamma + optional smoothing ----
if cLimD(2) > cLimD(1)
    Dnorm = (D - cLimD(1)) / (cLimD(2) - cLimD(1));
else
    Dnorm = zeros(size(D));
end
Dnorm = min(max(Dnorm,0),1);

cut  = opt.AlphaCutoff;             % ≤cut → fully transparent
mask = Dnorm > cut;

if opt.MinPatchArea > 0
    try
        mask = bwareaopen(mask, round(opt.MinPatchArea));
    catch
        % no IPT: skip cleanup
    end
end

r = zeros(size(Dnorm));
if cut < 1
    r(mask) = ((Dnorm(mask) - cut) ./ (1 - cut)) .^ opt.AlphaGamma;
end

if opt.AlphaSmoothSigma > 0
    try
        r = imgaussfilt(r, opt.AlphaSmoothSigma);
    catch
        % no IPT: skip smoothing
    end
end

A = opt.Alpha * min(max(r,0),1);
set(img2,'AlphaData',A);

% ---- Optional contours over density ----
if opt.ShowContours
    if isempty(opt.ContourLevels)
        levels = cLimD(1) + [0.60 0.85] * (cLimD(2)-cLimD(1));
    else
        if all(opt.ContourLevels >= 0 & opt.ContourLevels <= 1)  % normalized
            levels = cLimD(1) + opt.ContourLevels * (cLimD(2)-cLimD(1));
        else
            levels = opt.ContourLevels;
        end
    end
    [~, hc] = contour(ax2, x_den_mm, y_den_mm, D, levels, ...
                      'LineColor', opt.ContourColor, ...
                      'LineWidth', opt.ContourLineWidth);
    set(hc,'HitTest','off');
end

% -------------------- Sync & labels --------------------
linkaxes([ax1,ax2],'xy');
set(ax1,'XLim',[min([x_mri_mm(:); x_den_mm(:)]) max([x_mri_mm(:); x_den_mm(:)])], ...
        'YLim',[min([y_mri_mm(:); y_den_mm(:)]) max([y_mri_mm(:); y_den_mm(:)])]);
set(ax2,'XLim',get(ax1,'XLim'),'YLim',get(ax1,'YLim'),'Visible','off');

if ~opt.omitAxesLabels
    xlabel(ax1, sprintf('%s (mm)', axisLabelX(view)));
    ylabel(ax1, sprintf('%s (mm)', axisLabelY(view)));
else
    % Hide ticks, tick labels, and axis lines — keep the title visible
    set(ax1, 'XTick', [], 'YTick', [], ...
             'XTickLabel', [], 'YTickLabel', [], ...
             'XColor', 'none', 'YColor', 'none');
    box(ax1, 'off');   % optional: remove the box outline too
end

switch lower(string(opt.TitleStyle))
    case "none"
    case "currentslice"
        switch view
            case "axial",    title(ax1, sprintf('Z = %.1f mm', slice_mm));
            case "coronal",  title(ax1, sprintf('Y = %.1f mm', slice_mm));
            case "sagittal", title(ax1, sprintf('X = %.1f mm', slice_mm));
        end
    otherwise
        title(ax1, sprintf('Peak at [%.1f  %.1f  %.1f] mm, count=%.3g', ...
            peakXYZ_mm(1), peakXYZ_mm(2), peakXYZ_mm(3), peakCount));
end

% Colorbar (preserve axes positions to avoid shrinking)
if opt.ShowColorbar
    axPos = get(ax1,'Position');
    cb = colorbar(ax2,'Location','eastoutside'); %#ok<NASGU>
    ylabel(cb,'Dipole density (a.u.)');
    set(ax1,'Position',axPos); set(ax2,'Position',axPos);
end
set(ax1,'Layer','top'); grid(ax1,'off');

% -------------------- Outputs -------------------------
h = struct('figure',fig,'axMRI',ax1,'axDensity',ax2,'imgDensity',img2, ...
           'view',char(view),'slice_idx_density',idxD,'slice_mm',slice_mm,'CLim',cLimD);

end % ===== main =====

% ==================== Helpers =========================
function idx = resolveSlice(slc, mmvec, default_idx, units)
    if isempty(slc), idx = default_idx; return; end
    switch lower(string(units))
        case "mm",    [~, idx] = min(abs(mmvec - slc));
        case "index", idx = max(1, min(numel(mmvec), round(slc)));
        otherwise,    error('SliceUnits must be ''mm'' or ''index''.');
    end
end

function idx = nearestIndex(vec, val)
    [~, idx] = min(abs(vec - val));
end

function lab = axisLabelX(view)
    switch lower(string(view))
        case "axial",    lab = 'X';
        case "coronal",  lab = 'X';
        case "sagittal", lab = 'Y';
    end
end

function lab = axisLabelY(view)
    switch lower(string(view))
        case "axial",    lab = 'Y';
        case "coronal",  lab = 'Z';
        case "sagittal", lab = 'Z';
    end
end