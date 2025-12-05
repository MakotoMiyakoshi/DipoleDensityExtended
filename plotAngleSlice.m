function plotAngleSlice(S, varargin)
% plotAngleSlice(S, ...)  Visualize an axial slice of angle_mean_ijk
% REQUIRED
%   S : struct from makeDipoleAngleVolume
%
% OPTIONS
%   'z_mm'       : axial MNI z (mm). If given, overrides 'k'
%   'k'          : axial slice index (1..dim(3))
%   'clim'       : [min max] color limits in degrees (default [0 90] or [0 180])
%   'counts_min' : additional visibility threshold (default 0; S already masked)
%   'cmap'       : colormap name (default 'turbo')
%   'maskVol'    : logical [dim] mask to overlay (e.g., cortical-only), optional
%
% Notes: Uses imagesc with axis xy; orientation follows your AAL/mri
p = inputParser;
p.addParameter('z_mm', [], @(x)isnumeric(x)&&isscalar(x));
p.addParameter('k',    [], @(x)isnumeric(x)&&isscalar(x));
p.addParameter('clim', [], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('counts_min', 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('cmap','turbo', @(x)ischar(x)||isstring(x));
p.addParameter('maskVol',[], @(x)islogical(x)||isnumeric(x));
p.parse(varargin{:});
z_mm = p.Results.z_mm;
k    = p.Results.k;
clim = p.Results.clim;
cmin = p.Results.counts_min;
cmap = p.Results.cmap;
maskVol = p.Results.maskVol;

A   = S.angle_mean_ijk;
Cnt = S.count_ijk;
dim = S.dim;
T   = S.T_mni;
vox = S.voxel_mm;

% choose k
if ~isempty(z_mm)
    idx = round(([0 0 z_mm 1] / T.').');  % same trick as builder: row vec * inv(T)
    k = idx(3);
elseif isempty(k)
    k = round(dim(3)/2);
end
k = max(1, min(dim(3), k));

% slice and transparency by counts
slA = squeeze(A(:,:,k)).';     % transpose for axial display (rows=y, cols=x)
slC = squeeze(Cnt(:,:,k)).';

% apply extra count threshold if requested
slA(slC < cmin) = NaN;

% choose default clim
if isempty(clim)
    % guess whether your angles are folded or not
    mx = nanmax(A(:));
    clim = [0 (mx>100)*180 + (mx<=100)*90];  % 0–90 if small max, else 0–180
end

% plot
figure('Color','w','Position',[100 100 700 650])
imagesc(slA, clim);
axis image; axis xy; colormap(cmap); colorbar;
xlabel('i (left–right)'); ylabel('j (posterior–anterior)');
title(sprintf('Angle mean (deg) @ z-index k=%d (%.1f mm approx.)', k, estZmm(T,k)));

% overlay contours of a mask if given (e.g., cortical-only)
if ~isempty(maskVol)
    M = squeeze(maskVol(:,:,k)).';
    hold on; contour(M, [0.5 0.5], 'LineColor',[0 0 0], 'LineWidth',0.8);
end
end

function zmm = estZmm(T,k)
% estimate z (mm) of the center of voxel k along third axis
% world = T * [i j k 1]^T; for center, use i= (dimX+1)/2, j= (dimY+1)/2
zmm = ( [0 0 k 1] * T.' ); % not exact center, but good enough to label
zmm = zmm(3);
end
