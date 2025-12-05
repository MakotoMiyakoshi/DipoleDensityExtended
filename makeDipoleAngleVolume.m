function S = makeDipoleAngleVolume(dipPos_mm, angles_deg, T_mni, dim_ijk, varargin)
% Maps dipole angles (deg) into MNI ijk grid, smooths, returns mean angle & counts.
%
% REQUIRED
%   dipPos_mm : [N x 3] dipole XYZ in mm (same frame as T_mni)
%   angles_deg: [N x 1] angle per dipole (0–180 or 0–90 if folded)
%   T_mni     : [4 x 4] affine: ijk -> mm  (e.g., aal.transform)
%   dim_ijk   : [1 x 3] grid size (e.g., [91 109 91])
%
% OPTIONS (name/value)
%   'sigma_mm'  : [1x3], Gaussian sigma in mm (default [4 4 4])
%   'min_count' : scalar, min RAW count to report angle (default 1)
%   'mask_by'   : 'raw'|'smoothed' counts (default 'raw')
%
% OUTPUT
%   S.angle_mean_ijk : [dim] smoothed mean angle (NaN where masked)
%   S.count_ijk      : [dim] smoothed counts (float)
%   S.raw_count_ijk  : [dim] raw hit counts (integer)
%   S.inbounds_frac  : fraction of dipoles mapped inside the volume
%   S.T_mni, S.voxel_mm, S.dim

p = inputParser;
p.addParameter('sigma_mm',[4 4 4], @(x)isnumeric(x)&&numel(x)==3);
p.addParameter('min_count',1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('mask_by','raw', @(s)ischar(s)||isstring(s));
p.parse(varargin{:});
sigma_mm  = p.Results.sigma_mm;
min_count = p.Results.min_count;
mask_by   = lower(string(p.Results.mask_by));

% 1) Voxel size (norm of affine axes)
voxel_mm = [norm(T_mni(1:3,1)) norm(T_mni(1:3,2)) norm(T_mni(1:3,3))];

% 2) Map mm -> ijk (1-based). Use column homogeneous coords: ijk = inv(T)*mm
N = size(dipPos_mm,1);
X  = [dipPos_mm, ones(N,1)];        % [N x 4], mm as rows
ijk_h = (T_mni \ X.').';            % [N x 4], solves T*ijk = mm  (more stable than inv)
ijk   = ijk_h(:,1:3);               % 1-based floating ijk
ijk   = round(ijk);                 % nearest voxel

% 3) Keep only in-bounds
inb = ijk(:,1)>=1 & ijk(:,1)<=dim_ijk(1) & ...
      ijk(:,2)>=1 & ijk(:,2)<=dim_ijk(2) & ...
      ijk(:,3)>=1 & ijk(:,3)<=dim_ijk(3);
inbounds_frac = mean(inb);
ijk = ijk(inb,:); ang = angles_deg(inb);
if isempty(ijk)
    warning('All dipoles mapped out-of-bounds. Check coordinate frames / transform.');
end

% 4) Accumulate raw sum & raw counts
lin = sub2ind(dim_ijk, ijk(:,1), ijk(:,2), ijk(:,3));
volSum   = accumarray(lin, ang, [prod(dim_ijk),1], @sum, 0);
volCount = accumarray(lin, 1,   [prod(dim_ijk),1], @sum, 0);
volSum   = reshape(volSum,   dim_ijk);
volCount = reshape(volCount, dim_ijk);

% 5) Separable Gaussian smoothing (sum and counts)
sigma_vox = max(sigma_mm ./ voxel_mm, eps);
gx = gauss1d(sigma_vox(1)); gy = gauss1d(sigma_vox(2)); gz = gauss1d(sigma_vox(3));

sumS   = convn(convn(convn(volSum,   gx,'same'), permute(gy,[2 1 3]),'same'), permute(gz,[2 3 1]), 'same');
countS = convn(convn(convn(volCount, gx,'same'), permute(gy,[2 1 3]),'same'), permute(gz,[2 3 1]), 'same');

% 6) Mean angle
meanS = sumS ./ max(countS, eps);

% 7) Mask using RAW or SMOOTHED counts (default: RAW)
switch mask_by
    case "raw"
        mask = volCount < min_count;   % robust for sparse maps
    case "smoothed"
        mask = countS   < min_count;   % stricter; can over-mask
    otherwise
        error('mask_by must be ''raw'' or ''smoothed''.');
end
meanS(mask) = NaN;

% 8) Pack
S = struct();
S.angle_mean_ijk = meanS;
S.count_ijk      = countS;
S.raw_count_ijk  = volCount;
S.inbounds_frac  = inbounds_frac;
S.T_mni          = T_mni;
S.voxel_mm       = voxel_mm;
S.dim            = dim_ijk;

% --- quick console diagnostics ---
fprintf('[makeDipoleAngleVolume] in-bounds: %.1f%% | nonzero raw voxels: %d | min_count=%g (%s)\n', ...
    100*inbounds_frac, nnz(volCount), min_count, mask_by);

end

function g = gauss1d(sigma)
% centered, length ~ ceil(6*sigma)+1, normalized
L = max(3, ceil(6*sigma));
if mod(L,2)==0, L=L+1; end
x = ((1:L) - (L+1)/2);
g = exp(-(x.^2)/(2*sigma^2));
g = g / sum(g);
g = reshape(g, [numel(g) 1 1]);   % x-axis kernel
end
