function [peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm] = calc_peakDipoleDensity(XYZ_mm, affineMatrix, smoothSigma_mm)
% Build a 3-D spatial histogram on a voxel grid defined by a 4x4 affine
% (voxel->mm), find the peak, and return it in mm.
%
% INPUTS
%   XYZ_mm        : [N x 3] dipole positions in millimeters (same space as T_vox2mm)
%   T_vox2mm      : 4x4 affine mapping [i j k 1]^T (1-based voxel) -> [x y z 1]^T (mm)
%   smoothSigma_mm: Gaussian sigma in mm for optional smoothing (0/[] = no smoothing)
%
% OUTPUTS
%   peakXYZ_mm : [1 x 3] location (mm) of the histogram peak (bin center)
%   peakCount  : scalar peak count (after smoothing if applied)
%   ijk_hist   : [nx x ny x nz] 3-D histogram counts (double)
%   ijk_peak   : [1 x 3] voxel indices (1-based) of the peak
%
% Notes:
% - Uses inverse affine to map mm -> voxel index space, then bins by rounding.
% - Smoothing is applied in voxel units converted from mm via the affine.
% - Works with non-isotropic voxels; handles flips (e.g., -2 mm in X) correctly.
% 
% Example code:
% % Step 1. Prepare figure layout
% figure('Color','w');
% tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
% 
% % Pre-specify an axes where you want the overlay to go
% axTarget = nexttile(1);  % choose any tile or existing axes handle
% 
% % Step 2. Compute dipole density (example call)
% [peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm] = ...
%     calc_peakDipoleDensity([20 20 20], aal.transform, 8);
% 
% % Step 3. Generate the density overlay using vis_peakDipoleDensity as usual,
% %         but capture its handle (so we can grab the rendered objects)
% h = vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, ...
%     x_mm, y_mm, z_mm, 'View','axial','Slice',0,'Neurological',true,'ShowColorbar',false);
% 
% % Step 4. Copy the two overlayed images (MRI + density) to your target axes
% copyobj([h.axMRI.Children; h.axDensity.Children], axTarget);
% 
% % Step 5. Match visual settings
% axis(axTarget,'equal','xy'); grid(axTarget,'off');
% xlabel(axTarget,'X (mm)');
% ylabel(axTarget,'Y (mm)');
% title(axTarget,'Dipole Density (copied to specified axes)');
% 
% % Optional: delete the temporary figure that vis_peakDipoleDensity created
% delete(h.figure);

% Copyright (C) Makoto Miyakoshi, Cincinnati Children's Hospital Medical Center
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.
%
% History
% 10/23/2025 Makoto and ChatGPT5. Created.

if nargin < 3,              smoothSigma_mm = [];   end
if isempty(smoothSigma_mm), smoothSigma_mm = 0;    end

% Hard-code the grid size.
if isequal(affineMatrix(:,4), [92 -128 -74 1]')
    gridSize = [91 109 91];
else
    error('Non-MNI affine matrix detected. Not supported.')
end

% 1) Map mm -> voxel indices (1-based) with inverse affine
T_mm2vox = inv(affineMatrix);
n = size(XYZ_mm,1);
hom = [XYZ_mm, ones(n,1)];
ijk = (T_mm2vox * hom.').';
ijk = ijk(:,1:3);  % continuous voxel coordinates (1-based)

% 2) Round to nearest voxel and keep only in-bounds
ix = round(ijk(:,1));  iy = round(ijk(:,2));  iz = round(ijk(:,3));
inb = ix>=1 & ix<=gridSize(1) & iy>=1 & iy<=gridSize(2) & iz>=1 & iz<=gridSize(3);
ix = ix(inb);  iy = iy(inb);  iz = iz(inb);

% 3) Accumulate counts into a 3-D array
ijk_hist = zeros(gridSize, 'double');
if ~isempty(ix)
    % Use 3-D subscripts directly
    subs = [ix, iy, iz];
    ijk_hist = accumarray(subs, 1, gridSize, @sum, 0);
end

    %{
    % Visually confirm the 'tissue'.
    figure
    tiledlayout(7,13,"TileSpacing","none","Padding","tight")
    for zIdx = 1:gridSize(3)
        nexttile
        imagesc(ijk_hist(:,:,zIdx)'); axis xy
    end
    %}

% 4) Optional Gaussian smoothing (sigma given in mm → convert to voxels)
if smoothSigma_mm > 0
    % voxel sizes (mm) are the norms of the first 3 columns of T_vox2mm
    vx = norm(affineMatrix(1:3,1));
    vy = norm(affineMatrix(1:3,2));
    vz = norm(affineMatrix(1:3,3));
    sigma_vox = [smoothSigma_mm/vx, smoothSigma_mm/vy, smoothSigma_mm/vz];
    ijk_hist = imgaussfilt3(ijk_hist, sigma_vox);
end

% 5) Peak (argmax) in voxel space
[peakCount, linMax] = max(ijk_hist(:));
[pi, pj, pk] = ind2sub(size(ijk_hist), linMax);
ijk_peak = [pi pj pk];

% 6) Convert peak voxel index -> mm using forward affine
peakXYZ_mm = (affineMatrix * [pi pj pk 1].').';
peakXYZ_mm = peakXYZ_mm(1:3);

% 7) Generate xyz tick labels in mm.
nx = size(ijk_hist,1);
ny = size(ijk_hist,2);
nz = size(ijk_hist,3);
i = 1:nx;
j = 1:ny;
k = 1:nz;
x_mm = affineMatrix(1,1)*i + affineMatrix(1,4);  % -2*i + 92
y_mm = affineMatrix(2,2)*j + affineMatrix(2,4);  %  2*j - 128
z_mm = affineMatrix(3,3)*k + affineMatrix(3,4);  %  2*k - 74