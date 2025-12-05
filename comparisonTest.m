% 10/23/2025 Makoto. Created for validating dipole density tools.
clear all
clc

addpath('/srv/Makoto/Tools/dipoledensitytools')
aal = ft_read_atlas('/srv/TOOLKITS/fieldtrip-master/template/atlas/aal/ROI_MNI_V4.nii');
EEG = pop_loadset('filename','sub-002.set','filepath','/srv/Makoto/Dortmunt/p0100_imported/', 'loadmode', 'info');

% 1. Newly developed one (optimized, fast, verified)
[peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm] = calc_peakDipoleDensity([20 20 20], aal.transform, 8);

figure
imagesc(x_mm, y_mm, ijk_hist(:,:,ijk_peak(3))'); axis xy

% 2. Conventional one.
plotDipoleDensitySingle([20 20 20], 8, aal.transform, 'axial')

% 3. Figuring out how to use 'mri3dplot'
load(EEG.dipfit.mrifile);

close all
h = vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm, 'Colormap', jet(256), 'Neurological', true)
h = vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm, 'Colormap', jet(256), 'Neurological', false)
h = vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm, 'Colormap', jet(256), 'View', 'sagittal')
h = vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm, 'Colormap', jet(256), 'View', 'coronal')
h = vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm, 'Colormap', jet(256), 'view', 'axial', 'Slice', 40, 'TitleStyle', 'peakXYZ')
h = vis_peakDipoleDensity(mri, peakXYZ_mm, peakCount, ijk_hist, ijk_peak, x_mm, y_mm, z_mm, 'Colormap', jet(256), 'view', 'axial', 'Slice', 40, 'TitleStyle', 'currentSlice')

