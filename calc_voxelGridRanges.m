function [rangeCenters, rangeEdges] = calc_voxelGridRanges(T, gridSize)
% T: 4x4 voxel->mm affine (1-based indices)
% gridSize: [nx ny nz]

nx = gridSize(1); ny = gridSize(2); nz = gridSize(3);

% --- Ranges of voxel CENTERS (i,j,k = 1..nx etc.)
C = [ 1   1   1;
      nx  1   1;
      1   ny  1;
      1   1   nz;
      nx  ny  1;
      nx  1   nz;
      1   ny  nz;
      nx  ny  nz ];
C = [C, ones(8,1)];
XYZc = (T * C.').'; XYZc = XYZc(:,1:3);
rangeCenters = [min(XYZc); max(XYZc)];  % rows: [xmin ymin zmin; xmax ymax zmax]

% --- Full spatial EXTENT (edges), i.e., centers ± 0.5 voxel
E = [ 0.5     0.5     0.5;
      nx+0.5  0.5     0.5;
      0.5     ny+0.5  0.5;
      0.5     0.5     nz+0.5;
      nx+0.5  ny+0.5  0.5;
      nx+0.5  0.5     nz+0.5;
      0.5     ny+0.5  nz+0.5;
      nx+0.5  ny+0.5  nz+0.5 ];
E = [E, ones(8,1)];
XYZe = (T * E.').'; XYZe = XYZe(:,1:3);
rangeEdges = [min(XYZe); max(XYZe)];
end