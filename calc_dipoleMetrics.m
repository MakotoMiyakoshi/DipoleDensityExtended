function dipMetrics = calc_dipoleMetrics(dip_pos, dip_dir, mesh, varargin)
% DIPOLE_SURFACE_METRICS (with inside/outside + signed distance)
% Returns:
%   (1) distance along dipole axis to first surface hit (t_axis)
%   (2) angle between dipole axis and cortical normal (always computed)
%   (3) true nearest (axis-independent) distance to mesh
%   (4) inside/outside classification (ray-casting parity test)
%   (5) signed distance: negative if inside, positive if outside
%
% INPUTS
%   dip_pos : [N x 3] dipole origins (mm)
%   dip_dir : [N x 3] dipole directions (any length)
%   mesh.pos [V x 3], mesh.tri [F x 3]  % closed, watertight, manifold recommended
%
% OPTIONS (name-value)
%   'KNear'        : integer, candidate vertex count for pruning (default 64)
%   'InsideRayDir' : 1x3 ray direction for parity test (default [1 1e-4 2e-4])
%   'InsideEps'    : small tolerance for ray test (default 1e-12)
%
% OUTPUT struct fields:
%   t_axis, hit_point, hit_tri, u_hit, v_hit
%   angle_deg
%   nearest_dist, nearest_point, nearest_tri, u_bary, v_bary
%   is_inside              % logical Nx1
%   signed_nearest_dist    % nearest_dist with sign (inside -> negative)

p = inputParser;
p.addParameter('KNear', 64, @(x)isnumeric(x)&&isscalar(x)&&x>=3);
p.addParameter('InsideRayDir', [1 1e-4 2e-4], @(x)isnumeric(x)&&numel(x)==3);
p.addParameter('InsideEps', 1e-12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
K           = p.Results.KNear;
ray_dir_in  = p.Results.InsideRayDir(:).';
ray_dir_in  = ray_dir_in./max(norm(ray_dir_in), eps);
inside_eps  = p.Results.InsideEps;

V = double(mesh.pos);
F = double(mesh.tri);
nV = size(V,1); nF = size(F,1);

% Precompute adjacency & normals
VF   = vertexFaceAdjacency(F, nV);
TR   = triangulation(F, V);
vertN = vertexNormal(TR);

% Precompute global triangle arrays for inside/outside ray casting
A_all = V(F(:,1),:);
B_all = V(F(:,2),:);
C_all = V(F(:,3),:);

% Normalize dipole directions
d_norm = vecnorm(dip_dir,2,2); d_norm(d_norm==0) = eps;
D = dip_dir ./ d_norm;

numDips = size(dip_pos,1);
t_axis        = nan(numDips,1);
hit_point     = nan(numDips,3);
hit_tri       = zeros(numDips,1);
u_hit         = nan(numDips,1);
v_hit         = nan(numDips,1);
angle_deg     = nan(numDips,1);
nearest_dist  = nan(numDips,1);
nearest_point = nan(numDips,3);
nearest_tri   = zeros(numDips,1);
u_bary        = nan(numDips,1);
v_bary        = nan(numDips,1);
is_inside     = false(numDips,1);

for dipIdx = 1:numDips
    if mod(dipIdx, 1000) == 0
        disp(sprintf('%d/%d...', dipIdx, numDips))
    end

    p0 = dip_pos(dipIdx,:);
    d  = D(dipIdx,:);

    % ---- Candidate faces near dipole (for both ray + nearest) ----
    d2 = sum((V - p0).^2, 2);
    [~, nn_idx] = mink(d2, min(K, nV));
    candF = unique(cat(1, VF{nn_idx}));
    if isempty(candF), candF = (1:nF)'; end

    A = V(F(candF,1),:);
    B = V(F(candF,2),:);
    C = V(F(candF,3),:);

     % ---- (1) Ray intersection (Möller–Trumbore) ----
    [t_hit_i, u_i, v_i, idxLoc_hit] = rayIntersectMany(p0, d, A, B, C);

    if ~isnan(t_hit_i)
        % Hit data
        t_axis(dipIdx)       = t_hit_i;
        hit_tri(dipIdx)      = candF(idxLoc_hit);
        hit_point(dipIdx,:)  = p0 + t_hit_i * d;
        u_hit(dipIdx)        = u_i;
        v_hit(dipIdx)        = v_i;

        % >>> REMOVED: angle computation here <<<
        % (We no longer compute angle on hit_tri; we will always use nearest_tri.)
    end
    
    % ---- (2) True nearest distance & closest point ----
    [P_cand, dlist, u_nv, v_nv, idxLoc_min] = closestPointPointTriangleVec_bary(p0, A, B, C);
    [nearest_dist(dipIdx), jmin] = min(dlist);
    nearest_point(dipIdx,:) = P_cand(jmin,:);
    nearest_tri(dipIdx)     = candF(idxLoc_min(jmin));
    u_bary(dipIdx)          = u_nv(jmin);
    v_bary(dipIdx)          = v_nv(jmin);

    % ---- (NEW) Angle: ALWAYS from nearest point barycentrics on nearest_tri ----
    tri = F(nearest_tri(dipIdx),:);
    u = u_bary(dipIdx);
    v = v_bary(dipIdx);
    w = 1 - u - v;
    n  = w*vertN(tri(1),:) + u*vertN(tri(2),:) + v*vertN(tri(3),:);
    n  = n ./ max(norm(n), eps);
    angle_deg(dipIdx) = acosd( min(1,max(-1, dot(d, n))) );

    % ---- (3) Inside/Outside via robust ray-casting parity test ----
    % Use a fixed, slightly tilted direction to avoid edge/vertex degeneracy.
    is_inside(dipIdx) = pointInsideMesh_ParityTest(p0, ray_dir_in, A_all, B_all, C_all, inside_eps);
end

signed_nearest_dist = nearest_dist;
signed_nearest_dist(is_inside) = -signed_nearest_dist(is_inside);

dipMetrics = struct('t_axis', t_axis, ...
             'hit_point', hit_point, ...
             'hit_tri', hit_tri, ...
             'u_hit', u_hit, ...
             'v_hit', v_hit, ...
             'angle_deg', angle_deg, ...
             'nearest_dist', nearest_dist, ...
             'nearest_point', nearest_point, ...
             'nearest_tri', nearest_tri, ...
             'u_bary', u_bary, ...
             'v_bary', v_bary, ...
             'is_inside', is_inside, ...
             'signed_nearest_dist', signed_nearest_dist);
end

% ================== Helpers ==================

function VF = vertexFaceAdjacency(F, nV)
VF = cell(nV,1);
for fi = 1:size(F,1)
    v = F(fi,:);
    VF{v(1)} = [VF{v(1)}; fi];
    VF{v(2)} = [VF{v(2)}; fi];
    VF{v(3)} = [VF{v(3)}; fi];
end
end

function [tmin, u_best, v_best, idx_local] = rayIntersectMany(O, D, A, B, C)
% Ray→many triangles; nearest forward hit. Returns barycentric (u,v) at hit.
E1 = B - A; E2 = C - A;
h  = cross(repmat(D,size(E2,1),1), E2, 2);
a  = sum(E1 .* h, 2);
mask = abs(a) > 1e-12;
tmin = NaN; u_best = NaN; v_best = NaN; idx_local = NaN;
if ~any(mask), return; end
f  = 1 ./ a(mask);
S  = repmat(O, nnz(mask), 1) - A(mask,:);
u  = f .* sum(S .* h(mask,:), 2);
q  = cross(S, E1(mask,:), 2);
v  = f .* sum(repmat(D, nnz(mask), 1) .* q, 2);
t  = f .* sum(E2(mask,:) .* q, 2);
good = (u>=0) & (v>=0) & (u+v<=1) & (t>0);
if ~any(good), return; end
t_good = t(good);
[ tmin, j ] = min(t_good);
maskIdx = find(mask);
goodIdx = maskIdx(good);
idx_local = goodIdx(j);
u_best = u(good); u_best = u_best(j);
v_best = v(good); v_best = v_best(j);
end

function inside = pointInsideMesh_ParityTest(P, D, A, B, C, epsl)
% Robust parity (odd-even) test by casting a ray from P along D
% Counts triangle intersections with strict interior test to avoid double hits.
% Returns true if count is odd (inside), false if even (outside).
E1 = B - A; E2 = C - A;

% Möller–Trumbore components
h  = cross(repmat(D,size(E2,1),1), E2, 2);
a  = sum(E1 .* h, 2);

mask = abs(a) > epsl;                 % skip near-parallel triangles
if ~any(mask)
    inside = false; return;
end
f  = 1 ./ a(mask);
S  = repmat(P, nnz(mask), 1) - A(mask,:);
u  = f .* sum(S .* h(mask,:), 2);
q  = cross(S, E1(mask,:), 2);
v  = f .* sum(repmat(D, nnz(mask), 1) .* q, 2);
t  = f .* sum(E2(mask,:) .* q, 2);

% Strict interior (avoid counting along edges/vertices twice)
tol = 1e-9;
good = (u>tol) & (v>tol) & (u+v<1-tol) & (t>tol);

% Count intersections
cnt = nnz(good);
inside = mod(cnt,2)==1;
end

function [Pclosest, d, u_out, v_out, idx_local] = closestPointPointTriangleVec_bary(P, A, B, C)
% Closest point from SINGLE point P to MANY triangles (A,B,C) with barycentric (u,v).
% Returns arrays per triangle, plus local index mapping.
M  = size(A,1);
AB = B - A; AC = C - A; AP = bsxfun(@minus, P, A);

% Dot products
dABAB = sum(AB.*AB,2);
dABAC = sum(AB.*AC,2);
dACAC = sum(AC.*AC,2);
dABAP = sum(AB.*AP,2);
dACAP = sum(AC.*AP,2);

% Initial (unclamped) barycentric solve for projection to plane
D = dABAB.*dACAC - dABAC.^2; D(D==0) = eps;
u = (dACAC.*dABAP - dABAC.*dACAP) ./ D;
v = (dABAB.*dACAP - dABAC.*dABAP) ./ D;
w = 1 - u - v;

% Classify regions and clamp
isFACE = (u>=0) & (v>=0) & (w>=0);

% Vertex regions
isA = (dABAP <= 0) & (dACAP <= 0);
isB = ( (dABAP >= dABAB) & ( (dACAP - dABAC) <= 0 ) );
isC = ( (dACAP >= dACAC) & ( (dABAP - dABAC) <= 0 ) );

% Edge AB
vc = dABAP .* (dACAP - dABAC) - (dABAP - dABAB) .* dACAP;
isEAB = (vc <= 0) & (dABAP >= 0) & (dABAP <= dABAB) & ~(isFACE|isA|isB|isC);
% Edge AC
vb = (dABAP .* dACAC) - (dACAP .* dABAC);
isEAC = (vb <= 0) & (dACAP >= 0) & (dACAP <= dACAC) & ~(isFACE|isA|isB|isC);
% Edge BC
BP = bsxfun(@minus, P, B); CP = bsxfun(@minus, P, C);
dE = C - B; dF = B - C;
dbb = sum(dE.*dE,2);
tB  = sum(dE.*BP,2)./max(dbb,eps);
tC  = sum(dF.*CP,2)./max(dbb,eps);
isEBC = (tB>=0 & tB<=1) & ~(isFACE|isA|isB|isC|isEAB|isEAC);

% Assemble closest points per region
Pclosest = zeros(M,3);
u_out = zeros(M,1); v_out = zeros(M,1);

if any(isFACE)
  Pclosest(isFACE,:) = A(isFACE,:).*w(isFACE) + B(isFACE,:).*u(isFACE) + C(isFACE,:).*v(isFACE);
  u_out(isFACE) = u(isFACE); v_out(isFACE) = v(isFACE);
end
if any(isA)
  Pclosest(isA,:) = A(isA,:); u_out(isA)=0; v_out(isA)=0;
end
if any(isB)
  Pclosest(isB,:) = B(isB,:); u_out(isB)=1; v_out(isB)=0;
end
if any(isC)
  Pclosest(isC,:) = C(isC,:); u_out(isC)=0; v_out(isC)=1;
end
if any(isEAB)
  t = dABAP(isEAB)./max(dABAB(isEAB),eps);
  Pclosest(isEAB,:) = A(isEAB,:) + AB(isEAB,:).*t;
  u_out(isEAB) = t; v_out(isEAB) = 0;
end
if any(isEAC)
  t = dACAP(isEAC)./max(dACAC(isEAC),eps);
  Pclosest(isEAC,:) = A(isEAC,:) + AC(isEAC,:).*t;
  u_out(isEAC) = 0; v_out(isEAC) = t;
end
if any(isEBC)
  t = min(max(tB(isEBC),0),1);
  Pclosest(isEBC,:) = B(isEBC,:) + (C(isEBC,:)-B(isEBC,:)).*t;
  % On BC: weight for A is 0; adopt u=t, v=1-t.
  u_out(isEBC) = t; v_out(isEBC) = 1 - t;
end

d = sqrt(sum((Pclosest - P).^2, 2));

% local indices (1..M) for convenience
idx_local = (1:M).';
end
