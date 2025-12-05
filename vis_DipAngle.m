function vis_DipAngle(mesh, dipPos, dipMom, out, opts)
% VIS_DIPANGLE  Interactive QC of dipole–cortex angles (fast + transparent)
% mesh.pos [Vx3], mesh.tri [Fx3]
% dipPos [Nx3], dipMom [Nx3]
% out: struct with fields at least:
%   nearest_tri, u_bary, v_bary, nearest_point
%   (optional) hit_tri, u_hit, v_hit, hit_point, t_axis
%
% opts (all optional):
%   .useNearest     (true)     % nearest-tri + u_bary/v_bary; false=>use hit if available
%   .useFaceNorm    (false)    % false=vertex-interpolated; true=face normal
%   .signResolve    (true)     % flip dipole if opposed to normal (keep 0–180 but consistent)
%   .fold90         (false)    % mirror >90 -> 180-θ (polarity-invariant 0–90)
%   .indices        ([])       % which dipoles to browse (default: all)
%   .showMNIradial  (false)    % draw radial vector from [0 0 0]
%   .arcRadius      (5.0)      % mm, size of angle arc
%   .oneRingDepth   (1)        % 1 or 2 ring expansion for local patch
%   .alphaLocal     (0.55)     % opacity of local patch (will be capped to 0.35 if transparent)
%   .scaleDipole    (20)       % arrow length for dipole vector
%   .scaleNormal    (15)       % arrow length for normal vector
%   .bgMode         ('surface')% 'surface' | 'wire' | 'points' | 'off'
%   .bgTargetFaces  (20000)    % decimation target faces for background
%   .bgTargetVerts  (25000)    % point count for 'points' mode
%   .bgAlpha        (0.20)     % transparency for 'surface' background (0.1–0.3 good)
%   .camZoom        (1.0)      % global camera backoff multiplier (not auto-aligned)
%
% Makoto + ChatGPT, 2025

if nargin<5, opts = struct; end
opts = filldefaults(opts, struct( ...
  'useNearest',true,'useFaceNorm',false,'signResolve',true,'fold90',false, ...
  'indices',[],'showMNIradial',false,'arcRadius',5.0,'oneRingDepth',1, ...
  'alphaLocal',0.55,'scaleDipole',20,'scaleNormal',15, ...
  'bgMode','surface','bgTargetFaces',20000,'bgTargetVerts',25000, ...
  'bgAlpha',0.20,'camZoom',1.0 ...
));

V  = double(mesh.pos);
F  = double(mesh.tri);
TR = triangulation(F,V);

% Outward-enforced normals (once)
[faceN, vertN] = outwardNormals(TR);

% Which dipoles to browse
N = size(dipPos,1);
if isempty(opts.indices)
    idxUse = 1:N;
else
    idxUse = opts.indices(:)';   % row vector
end

% Precompute vertex->faces adjacency for local patch
FA = vertexAttachments(TR); % cell per vertex

% =================== Figure & static background ===================
hf = figure('Name','Dipole–Cortex Angle Inspector (fast+transparent)', ...
            'Color','w','Position',[80 60 1200 800]);
set(hf,'Renderer','opengl');         % needed for transparency (software OpenGL is fine)
set(hf,'GraphicsSmoothing','off');   % faster remotely

ax = axes('Parent',hf); hold(ax,'on'); axis(ax,'equal'); axis(ax,'vis3d'); grid(ax,'on');
xlabel(ax,'X (mm)'); ylabel(ax,'Y (mm)'); zlabel(ax,'Z (mm)');
set(hf,'ToolBar','figure');
view(ax,[135 20]);                   % default view; user can rotate freely
rotate3d(ax,'on');

% Headlight that follows view
hLight = camlight(ax,'headlight'); material(ax,'dull');
hRot = rotate3d(ax.Parent);
hRot.ActionPostCallback = @(~,~) refreshHeadlight();

% Lightweight background cortex (built once)
pGlobal = [];  % handle for toggling alpha
switch lower(opts.bgMode)
  case 'surface'
    [Fbg, Vbg] = tryReduce(F, V, opts.bgTargetFaces);
    pGlobal = patch('Faces',Fbg,'Vertices',Vbg,'Parent',ax, ...
      'FaceColor',[0.86 0.88 0.95], ...
      'EdgeColor','none', ...
      'FaceLighting','none', ...
      'FaceAlpha', opts.bgAlpha);   % transparency so deep dipoles are visible
  case 'wire'
    [Fbg, Vbg] = tryReduce(F, V, opts.bgTargetFaces);
    pGlobal = patch('Faces',Fbg,'Vertices',Vbg,'Parent',ax, ...
      'FaceColor','none','EdgeColor',[0.75 0.8 0.9],'EdgeAlpha',0.4,'LineWidth',0.5);
  case 'points'
    nv  = size(V,1); k = min(opts.bgTargetVerts, nv); idx = round(linspace(1,nv,k));
    pGlobal = scatter3(ax, V(idx,1),V(idx,2),V(idx,3), 2, [0.75 0.8 0.9], 'filled');
  case 'off'
    % nothing
  otherwise
    [Fbg, Vbg] = tryReduce(F, V, opts.bgTargetFaces);
    pGlobal = patch('Faces',Fbg,'Vertices',Vbg,'Parent',ax, ...
      'FaceColor',[0.86 0.88 0.95], ...
      'EdgeColor','none', ...
      'FaceLighting','none', ...
      'FaceAlpha', opts.bgAlpha);
end

% =================== UI ===================
uicontrol('Style','text','String','Index','Units','normalized','Position',[0.01 0.93 0.05 0.03],'BackgroundColor','w');
hIdx = uicontrol('Style','edit','String',num2str(idxUse(1)),'Units','normalized','Position',[0.06 0.93 0.06 0.04],...
    'Callback',@(h,~) redraw(str2double(h.String)));
uicontrol('Style','pushbutton','String','Prev','Units','normalized','Position',[0.13 0.93 0.06 0.04],...
    'Callback',@(~,~) step(-1));
uicontrol('Style','pushbutton','String','Next','Units','normalized','Position',[0.20 0.93 0.06 0.04],...
    'Callback',@(~,~) step(+1));
uicontrol('Style','checkbox','String','Use Nearest','Value',opts.useNearest,'Units','normalized','Position',[0.28 0.93 0.11 0.04],...
    'BackgroundColor','w','Callback',@(h,~) setflag('useNearest',logical(h.Value)));
uicontrol('Style','checkbox','String','Face Normals','Value',opts.useFaceNorm,'Units','normalized','Position',[0.40 0.93 0.11 0.04],...
    'BackgroundColor','w','Callback',@(h,~) setflag('useFaceNorm',logical(h.Value)));
uicontrol('Style','checkbox','String','Sign-resolve','Value',opts.signResolve,'Units','normalized','Position',[0.52 0.93 0.11 0.04],...
    'BackgroundColor','w','Callback',@(h,~) setflag('signResolve',logical(h.Value)));
uicontrol('Style','checkbox','String','Fold 0–90','Value',opts.fold90,'Units','normalized','Position',[0.64 0.93 0.10 0.04],...
    'BackgroundColor','w','Callback',@(h,~) setflag('fold90',logical(h.Value)));

% Transparency toggle for background cortex (only does something in 'surface' mode)
hTransp = uicontrol('Style','checkbox','String','Transparent cortex','Value',true, ...
    'Units','normalized','Position',[0.76 0.82 0.20 0.04],'BackgroundColor','w', ...
    'Callback',@(h,~) setCortexAlpha(logical(h.Value)));

hTxt = annotation(hf,'textbox',[0.76 0.92 0.22 0.07],'String','','EdgeColor','k','BackgroundColor',[1 1 1 0.85]);

% =================== Dynamic graphics placeholders (all tracked) ===================
hLocal   = gobjects(1,1);
hDip     = gobjects(1,1);
hDipDot  = gobjects(1,1);
hNorm    = gobjects(1,1);
hNormDot = gobjects(1,1);
hRad     = gobjects(1,1);
hArc     = gobjects(1,1);
hHitDot  = gobjects(1,1);

% Initial draw
redraw(idxUse(1));

% =================== nested helpers ===================
    function step(dir)
        k = str2double(hIdx.String);
        pos = find(idxUse==k,1); if isempty(pos), pos = 1; end
        pos = max(1,min(numel(idxUse), pos+dir));
        hIdx.String = num2str(idxUse(pos));
        redraw(idxUse(pos));
    end

    function setflag(name, val)
        opts.(name) = val;
        k = str2double(hIdx.String);
        redraw(k);
    end

    function setCortexAlpha(tf)
        if strcmpi(opts.bgMode,'surface') && isgraphics(pGlobal) && strcmpi(get(pGlobal,'Type'),'patch')
            if tf
                set(pGlobal,'FaceAlpha', opts.bgAlpha);
            else
                set(pGlobal,'FaceAlpha', 1.0);
            end
            drawnow limitrate;
        end
    end

    function refreshHeadlight()
        if isgraphics(hLight), delete(hLight); end
        hLight = camlight(ax,'headlight');
    end

    function redraw(k)
        if ~ismember(k, idxUse), return; end

        % delete ONLY dynamic overlays; keep background
        del(hLocal, hDip, hDipDot, hNorm, hNormDot, hRad, hArc, hHitDot);

        P0 = dipPos(k,:); 
        D  = dipMom(k,:); D = D ./ max(norm(D),eps);

        useNearest = opts.useNearest;
        haveHit = isfield(out,'t_axis') && ~isempty(out.t_axis) && k<=numel(out.t_axis) && ~isnan(out.t_axis(k));
        if ~useNearest && haveHit
            triIdx = out.hit_tri(k);
            if isfield(out,'u_hit') && isfield(out,'v_hit')
                u = out.u_hit(k); v = out.v_hit(k);
                NP = out.hit_point(k,:);
            else
                u = 1/3; v = 1/3; NP = incenter(TR, triIdx);
            end
        else
            triIdx = out.nearest_tri(k);
            u = out.u_bary(k); v = out.v_bary(k); NP = out.nearest_point(k,:);
        end
        w = 1 - u - v;
        triVerts = F(triIdx,:);

        % local patch (slightly more transparent to see arrows)
        facesLocal = localPatchFaces(triVerts, FA, opts.oneRingDepth, F);
        if ~isempty(facesLocal)
            hLocal = patch('Faces', F(facesLocal,:), 'Vertices', V, ...
                           'Parent', ax, 'FaceAlpha', min(opts.alphaLocal,0.35), ...
                           'EdgeColor', [0.6 0.6 0.7], 'FaceColor', [0.85 0.88 1.00], ...
                           'FaceLighting','none');
        end

        % normal (vertex or face)
        if opts.useFaceNorm
            n = faceN(triIdx,:);
        else
            n = w*vertN(triVerts(1),:) + u*vertN(triVerts(2),:) + v*vertN(triVerts(3),:);
            n = n ./ max(norm(n),eps);
        end

        % sign resolution / folding
        dVis = D;
        if opts.signResolve && dot(dVis,n) < 0, dVis = -dVis; end
        theta = acosd( min(1,max(-1, dot(dVis,n))) );
        if opts.fold90 && theta > 90, theta = 180 - theta; end

        % dipole arrow + origin dot (red)
        hDip = quiver3(ax, P0(1),P0(2),P0(3), dVis(1),dVis(2),dVis(3), opts.scaleDipole, ...
                       'LineWidth',2.5,'Color',[0.8 0.2 0.2],'MaxHeadSize',0.6);
        hDipDot = plot3(ax, P0(1),P0(2),P0(3), 'o', 'MarkerSize',5, ...
                        'MarkerFaceColor',[0.8 0.2 0.2], 'MarkerEdgeColor','none');

        % normal arrow + NP dot (blue)
        hNorm = quiver3(ax, NP(1),NP(2),NP(3), n(1),n(2),n(3), opts.scaleNormal, ...
                        'LineWidth',2.5,'Color',[0.2 0.2 0.8],'MaxHeadSize',0.6);
        hNormDot = plot3(ax, NP(1),NP(2),NP(3), 'o', 'MarkerSize',5, ...
                         'MarkerFaceColor',[0.2 0.2 0.8], 'MarkerEdgeColor','none');

        % optional MNI radial (green)
        if opts.showMNIradial
            r = P0 ./ max(norm(P0),eps);
            hRad = quiver3(ax, P0(1),P0(2),P0(3), r(1),r(2),r(3), opts.scaleNormal, ...
                           'LineWidth',2,'Color',[0.1 0.6 0.1],'MaxHeadSize',0.6);
        end

        % angle arc between dVis and n at NP
        hArc = drawAngleArc3D(ax, NP, n, dVis, opts.arcRadius, [0 0 0], 64);

        % hit point marker (if exists)
        if haveHit
            hp = out.hit_point(k,:);
            hHitDot = plot3(ax, hp(1),hp(2),hp(3), 's', 'MarkerSize',6, ...
                            'MarkerFaceColor',[1 0.7 0.3], 'MarkerEdgeColor','k');
        end

        % annotation
        defStr = sprintf('%s | %s | %s', ...
            tern(opts.useNearest,'nearest','hit'), ...
            tern(opts.useFaceNorm,'face-N','vertex-N'), ...
            tern(opts.signResolve,'sign-resolved', tern(opts.fold90,'folded 0–90','raw 0–180')));
        hTxt.String = sprintf('dipole #%d\n\\theta = %.2f°\n%s', k, theta, defStr);

        % legend (clean + informative)
        delete(findobj(ax,'Type','Legend'));
        if isgraphics(hDip) && isgraphics(hNorm)
            legend(ax, [hDip, hNorm], {'Dipole moment (red)', 'Surface normal (blue)'}, ...
                'TextColor','k', 'Location','northeast', 'Box','off', 'fontsize', 12);
        end

        drawnow limitrate;
    end
end

% =================== helpers ===================

function del(varargin)
for h = varargin
    if isgraphics(h{1}), delete(h{1}); end
end
end

function [faceN, vertN] = outwardNormals(TR)
V = TR.Points; F = TR.ConnectivityList;
faceN = faceNormal(TR);
vertN = vertexNormal(TR);
ctr   = mean(V,1);
if mean(sum((incenter(TR)-ctr).*faceN,2) > 0) < 0.5, faceN = -faceN; end
if mean(sum((V-ctr).*vertN,2) > 0) < 0.5, vertN = -vertN; end
faceN = faceN ./ max(vecnorm(faceN,2,2), eps);
vertN = vertN ./ max(vecnorm(vertN,2,2), eps);
end

function faces = localPatchFaces(triVerts, vertexAttachments, depth, F)
faces = unique([vertexAttachments{triVerts(1)}(:); ...
                vertexAttachments{triVerts(2)}(:); ...
                vertexAttachments{triVerts(3)}(:)]);
if depth<=1, return; end
ring1Verts = unique(F(faces,:));
for v = ring1Verts'
    faces = unique([faces; vertexAttachments{v}(:)]);
end
end

function h = drawAngleArc3D(ax, origin, v1, v2, radius, color, nseg)
if nargin<6 || isempty(color), color = [0 0 0]; end
if nargin<7, nseg = 64; end
v1 = v1 ./ max(norm(v1),eps);
v2 = v2 ./ max(norm(v2),eps);
u = v1;
w = v2 - dot(v2,u)*u;
if norm(w) < 1e-12
    p = origin + radius*u;
    h = plot3(ax, [origin(1) p(1)], [origin(2) p(2)], [origin(3) p(3)], '-', 'LineWidth',2, 'Color',color);
    return
end
w = w ./ norm(w);
theta = acos( min(1,max(-1, dot(v1,v2))) );
t = linspace(0,theta,nseg);
pts = origin + radius*(cos(t)'.*u + sin(t)'.*w);
h = plot3(ax, pts(:,1), pts(:,2), pts(:,3), 'k-', 'LineWidth',2, 'Color',color);
end

function [F2,V2] = tryReduce(F,V,targetFaces)
try
    [F2,V2] = reducepatch(F,V,targetFaces);
catch
    warning('reducepatch not available; using full-resolution background.');
    F2 = F; V2 = V;
end
end

function s = tern(cond,a,b), if cond, s=a; else, s=b; end, end

function o = filldefaults(o, d)
fn = fieldnames(d);
for i=1:numel(fn)
    if ~isfield(o,fn{i}) || isempty(o.(fn{i})), o.(fn{i}) = d.(fn{i}); end
end
end
