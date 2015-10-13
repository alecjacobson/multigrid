prefix = 'Model10_';
OC = readDMAT([prefix '4-points.dmat']);
dt = [];
mt = [];
pt = [];
N = [];
% one-linear to use pcg with same prototype:

CTV = [];
CTT = [];
% Prepare tetmesh input
for r = 1:5
  suffix = sprintf('%d',r);
  if r == 5
    tri_filename = [prefix 'input' '.obj'];
  else
    tri_filename = [prefix suffix '.obj'];
  end
  [V,F] = load_mesh(tri_filename);
  tet_filename = [prefix sprintf('%d.mesh',r)];
  if ~exist(tet_filename,'file')
    max_vol = max(doublearea(V,F))^2;
    [TV,TT,TF] = tetgen(V,F,'Flags',sprintf('-q2a%0.17f',max_vol));
    writeMESH(tet_filename,TV,TT,TF);
  else
    [TV,TT,TF] = readMESH(tet_filename);
  end

  N = [N;size(TV,1)];

  % boudary conditions
  %b = snap_points(OC,TV);
  %bc = kron(ones(1,ceil(numel(b)/2)),[1 0])';
  %bc = bc(1:numel(b));
  b = unique(boundary_faces(TT));
  bc = zeros(numel(b),1);
  tic;
  L = cotmatrix3(TV,TT);
  M = massmatrix3(TV,TT,'barycentric');
  sA = -L;
  %sA = L*(M\L);
  sB = M*ones(size(L,1),1);
  fprintf('Setup: %gs\n',toc);

  % direct solve
  tic;
  Zd = min_quad_with_fixed(sA,sB,b,bc);
  dt = [dt;toc];
  fprintf('Direct solve: %gs\n',dt(end));

  iters = 1000000;
  tol = 1e-4;
  % Set up interpolation
  if ~isempty(CTV)
    tic;
      Z0 = prolong(CTV,CTT,CZ,TV);
      Z0(b) = bc;
      fprintf('Prolong: %gs\n',toc);
    tic;
      Zm = min_quad_with_fixed_pcg(sA,sB,b,bc,tol,iters,Z0);
      norm(Z0-Zm)/norm(Zm)
      medit(TV,TT,TF,'Data',Z0-Zm);
      mt = [mt;toc];
      fprintf('pcg (warm): %gs\n',mt(end));
  else
    Zm = Zd;
    mt = dt;
  end

  tic;
    Zp = min_quad_with_fixed_pcg(sA,sB,b,bc,tol,iters,[]);
    pt = [pt;toc];
  fprintf('pcg (cold): %gs\n',pt(end));

  a = 1;
  nsp = 3;
  subplot(1,nsp,1);
  tsurf(TF,TV,'FaceAlpha',a,'EdgeAlpha',a,'EdgeColor','none','CData',Zd,fphong);
  %hold on;scatter3(OC(:,1),OC(:,2),OC(:,3),'SizeData',100,'LineWidth',4);hold off
  %hold on;scatter3(TV(b,1),TV(b,2),TV(b,3),'SizeData',100,'LineWidth',4);hold off
  view(2);
  axis equal;
  subplot(1,nsp,2);
  tsurf(TF,TV,'FaceAlpha',a,'EdgeAlpha',a,'EdgeColor','none','CData',Zm,fphong);
  view(2);
  axis equal;
  drawnow;
  subplot(1,nsp,3);
  loglog(N,[dt mt pt],'LineWidth',4)
  axis square;

  CTV = TV;
  CTT = TT;
  CZ = Zm;
end
