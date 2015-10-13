dim = 3;
switch dim
case 2
  [V,T] = load_mesh('woody.obj');
  [V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  V = V(:,1:2);
  F = [];
  M = massmatrix(V,T,'voronoi');
  B = M*ones(size(V,1),1);
  system_fun = @(V,T,A,B) ...
    poisson_system( ...
      V,T,A,B, ...
      unique(outline(T)), ...
      zeros(numel(unique(outline(T))),1));
  max_size = 100;
case 3
  %[V,F] = load_mesh('knight.ply');
  %[V,F] = decimate_cgal(V,F,0.25);
  %[V,F] = meshfix(V,F);
  %h = mean(mean(edge_lengths(V,F)));
  %max_vol = h^3;
  %[V,T] = tetgen(V,F,'Flags',sprintf('-a%0.17f -q2',max_vol));
  %%writeMESH('knight.mesh',V,T,F);
  %writeMESH('knight-half.mesh',V,T,F);
  %[V,T,F] = readMESH('knight.mesh');
  [V,T,F] = readMESH('knight-half.mesh');
  M = massmatrix(V,T,'barycentric');
  system_fun = @(V,T,A,B) ...
    poisson_system( ...
      V,T,A,B, ...
      unique(boundary_faces(T)), ...
      zeros(numel(unique(boundary_faces(T))),1));
  B = M*ones(size(V,1),1);
  %[V,F] = load_mesh('knight.ply');
  %[V,F] = decimate_cgal(V,F,0.25);
  %medit(V,T,F);
  max_size = 1000;
end

% jacobi doesn't work for bilaplacian because it is not diagonally dominant
relax_method = 'sor';
relax_weight = 1.0;
% Bilaplacian gets confused with too many pre-smooth iterations...
pre_jac = 2^0;
% Bilaplacian requires a lot more iterations...
post_jac = 2^5;

tic;
[Aex,Bex] = system_fun(V,T,[],B);
Zex = Aex\Bex;
fprintf('backslash: %g\n',toc);

s = semilogy([nan nan;nan nan],[nan nan;nan nan],'-o','LineWidth',3);
Z = zeros(size(V,1),1);
E = [];
min_E = 1e-12;
data = [];

tic;
while true

  [Z(:,1),data] = multigrid(system_fun,B,V,T,'Z0',Z(:,1), ...
    'BoundaryFacets',F, ...
    'PreJacobiIterations',pre_jac, ...
    'PostJacobiIterations',post_jac, ...
    'MaxSize',max_size, ...
    'RelaxWeight',relax_weight, ...
    'RelaxMethod',relax_method, ...
    'Data',data);

  E = [E;max(abs(Aex*Z-Bex),[],1)/max(abs(Bex(:)))];
  if size(T,1) > 100
    for e = 1:size(E,2)
      set(s(e),'XData',1:size(E,1),'YData',E(:,e));
    end
    drawnow;
  end
  if E(end) < min_E
    %break;
  end
end
fprintf('MG: %g\n',toc);
