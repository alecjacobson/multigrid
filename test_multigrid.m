dim = 3;
method = 'fem';
addpath(genpath('~/Dropbox/SuiteSparse/CHOLMOD/'));

switch dim
case 2
  [V,T] = load_mesh('woody.obj');
  [V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  %[V,T] = upsample(V,T);
  V = V(:,1:2);
  F = [];
case 3
  %[V,F] = load_mesh('~/Documents/nested_cages/Meshes/Results/octopus-500k-disppath_final/octopus_0.obj');
  %%[V,F] = load_mesh('knight.ply');
  %%[V,F] = decimate_cgal(V,F,0.25);
  %%[V,F] = meshfix(V,F);
  %%[V,F] = load_mesh('~/Dropbox/models/armadillo.off');
  %%[V,F] = load_mesh('~/Dropbox/models/basic shapes/sphere-hires.obj');
  %h = mean(mean(edge_lengths(V,F)));
  %max_vol = h^3;
  %tic;
  %[V,T] = tetgen(V,F,'Flags',sprintf('-a%0.17f -q2',max_vol));
  %fprintf('Initial tetrahedralization: %gs %d tets\n',toc,size(T,1));
  [V,T,F] = readMESH('octopus.mesh');
  %writeMESH('knight.mesh',V,T,F);
  %writeMESH('knight-half.mesh',V,T,F);
  %[V,T,F] = readMESH('knight.mesh');
  %[V,T,F] = readMESH('knight-half.mesh');
  %medit(V,T,F);
  %writeMESH('armadillo.mesh',V,T,F);
  %[V,T,F] = readMESH('armadillo.mesh');
end

dim = size(T,2)-1;
M = massmatrix(V,T);

B = M*ones(size(V,1),1);
system_fun = @poisson_system;

%B = M*sin(V(:,1)/(max(V(:,1))-min(V(:,1)))*10);
%lambda = 1e4;
%system_fun = @(V,T,A,B,data) heat_system(V,T,A,B,lambda,data);

[Aex,Bex] = system_fun(V,T,[],B,[]);

switch method
case 'galerkin'
  Ain = Aex;
  B = Bex;
case 'fem'
  Ain = system_fun;
end

%% CHOLMOD form SuiteSparse toolbox
%tic;
%[cmL,p,cmS] = lchol(Aex);
%cmS = sparse(cmS,1:numel(cmS),1,size(Aex,1),size(Aex,1));
%fprintf('lchol precomp: %g\n',toc);
%tic;
%Zcm = cmS*(cmL'\(cmL\(cmS'*Bex)));
%fprintf('lchol solve: %g\n',toc);
%fprintf('lchol error to beat: %g\n',max(abs(Aex*Zcm-Bex))/max(abs(Bex)));

% MATLAB built ins
tic;
Zex = Aex\Bex;
min_E = max(abs(Aex*Zc-Bex))/max(abs(Bex));
fprintf('backslash: %g\n',toc);
%tic;
%[cL,p,cS] = chol(Aex);
%assert(p==0);
%fprintf('chol precomp: %g\n',toc);
%tic;
%Zc = cS*(cL\(cL'\(cS'*Bex)));
%fprintf('chol solve: %g\n',toc);
%min_E = max(abs(Aex*Zc-Bex))/max(abs(Bex));

fprintf('chol error to beat: %g\n',min_E);

relax_weight = 1.0;
relax_method = 'sor';
switch dim
case 2
  max_size = 100;
case 3
  %max_size = size(V,1)*0.006;
  max_size = 5000;
end

extrapolation = 'linear';
data = {};
% precomputation
tic;
[~,data] = multigrid( ...
  Ain,B, ...
  V,T, ...
  'Extrapolation',extrapolation, ...
  'RelaxWeight',relax_weight, ...
  'RelaxMethod',relax_method, ...
  'Data',data, ...
  'MaxSize',max_size, ...
  'PreJacobiIterations',1, ... % do at least one to set data
  'PostJacobiIterations',0, ...
  'BoundaryFacets',F, ...
  'Visualize',true);
fprintf('Pre-computation: %g, levels: %d\n',toc,numel(data));
save('mg-precompute.mat','data');

pre_jac = 2^0;
% This parameter should be tuned for the given example
post_jac = 2^1;

% Warm iteration
tic;
[~,data] = multigrid( ...
  Ain,B, ...
  [],[], ...
  'RelaxWeight',relax_weight, ...
  'RelaxMethod',relax_method, ...
  'Data',data, ...
  'MaxSize',max_size, ...
  'PreJacobiIterations',pre_jac, ...
  'PostJacobiIterations',post_jac, ...
  'Visualize',false);
fprintf('Warm iteration: %g\n',toc);
s = semilogy([nan nan;nan nan],[nan nan;nan nan],'-o','LineWidth',3);

Z = zeros(size(B));
%Z = B;
E = [];
tic;
while true

  [Z,data] = multigrid( ...
    Ain,B, ...
    [],[], ...
    'RelaxWeight',relax_weight, ...
    'RelaxMethod',relax_method, ...
    'Data',data, ...
    'MaxSize',max_size, ...
    'PreJacobiIterations',pre_jac, ...
    'PostJacobiIterations',post_jac, ...
    'Visualize',false, ...
    'Z0',Z);

  E = [E;max(abs(Aex*Z-Bex),[],1)/max(abs(Bex(:)))];
  if size(T,1) > 100
    for e = 1:size(E,2)
      set(s(e),'XData',1:size(E,1),'YData',E(:,e));
    end
    drawnow;
  end
  if E(end) < min_E
    break;
  end
end
fprintf('MG: %g\n',toc);
