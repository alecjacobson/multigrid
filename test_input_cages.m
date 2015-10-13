%prefix = '~/Documents/nested_cages/Meshes/Results/Model9_varap/Model9_varap_';
%prefix = '~/Documents/nested_cages/Meshes/Results/horse_varap/horse_varap_';
%prefix = '/Users/ajx/Documents/nested_cages/Meshes/Results/couplingdown_volume/couplingdown_';
prefix = '/Users/ajx/Documents/nested_cages/Meshes/Results/octopus_volume/octopus_';
%prefix = '/Users/ajx/Documents/nested_cages/Meshes/Results/pelvis_volume_new/pelvis_';
%prefix =  '/Users/ajx/Documents/nested_cages/Meshes/Results/arma_volume/arma_';
%prefix = '/Users/ajx/Documents/nested_cages/Meshes/Results/noisey_bunny-50k_varap/noisey_bunny-50k_varap_';
%prefix = '/Users/ajx/Documents/nested_cages/Meshes/Results/bunny-50k_varap/bunny-50k_varap_';
%prefix = '~/Documents/nested_cages/Meshes/Results/octopus-500k-disppath_final/octopus_';

%precompute = false;
precompute = true;

if precompute
  % load our cages
  cages = {};
  for iter = 0:100
    mesh = [];
    filename = sprintf('%s%d.obj',prefix,iter);
    if ~exist(filename,'file')
      break;
    end
    [mesh.V,mesh.F] = load_mesh(filename);
    cages{end+1} = mesh;
  end
  
  % create naive decimations to match
  naive = cell(numel(cages),1);
  naive{1}.V = cages{1}.V;
  naive{1}.F = cages{1}.F;
  for iter = 2:numel(cages)
    mesh = [];
    [mesh.V,mesh.F] = ...
      decimate_cgal( ...
        naive{iter-1}.V, ...
        naive{iter-1}.F, ...
        size(cages{iter}.F,1)/size(naive{iter-1}.F,1));
    [mesh.V,mesh.F] = meshfix(mesh.V,mesh.F);
    naive{iter} = mesh;
    clf;
    hold on;
    tsurf(cages{iter}.F,cages{iter}.V);
    tsurf(naive{iter}.F,bsxfun(@plus,naive{iter}.V, ...
      [1.2*(max(cages{iter}.V(:,1))-min(cages{iter}.V(:,1))) 0 0]));
    hold off;
    axis equal;
    drawnow;
    %pause;
  end
  
  % strip off first layer as input fine mesh
  V = cages{1}.V;
  F = cages{1}.F;
  cages = cages(2:end);
  naive = naive(2:end);

  [~,~,IF,~,~] = selfintersect(V,F,'DetectOnly',true,'FirstOnly',true);
  if ~isempty(IF)
    fprintf('meshfix''ing input...\n');
    [V,F] = meshfix(V,F);
  end
  
  h = mean(mean(edge_lengths(V,F)));
  max_vol = h^3;
  fprintf('tetgen''ing input...\n');
  [V,T] = tetgen(V,F,'Flags',sprintf('-a%0.17f -q2',max_vol));
  
  dim = size(T,2)-1;

  M = massmatrix(V,T);

  B = M*ones(size(V,1),1);
  system_fun = @poisson_system;
  
  %ease = @(t) 3*t.^2 - 2*t.^3;
  %B = M*ease(sin(V(:,1)/(max(V(:,1))-min(V(:,1)))*100)*0.5+0.5);
  %lambda = 1e-4;
  %[G,u,X,div_X,phi,pre,B,lambda] = heat_geodesic(V,T,c,[],'BoundaryConditions','neumann');
  %%B = div_X;
  %%lambda = 1e4;
  %warning('Using B from workspace...\n');
  %system_fun = @(V,T,A,B,data) heat_system(V,T,A,B,lambda,data);

  Ain = system_fun;
  [Aex,Bex] = system_fun(V,T,[],B,[]);
  
  fprintf('Direct solve...\n');
  Zex = Aex\Bex;
  
  min_E = max(abs(Aex*Zex-Bex))/max(abs(Bex));
  
  relax_weight = 1.0;
  relax_method = 'sor';
  max_size = -1; %indicate that we're passing layers
  extrapolation = 'linear';
  %extrapolation = 'constant';
  
  cages_data = {};
  fprintf('Precomputing cages mg...\n');
  % precomputation
  [~,cages_data] = multigrid( ...
    Ain,B, ...
    V,T, ...
    'Extrapolation',extrapolation, ...
    'Hierarchy',cages, ...
    'RelaxWeight',relax_weight, ...
    'RelaxMethod',relax_method, ...
    'Data',cages_data, ...
    'MaxSize',max_size, ...
    'PreJacobiIterations',1, ... % do at least one to set data
    'PostJacobiIterations',0, ...
    'BoundaryFacets',F);
  
  fprintf('Precomputing naive mg...\n');
  naive_data = {};
  % precomputation
  [~,naive_data] = multigrid( ...
    Ain,B, ...
    V,T, ...
    'Hierarchy',naive, ...
    'RelaxWeight',relax_weight, ...
    'RelaxMethod',relax_method, ...
    'Data',naive_data, ...
    'MaxSize',max_size, ...
    'PreJacobiIterations',1, ... % do at least one to set data
    'PostJacobiIterations',0, ...
    'BoundaryFacets',F);
end

pre_jac = 2^0;
% This parameter should be tuned for the given example
post_jac = 20;

s = semilogy([nan nan;nan nan],[nan nan;nan nan],'-o','LineWidth',3);

set(gca,'FontSize',30);
title('Residual Error','FontSize',30);
l = legend('Ours','Naive','FontSize',30);

tic;
Z = zeros(size(B,1),2);
E = [];
for iter = 1:30
  [Z(:,1),cages_data] = multigrid( ...
    Ain,B, ...
    [],[], ...
    'Z0',Z(:,1), ...
    'Hierarchy',cages, ...
    'RelaxWeight',relax_weight, ...
    'RelaxMethod',relax_method, ...
    'Data',cages_data, ...
    'MaxSize',max_size, ...
    'PreJacobiIterations',pre_jac, ... % do at least one to set data
    'PostJacobiIterations',post_jac);
  %if mod(iter,4) == 1
  %  medit(V,T,F,'Data',Z(:,1))
  %end

  [Z(:,2),naive_data] = multigrid( ...
    Ain,B, ...
    [],[], ...
    'Z0',Z(:,2), ...
    'Hierarchy',naive, ...
    'RelaxWeight',relax_weight, ...
    'RelaxMethod',relax_method, ...
    'Data',naive_data, ...
    'MaxSize',max_size, ...
    'PreJacobiIterations',pre_jac, ... % do at least one to set data
    'PostJacobiIterations',post_jac);

  E = [E;max(abs(bsxfun(@minus,Aex*Z,Bex)),[],1)/max(abs(Bex(:)))];
  if size(T,1) > 10000
    for e = 1:size(E,2)
      set(s(e),'XData',1:size(E,1),'YData',E(:,e));
    end
    drawnow;
  end
  %if (min(E(end,:)) < min_E && max(E(end,:))>1) || max(E(end,:))<min_E
  %  break;
  %end

end
fprintf('mg solve: %g\n',toc);
