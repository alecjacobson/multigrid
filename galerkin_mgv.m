function [Z,data] = galerkin_mgv(A,B,V,T,varargin)
  % GALERKIN_MGV Solve A Z = B using Ritz-Galerkin style multigrid.
  %
  % Z = galerkin_mgv(A,B,V,T)
  % [Z,data] = ...
  %   galerkin_mgv(A,B,V,T,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   A  #V by #V system matrix
  %   B  #V by 1 right hand side vector
  %   V  #V by dim list of vertex positions at finest level
  %   T  #T by dim+1 list of element indices into V
  %   Optional:
  %     'PreJacobiIterations' followed by number of pre-jacobi iterations {4}
  %     'PostJacobiIterations' followed by number of post-jacobi iterations.
  %       Seems this should be set greater than PreJacobiIterations for best
  %       performance {10}
  %     'BoundaryFacets'  followed by #F by 3 list of surface triangles. These
  %       triangles could be super-triangles of the actually boundary of T
  %       (only makes sense for dim==3)
  %     'RelaxWeight' followed by Jacobi relaxing weight {0.2}
  %     'MaxSize'  followed by maximum subproblem size before using direct
  %       solver {1000}
  %     'Visualize'  followed by whether to visualize progress {false}
  %     'Z0'  followed by #V by 1 intial guess (warm start) {[] --> zeros}
  %     'Data'  #levels of precomputed data (see output) {[]}
  % Outputs:
  %   Z  #V by 1 solution vector
  %   data  #levels of precomputed data: meshes, prolongation operators,
  %     reduced system matrices, relaxing precomputation
  %
  % See also: prolongation
  %    

  % Would be interesting to try
  % THE MULTIGRID PRECONDITIONED CONJUGATE GRADIENT METHOD
  % [Osamu Tatebe 1993] 

  % default values
  pre_jac = 4;
  post_jac = 10;
  w = 0.20;
  max_n = 1000;
  Z = [];
  vis = false;
  tetgen_flags = '-q10';
  restrict_scalar = 1;
  relax_method = 'jacobi';
  BF = [];
  data = {};

  rec_params = { ...
    'PreJacobiIterations','PostJacobiIterations','RelaxWeight', ...
    'Visualize','MaxSize','TetgenFlags','RestrictScalar','RelaxMethod'};
  rec_vars = { ...
    'pre_jac','post_jac','w', ...
    'vis','max_n','tetgen_flags','restrict_scalar','relax_method'};
  init_params = {'Z0','Data','BoundaryFacets'};
  init_vars = {'Z','data','BF'};
  params = {rec_params{:},init_params{:}};
  vars = {rec_vars{:},init_vars{:}};
  % Map of parameter names to variable names
  params_to_variables = containers.Map(params,vars);
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this
      % workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  if isempty(data)
    data = {[]};
  end
  assert(numel(data)>=1);

  if isempty(V)
    assert(isempty(T),'T should also be empty');
    V = data{1}.V;
    T = data{1}.T;
  else
    data{1}.V = V;
    data{1}.T = T;
  end

  dim = size(T,2)-1;
  switch dim
  case 2
    if size(V,2)==3
      warning('Only considering xy coordinates of triangle mesh');
      % must be 2d
      V = V(:,1:2);
    end
    assert(size(V,2) == dim);
  case 3
    assert(size(V,2) == dim);
  otherwise
    error('Unsupported dimension (%d)',dim);
  end


  % number of nodes
  n = size(A,1);
  % number of right-hand sides
  ncols = size(B,2);
  assert(n == size(B,1),'B should match A');
  assert(n == size(V,1),'V should match A');
  % Initial guess
  if isempty(Z)
    Z = zeros(n,ncols);
  end
  % number of boundary nodes
  %b = unique(outline(T));
  %nb = numel(b);
  %if nb > (n-nb)

  % Rescaling does seem to help, but I'm not convinced it's worth it
  if ~isfield(data{1},'maxA') || isempty(data{1}.maxA)
    data{1}.maxA = max(abs(A(:)));
  end
  %B = B/data{1}.maxA;
  %A = A/data{1}.maxA;

  if n < max_n
    %[cond(full(A))]
    Z = A \ B;
    return;
  end

  % prerelax
  [Z,data{1}] = relax(A,B,Z,pre_jac, ...
    'Method',relax_method,'Data',data{1},'Weight',w);

  % Build coarse mesh
  CF = [];
  if numel(data)<2
    data{2} = [];
    [data{2}.V,data{2}.T,CF] = coarsen( ...
      V,T, ...
      'BoundaryFacets',BF, ...
      'TetgenFlags',tetgen_flags);
  end

  % Build prolongation operator
  if ...
    ~isfield(data{1},'P') || isempty(data{1}.P) || ...
    ~isfield(data{1},'R') || isempty(data{1}.R)
    data{1}.P = prolongation(data{2}.V,data{2}.T,V);
    data{1}.R = restrict_scalar*data{1}.P';
  end
  % Restriction
  %R = prolongation(V,T,VC);
  %D = diag(sparse(sqrt(diag(R*P))));
  %P = P*D;
  %R = P';

  % https://computation-rnd.llnl.gov/linear_solvers/pubs/nongalerkin-2013.pdf
  % page 2
  if ~isfield(data{1},'RAP') || isempty(data{1}.RAP)
    data{1}.RAP = data{1}.R*A*data{1}.P;
  end

  % compute residual
  res = B-A*Z;
  Rres = data{1}.R*res;

  rec_vals = cellfun( ...
    @(name)evalin('caller',name),rec_vars,'UniformOutput',false);
  rec_params_vals = reshape({rec_params{:};rec_vals{:}},1,numel(rec_params)*2);
  [CZ,data_child] = ...
    galerkin_mgv( ...
    data{1}.RAP, ...
    Rres, ...
    [], ...
    [], ...
    'Data',data(2:end), ...
    'BoundaryFacets', CF, ...
    rec_params_vals{:});
  data(1+(1:numel(data_child))) = data_child;

  % Add correction
  Z = Z+data{1}.P*CZ;

  % postrelax
  [Z,data{1}] = relax(A,B,Z,post_jac, ...
    'Method',relax_method,'Data',data{1},'Weight',w);
end
