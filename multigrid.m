function [Z,data] = multigrid(Ain,B,V,T,varargin)
  % MULTIGRID  Solve A Z = B using either FEM-style (re-discretization)
  % multigrid or Galerkin (restrict-the-problem) style.
  %
  % Z = multigrid(A,B,V,T) % galerkin
  % Z = multigrid(A_fun,B,V,T) % fem
  % [Z,data] = ...
  %   multigrid(A,B,V,T,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   A  #V by #V system matrix
  %     or
  %   A_fun  handle to function computing [A,B] given a mesh (at any
  %     resolution (V,T) and a _vanilla_ right-hand side B
  %   B  #V by 1 initial right-hand side vector
  %   V  #V by dim list of vertex positions at finest level
  %   T  #T by dim+1 list of element indices into V
  %   Optional:
  %     'PreJacobiIterations' followed by number of pre-jacobi iterations {4}
  %     'PostJacobiIterations' followed by number of post-jacobi iterations.
  %       Seems this should be set greater than PreJacobiIterations for best
  %       performance {10}
  %     'BoundaryFacets'  followed by #F by 3 list of surface triangles. These
  %       triangles could be super-triangles of the actually boundary of T (only
  %       makes sense for dim==3)
  %     'Hierarchy' followed by #levels list of {CV,CF} pairs of surface meshes,
  %     'MaxSize' should be set the size of the last CF.
  %     'RelaxWeight' followed by Jacobi relaxing weight {0.2}
  %     'MaxSize'  followed by maximum subproblem size before using direct
  %       solver {1000}
  %     'Visualize'  followed by whether to visualize progress {false}
  %     'Z0'  followed by #V by 1 intial guess (warm start) {[] --> zeros}
  %     'Data'  #levels of precomputed data (see output) {[]}
  %     'Extrapolation'  followed by extrapolation method used for prolongation
  %       operator (see prolongation.m optional input)
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
  tetgen_flags = '-q2';
  restrict_scalar = 1;
  relax_method = 'jacobi';
  BF = [];
  levels = {};
  data = {};
  extrapolation = 'linear';

  rec_params = { ...
    'PreJacobiIterations','PostJacobiIterations','RelaxWeight', ...
    'Visualize','MaxSize','TetgenFlags','RestrictScalar','RelaxMethod', ...
    'Extrapolation'};
  rec_vars = { ...
    'pre_jac','post_jac','w', ...
    'vis','max_n','tetgen_flags','restrict_scalar','relax_method', ...
    'extrapolation'};
  init_params = {'Z0','Data','BoundaryFacets','Hierarchy'};
  init_vars = {'Z','data','BF','levels'};
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

  if ~isempty(levels)
    max_n = -1;
  end

  if isempty(data)
    data = {[]};
  end

  if isnumeric(Ain)
    A = Ain;
    method = 'galerkin';
  else
    assert(isa(Ain,'function_handle'));
    method = 'fem';
    if ~isfield(data{1},'A') 
      data{1}.A = [];
      data{1}.A_fun_data = [];
    end
    [data{1}.A,B,data{1}.A_fun_data] = ...
      Ain( ...
        V,T, ...
        data{1}.A,B,data{1}.A_fun_data);
    A = data{1}.A;
  end
  
  if isempty(data)
    data = {[]};
  end
  assert(numel(data)>=1);

  n = size(A,1);

  % number of right-hand sides
  ncols = size(B,2);
  assert(n == size(B,1),'B should match A');
  assert(isempty(V) || n == size(V,1),'V should match A');
  % Initial guess
  if isempty(Z)
    Z = zeros(n,ncols);
  end

  if n <= max_n || (max_n == -1 && numel(data) == 1 && numel(levels) == 0)
    %[cond(full(A))]
    Z = A \ B;
    return;
  end

  % prerelax
  [Z,data{1}] = relax(A,B,Z,pre_jac, ...
    'Method',relax_method,'Data',data{1},'Weight',w);

  % Build coarse mesh
  if isempty(levels)
    CV = [];
    CF = [];
  else
    CV = levels{1}.V;
    CF = levels{1}.F;
  end
  if numel(data)<2
    data{2} = [];
    [CV,CT,CF] = coarsen( ...
      V,T, ...
      'CV',CV, ...
      'CF',CF, ...
      'BoundaryFacets',BF, ...
      'TetgenFlags',tetgen_flags);
  else
    CV = [];
    CT = [];
  end

  % Build prolongation operator
  if ...
    ~isfield(data{1},'P') || isempty(data{1}.P) || ...
    ~isfield(data{1},'R') || isempty(data{1}.R)
    data{1}.P = prolongation(CV,CT,V,'Extrapolation',extrapolation);
    data{1}.R = restrict_scalar*data{1}.P';
  end
  % Restriction
  %R = prolongation(V,T,VC);
  %D = diag(sparse(sqrt(diag(R*P))));
  %P = P*D;
  %R = P';

  % https://computation-rnd.llnl.gov/linear_solvers/pubs/nongalerkin-2013.pdf
  % page 2
  switch method
  case 'galerkin'
    if ~isfield(data{1},'RAP') || isempty(data{1}.RAP)
      data{1}.RAP = data{1}.R*A*data{1}.P;
    end
    CAin = data{1}.RAP;
  case 'fem'
    CAin = Ain;
  end

  % compute residual
  res = B-A*Z;
  Rres = data{1}.R*res;

  % Trick to pass on parameters
  rec_vals = cellfun( ...
    @(name)evalin('caller',name),rec_vars,'UniformOutput',false);
  rec_params_vals = reshape({rec_params{:};rec_vals{:}},1,numel(rec_params)*2);
  % Solve for residual on coarse mesh recursively
  [CZ,data_child] = ...
    multigrid( ...
    CAin, ...
    Rres, ...
    CV, ...
    CT, ...
    'Hierarchy', levels(2:end), ...
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
