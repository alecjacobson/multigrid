function [Z,data] = fem_mgv(system_fun,B,V,T,varargin)
  % GALERKIN_MGV Solve A Z = B using FEM-style multigrid.
  %
  % Z = fem_mgv(A,B,V,T)
  % [Z,data] = ...
  %   fem_mgv(A,B,V,T,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   system_fun  handle to function computing [A,B] given a mesh (at any
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

  % default parameters
  pre_jac = 4;
  post_jac = 10;
  w = 0.2;
  vis = false;
  relax_method = 'jacobi';
  tetgen_flags = '-q2';
  restrict_scalar = 1;
  BF = [];
  Z = [];
  data = [];
  max_n = 100;

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

  % Construct system
  [A,B] = system_fun(V,T,[],B);
  % number of nodes
  n = size(A,1);
  dim = size(T,2)-1;
  % number of right-hand sides
  ncols = size(B,2);
  assert(n == size(B,1),'B should match A');
  assert(n == size(V,1),'V should match A');
  % Initial guess
  if isempty(Z)
    Z = zeros(n,ncols);
  end

  % Base case
  if n < max_n
    %[cond(full(A))]
    Z = A \ B;
    return;
  end

  % Pre-smooth 
  Z = relax(A,B,Z,pre_jac, ...
    'Method',relax_method,'Weight',w,'V',V,'T',T);

  % Compute _pointwise_ residual
  res = B-A*Z;

  % Build coarse mesh
  [CV,CT,CF] = coarsen( ...
    V,T, ...
    'BoundaryFacets',BF, ...
    'TetgenFlags',tetgen_flags);
  %[size(T,1) size(CT,1) size(T,1)/size(CT,1)]

  % Build linear prolongation operator
  P = prolongation(CV,CT,V);
  R = P';

  Rres = R*res;


  rec_vals = cellfun( ...
    @(name)evalin('caller',name),rec_vars,'UniformOutput',false);
  rec_params_vals = reshape({rec_params{:};rec_vals{:}},1,numel(rec_params)*2);
  % Solve for residual on coarse mesh recursively
  CZ = fem_mgv(system_fun,Rres,CV,CT, ...
    'BoundaryFacets',CF, ...
    rec_params_vals{:});

  % Prolong _pointwise_ residual from coarse mesh
  Z = Z+P*CZ;
  % Post-smooth
  Z = relax(A,B,Z,post_jac, ...
    'Method',relax_method,'Weight',w);
end
