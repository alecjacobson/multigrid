function [CV,CT,CF] = coarsen(V,T,varargin)
  % COARSEN Coarse a given volumetric mesh (V,T)
  %
  % [CV,CT,CF] = coarsen(V,T)
  % [CV,CT,CF] = coarsen(V,T,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by dim list of fine mesh vertex positions
  %   T  #T by dim+1 list of fine mesh element indices
  %   Optional:
  %     'TetgenFlags'  followed by tetgen flags {'-q2'}
  %     'BoundaryFacets'  followed by #BF by dim list of boundary facets
  %       indexing V (not necessary the boundary of the volumetric mesh T, but
  %       should probably be a superset of those facets)
  %     'CV'  followed by #CV by dim list of vertex positions of boundary
  %       vertices for coarse mesh {computed using uniformly_sample in 2D or
  %       decimate_cgal in 3d}
  %     'CF'  followed by #CF by dim list of boundary mesh facets indexing CV
  % Outputs:
  %   CV  #CV by dim list of output mesh vertices: input CV will always come
  %     first
  %   CT  #CT by dim+1 list of output mesh element indices into CV
  %   CF  #CF by dim list of boundary facets defining (at least as a superset)
  %     the boundary of (CV,CT). Will coincide with input CF if it was
  %     non-empty.
  %      
  
  tetgen_flags = '-q2';
  BF = [];
  CV = [];
  CF = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'BoundaryFacets','TetgenFlags','CV','CF'}, ...
    {'BF','tetgen_flags','CV','CF'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  dim = size(T,2)-1;
  % Generate mesh at next level
  switch dim
  case 3 
    if isempty(CV) || isempty(CF)
      if isempty(BF)
        BF = boundary_faces(T);
      end
      [BV,IM] = remove_unreferenced(V,BF);
      BF = IM(BF);
      ratio = (2^(-2/3));
      [CV,CF] = decimate_cgal(BV,BF,ratio);
      [~,~,IF] = selfintersect(CV,CF,'FirstOnly',true,'DetectOnly',true);
      O = outline(CF);
      if ~isempty(O) || ~isempty(IF)
        warning('meshfixing');
        [CV,CF] = meshfix(CV,CF);
      end
   end
   max_vol = mean(volume(V,T))*8;
   area_flags = sprintf('-a%0.17f ',max_vol);
   [CV,CT,FF] = tetgen(CV,CF,'Flags',[area_flags tetgen_flags]);
   %medit(CV,CT,CF);
  case 2 
    Vloop = V(outline_loop(T),:);
    s = size(Vloop,1);
    sc = ceil(s/sqrt(2));
    VCloop = uniformly_sample(Vloop,sc,'Loop',true);
    CF = [1:size(VCloop,1);2:size(VCloop,1) 1]';
    %hC = mean(edge_lengths(VCloop,CF));
    %max_area = hC^2;
    max_area = mean(doublearea(V,T)/2)*4;
    [CV,CT] = triangle(VCloop,CF,[], ...
      'Flags',sprintf('-q32a%0.17f',max_area),'Quiet');
    CF = []; % not supported for 2d
  end
end
