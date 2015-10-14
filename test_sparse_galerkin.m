[V,T] = regular_tetrahedral_mesh(20);

L = cotmatrix(V,T);
M = massmatrix(V,T);

A = -L;
b = unique(boundary_faces(T));
A(boundary,:) = I(boundary,:);
A(:,boundary) = I(:,boundary);
B = M*ones(size(V,1),1);
B(boundary) = 0;
Zex = A\B;

F = boundary_faces(T);


fprintf('Precomputing naive mg...\n');
naive_data = {};
% precomputation
[~,naive_data] = multigrid( ...
  A,B, ...
  V,T, ...
  'Data',naive_data, ...
  'RelaxMethod','sor', ...
  'RelaxWeight',1.0, ...
  'PreJacobiIterations',1, ... % do at least one to set data
  'PostJacobiIterations',0, ...
  'BoundaryFacets',F);

s = semilogy([nan ;nan],[nan ;nan],'-o','LineWidth',3);
set(gca,'FontSize',15);

tic;
Z = zeros(size(B,1),1);
E = [];
for iter = 1:100

  [Z,naive_data] = multigrid( ...
    A,B, ...
    [],[], ...
    'Z0',Z, ...
    'PreJacobiIterations',10, ...
    'PostJacobiIterations',10, ...
    'RelaxMethod','sor', ...
    'RelaxWeight',1.0, ...
    'Data',naive_data);

  E = [E;max(abs(bsxfun(@minus,A*Z,B)),[],1)/max(abs(B(:)))];
  if size(T,1) > 10000
    for e = 1:size(E,2)
      set(s(e),'XData',1:size(E,1),'YData',E(:,e));
    end
    drawnow;
  end

end
