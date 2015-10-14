%[V,T] = regular_tetrahedral_mesh(13);

L = cotmatrix(V,T);
M = massmatrix(V,T);

A = -L;
I = speye(size(A));
boundary = unique(boundary_faces(T));
A(boundary,:) = I(boundary,:);
A(:,boundary) = I(:,boundary);
B = M*ones(size(V,1),1);
B(boundary) = 0;
Zex = A\B;

F = boundary_faces(T);

% relaxation weight
rw = 1.0;
sparse_P = 1.0;

algebraic = true;
if algebraic
  effV = [];
  effT = [];
else
  effV = V;
  effT = T;
end

naive_data = {};
% precomputation
tic;
[~,naive_data] = multigrid( ...
  A,B, ...
  effV,effT, ...
  'Algebraic',algebraic, ...
  'SparseP',sparse_P, ...
  'Data',naive_data, ...
  'RelaxMethod','sor', ...
  'RelaxWeight',rw, ...
  'PreJacobiIterations',1, ... % do at least one to set data
  'PostJacobiIterations',0, ...
  'BoundaryFacets',F);
fprintf('Precomputation: %g secs\n',toc);

s = semilogy([nan ;nan],[nan ;nan],'-o','LineWidth',3);
set(gca,'FontSize',15);

tic;
Z = zeros(size(B,1),1);
E = [];
tic;
for iter = 1:100

  [Z,naive_data] = multigrid( ...
    A,B, ...
    [],[], ...
    'Z0',Z, ...
    'Algebraic',algebraic, ...
    'SparseP',sparse_P, ...
    'PreJacobiIterations',2, ...
    'PostJacobiIterations',2, ...
    'RelaxMethod','sor', ...
    'RelaxWeight',rw, ...
    'Data',naive_data);

  E = [E;max(abs(bsxfun(@minus,A*Z,B)),[],1)/max(abs(B(:)))];
  if size(T,1) > 10000
    for e = 1:size(E,2)
      set(s(e),'XData',1:size(E,1),'YData',E(:,e));
    end
    drawnow;
  end
  if E(end)<1e-13
    break;
  end

end
fprintf('Converged in: %g secs\n',toc);
