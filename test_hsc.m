A = adjacency_matrix(T);
A = diag(sparse(sum(A,2)))-A;
b = unique(boundary_faces(T));
I = speye(size(A,1));
A(b,:) = I(b,:);
A(:,b) = I(:,b);
int = setdiff(1:size(V,1),b);
B = zeros(size(V,1),1);
B(int) = 1- A(int,b)*zeros(numel(b),1);
tic;
Z = A\B;
toc
hsc_fun = hsc_setup(A,A);
tic;
Zpcg = pcg(A,B, 1e-14, 100, hsc_fun, []);
toc
