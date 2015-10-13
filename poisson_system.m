function [A,B,data] = poisson_system(V,T,A,B,data)
  % POISSON_SYSTEM  Build a Poisson equation system with homogeneous 0
  % Dirichlet boundary conditions from a given mesh (V,T),
  % right-hand.
  %
  % [A,B,data] = poisson_system(V,T,A,B,data)
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   T  #T by dim+1 list of mesh element indices into V
  %   A  #V by #V system matrix (**with** boundary conditions, e.g. from
  %     previous call to this function) or []
  %   B  #V by ncols rhs (without boundary conditions)
  %   data  struct with precomputed data (see output)
  % Output:
  %   A  #V by #V system matrix **with** boundary conditions
  %   B  #V by ncols rhs adjust to account for boundary conditions
  %   data  struct containing precomputed data:
  %     .b  followed by list of indices to boundary vertices 
  %     .int  followed by list of indices to interior vertices
  %     .A_int_b  portion of original A corresponding to A(int,b)
  %
  if nargin<5 || isempty(data)
    switch size(T,2)
    case 4
      data.b = unique(boundary_faces(T));
    case 3
      data.b = unique(outline(T));
    end
    n = size(V,1);
    data.int = setdiff(1:n,data.b);
  end
  bc = zeros(numel(data.b),1);

  % DO NOT DO: A(b,:) = M(b,:) this will lead to ill-condition matrices as mesh
  % resollution increases and diag(M) ~ 0
  if isempty(A)
    L = cotmatrix(V,T);
    M = massmatrix(V,T);
    A = -L;
    %A = (L*(M\L));
    % Set known rows of system matrix to identity
    I = speye(n);
    A(data.b,:) = I(data.b,:);
    % Set known columns
    A(:,data.b) = I(:,data.b);
    data.A_int_b = A(data.int,data.b);
  end

  % Replace known values in rhs
  B(data.b,:) = bc;
  % Move contribution of known values to unknown portion of rhs
  if any(bc(:))
    B(data.int,:) = B(data.int,:)-data.A_int_b*bc;
  end
end
