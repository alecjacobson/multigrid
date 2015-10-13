function [A,B,data] = heat_system(V,T,A,B,lambda,data)
  % HEAT_SYSTEM  Build a Heat equation system with homogeneous 0 Neumann
  % boundary conditions from a given mesh (V,T), right-hand.
  %
  % [A,B,data] = heat_system(V,T,A,B,data)
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

  % DO NOT DO: A(b,:) = M(b,:) this will lead to ill-condition matrices as mesh
  % resollution increases and diag(M) ~ 0
  if isempty(A)
    L = cotmatrix(V,T);
    M = massmatrix(V,T);
    A = -lambda*L+M;
    B = B;
  end

end
