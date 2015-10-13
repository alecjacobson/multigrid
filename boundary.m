function b = boundary(n,m)
  % BOUNDARY boundary indices on a regular grid
  %
  % Inputs:
  %   n  number of rows (y)
  %   m  number of columns (x)
  % Outputs
  %   b  2*n+2*m-4 list of boundary indices
  %
  b = [1:n n+1:n:(m-2)*n+1 2*n:n:(m-1)*n n*m-n+1:n*m];
end
