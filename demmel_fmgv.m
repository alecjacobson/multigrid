% Matlab code demmel_fmgv.m
% For "Applied Numerical Linear Algebra",  Question 6.16
% Written by James Demmel, Jul 10, 1993
%                Modified, Jun  2, 1997
%
% Full Multigrid V-Cycle for Poisson's equation on a square grid
% with zero Dirichlet boundary conditions.
% Algorithm from Brigg's ``Multigrid Tutorial''
% Include zero boundary values in arrays for ease of programming.
% Assume dimension n = 2^k + 1
% Inputs:
%   x = initial guess (n by n matrix with zeros on boundary)
%   b = right hand side (n by n matrix)
%   jac1, jac2 = number of weighted Jacobi steps to do before and
%                after recursive call to mgv
% Outputs:
%   z = improved solution (n by n matrix with zeros on boundary)
%
% function z=demmel_fmgv(x,b,jac1,jac2);
%
function z=demmel_fmgv(x,b,jac1,jac2);
[nn,m]=size(b);
%  Compute residual
   r=zeros(nn,nn);
   tmp = b(2:nn-1,2:nn-1) - ( 4*x(2:nn-1,2:nn-1) ...
       - x(1:nn-2,2:nn-1) - x(3:nn,2:nn-1) ...
       - x(2:nn-1,1:nn-2) - x(2:nn-1,3:nn) );
   r(2:nn-1,2:nn-1) = tmp;
% Get right hand side for coarsest grid by repeated restriction of
% original right hand side to coarser grids
rc=demmel_mgvrhs(nn,3,r);
% Solve 1 by 1 problem
z=zeros(3,3);
z(2,2)=rc(2,2)/4;
%  Assume nn = 2^k+1
k = round(log(nn-1)/log(2)); 
% Loop from the next-to-coarsest to finest grids
for i=2:k,
   n=2^i+1;
   m=2^(i-1)+1;
%  Interpolate coarser solution to this level
   zstrt=zeros(n,n);
   zstrt(3:2:n-2,3:2:n-2)=z(2:m-1,2:m-1);
   zstrt(2:2:n-1,3:2:n-2)=.5*(z(1:m-1,2:m-1)+z(2:m,2:m-1));
   zstrt(3:2:n-2,2:2:n-1)=.5*(z(2:m-1,1:m-1)+z(2:m-1,2:m));
   zstrt(2:2:n-1,2:2:n-1)=.25*(z(1:m-1,1:m-1)+z(1:m-1,2:m)+ ...
                              z(2:m,1:m-1)+z(2:m,2:m) );
% Get right hand side for current grid by repeated restriction of
% original right hand side to coarser grids
   rhs=demmel_mgvrhs(nn,n,r);
%  Do Multigrid V-cycle
   z=mgv(zstrt,rhs,jac1,jac2);
%  mesh(z),title(int2str(n)), pause
end
% Add final correction to initial guess
z=x+z;
return
