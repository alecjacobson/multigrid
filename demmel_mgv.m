% Matlab code demmel_mgv.m
% For "Applied Numerical Linear Algebra",  Question 6.16
% Written by James Demmel, Jul 10, 1993
%                Modified, Jun  2, 1997
%
% Multigrid V-cycle for Poisson's equation on a square grid
% with zero Dirichlet boundary conditions.
% Algorithm from Brigg's ``Multigrid Tutorial''
% Include zero boundary values in arrays for ease of programming.
% Assume dimension n = 2^k + 1
% Inputs:
%   x = initial guess (n by n matrix with zeros on boundary)
%   b = right hand side (n by n matrix)
%   jac1, jac2 = number of weighted Jacobi steps to do before and
%                after recursive call to demmel_mgv
% Outputs:
%   z = improved solution (n by n matrix with zeros on boundary)
%
% function z=demmel_mgv(x,b,jac1,jac2)
%
function z=demmel_mgv(x,b,jac1,jac2)
[n,m]=size(b);
if n == 3,
%  Solve 1 by 1 problem explicitly
   z=zeros(3,3);
   z(2,2)=b(2,2)/4;
else
%  Perform jac1 steps of weighted Jacobi with weight w
   w=2/3;
   for i=1:jac1;
      x(2:n-1,2:n-1)=(1-w)*x(2:n-1,2:n-1) ...
       + (w/4)*( x(1:n-2,2:n-1) + x(3:n,2:n-1) ...
               + x(2:n-1,1:n-2) + x(2:n-1,3:n) ...
               + b(2:n-1,2:n-1));
      %nsp = 2;
      %subplot(nsp,1,1);
      %surf(b,'EdgeAlpha',1/n);
      %subplot(nsp,1,2);
      %surf(x,'EdgeAlpha',1/n);
      %drawnow;
   end,
%  Compute residual on current grid
   r=zeros(n,n);
   tmp = b(2:n-1,2:n-1) - ( 4*x(2:n-1,2:n-1) ...
       - x(1:n-2,2:n-1) - x(3:n,2:n-1) ...
       - x(2:n-1,1:n-2) - x(2:n-1,3:n) );
   r(2:n-1,2:n-1) = tmp;
%  Restrict residual to coarse grid
   m=(n+1)/2;
   rhat=zeros(m,m);
   rhat(2:m-1,2:m-1)=.25*r(3:2:n-2,3:2:n-2) ...
      + .125*(r(2:2:n-3,3:2:n-2)+r(4:2:n-1,3:2:n-2) ...
      +       r(3:2:n-2,2:2:n-3)+r(3:2:n-2,4:2:n-1) ) ...
      + .0625*(r(2:2:n-3,2:2:n-3)+r(4:2:n-1,2:2:n-3) ...
      +        r(2:2:n-3,4:2:n-1)+r(4:2:n-1,4:2:n-1) );
%  Solve recursively with zero initial guess
%  Multiply residual by 4 for coarser grid
   xhat = demmel_mgv(zeros(size(rhat)),4*rhat,jac1,jac2);
%  Interpolate coarse solution to fine grid, and add to x
   xcor=zeros(size(x));
   xcor(3:2:n-2,3:2:n-2)=xhat(2:m-1,2:m-1);
   xcor(2:2:n-1,3:2:n-2)=.5*(xhat(1:m-1,2:m-1)+xhat(2:m,2:m-1));
   xcor(3:2:n-2,2:2:n-1)=.5*(xhat(2:m-1,1:m-1)+xhat(2:m-1,2:m));
   xcor(2:2:n-1,2:2:n-1)=.25*(xhat(1:m-1,1:m-1)+xhat(1:m-1,2:m)+ ...
                              xhat(2:m,1:m-1)+xhat(2:m,2:m) );
   z=x+xcor;
%  Perform jac2 steps of weighted Jacobi with weight w
   for i=1:jac2;
      z(2:n-1,2:n-1)=(1-w)*z(2:n-1,2:n-1) ...
       + (w/4)*( z(1:n-2,2:n-1) + z(3:n,2:n-1) ...
               + z(2:n-1,1:n-2) + z(2:n-1,3:n) ...
               + b(2:n-1,2:n-1));
      %subplot(nsp,1,1);
      %surf(b,'EdgeAlpha',1/n);
      %subplot(nsp,1,2);
      %surf(z,'EdgeAlpha',1/n);
      %drawnow;
   end,
end
return
end
