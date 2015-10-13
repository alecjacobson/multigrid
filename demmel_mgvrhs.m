% Matlab code demmel_mgvrhs.m
% For "Applied Numerical Linear Algebra",  Question 6.16
% Written by James Demmel, Jul 10, 1993
%                Modified, Jun  2, 1997
%
% Restrict a matrix from grid n to grid m, where n and m 
% are (powers of 2) minus 1
% used to get right-hand-side for multigrid solve
%
% function nrhs=demmel_mgvrhs(n,m,rhs)
% 
function nrhs=demmel_mgvrhs(n,m,rhs)
nk=round(log(n-1)/log(2));
mk=round(log(m-1)/log(2));
nrhs=rhs;
r=rhs;
for i=nk-1:-1:mk,
    mi=2^i+1;
    ni=2^(i+1)+1;
    nrhs=zeros(mi,mi);
    nrhs(2:mi-1,2:mi-1)=.25*r(3:2:ni-2,3:2:ni-2) ...
      + .125*(r(2:2:ni-3,3:2:ni-2)+r(4:2:ni-1,3:2:ni-2) ...
      +       r(3:2:ni-2,2:2:ni-3)+r(3:2:ni-2,4:2:ni-1) ) ...
      + .0625*(r(2:2:ni-3,2:2:ni-3)+r(4:2:ni-1,2:2:ni-3) ...
      +        r(2:2:ni-3,4:2:ni-1)+r(4:2:ni-1,4:2:ni-1) );
    r=nrhs;
end
return
