function [z,M]=mgv(x,b,jac1,jac2,V,F)
  % MGV  Multigrid v-cycle.
  %
  % Adpated from
  % http://www.cs.berkeley.edu/~demmel/ma221_Fall09/Matlab/MG_README.html
  % 
  % z=mgv(x,b,jac1,jac2)
  %
  % Inputs:
  %   x  n by m initial guess, n and m are expected to be of the form 2^k+1
  %   b  n by m right hand side
  % Outputs:
  %   z  n by m solution
  %   M  n*m by n*m diagonal mass matrix
  %%   V  n*m by 2 mesh vertex positions
  %%   F  #F by 3 list of face indices into V
  %

  function B = vec(A)
    B = A(:);
  end
  function B = sqr(A)
    B = reshape(A,n,m);
  end
  function [V,F] = grid_mesh(n,m,irregular)
    if irregular
      max_area = (1/(n-1))*(1/(m-1))*2;
      [V,F] = triangle( ...
        [0 0;m 0;m n;0 n]/n,[1 2;2 3;3 4;4 1],[], ...
        'Flags',sprintf('-q32a%0.17f',max_area),'Quiet');
    else
      [V,F] = create_regular_grid(m,n,0,0);
      V = bsxfun(@times,[m n],V)/n;
    end
  end

  vis = true;
  global f;
  function append_gif()
    %filename = 'mgv.gif';
    %frame = getframe(gcf);
    %[SIf,cm] = rgb2ind(frame.cdata,256);
    %if f == 1
    %  imwrite(SIf,cm,filename,'Loop',Inf,'Delay',0);
    %else
    %  imwrite(SIf,cm, filename,'WriteMode','append','Delay',0);
    %end
    %f=f+1;
  end
  disc_mode = 'fem';
  irregular = true;

  switch disc_mode
  case 'fd'
    [n,m]=size(x);
    w=0.25;
    A = -fd_laplacian(n,m);
    h = 1/(n-1);
    M = diag(sparse(repmat(h^2,n*m,1)));
    bnd = boundary(n,m);
    nv = n*m;
    V = [];
    F = [];
  case 'fem'
    % FEM still converges with 0.25 but is a bit more accurate with 0.15
    if irregular
      w=0.20;
    else
      w=0.15;
    end
    if ~exist('V','var')
      [n,m]=size(x);
      % First mesh is always regular
      [V,F] = grid_mesh(n,m,false);
    else
      assert(all(x(:)==0),'x should be all zero');
      x = zeros(size(V,1),1);
    end
    nv = size(V,1);
    A = -cotmatrix(V,F);
    M = massmatrix(V,F,'voronoi');
    bnd = unique(outline(F));
  end
  A(bnd,:) = sparse(1:numel(bnd),bnd,1,numel(bnd),nv);
  b(bnd) = 0;

  % these signs are correct if we're solving A x = b

  if numel(bnd) > nv-numel(bnd)
  %  Solve small problem explicitly
     %z=zeros(3,3);
     %z(2,2)=b(2,2)/4;
     z = A\b(:);
  else
  %  Perform jac1 steps of weighted Jacobi with weight w
     % Demmel originally sets this to 2/3 then divides by 4. It appears that
     % 0.25 is safe.
     tb = [];
     tx = [];
     for i=1:jac1;
        %x(2:n-1,2:n-1)=(1-w)*x(2:n-1,2:n-1) ...
        % + (w/4)*( x(1:n-2,2:n-1) + x(3:n,2:n-1) ...
        %         + x(2:n-1,1:n-2) + x(2:n-1,3:n) ...
        %         + b(2:n-1,2:n-1));
        x = x(:)+w*(b(:)-A*x(:));
        if vis
          nsp = 2;
          switch disc_mode
          case 'fd'
            subplot(nsp,1,1);
            surf(sqr(b),'EdgeAlpha',n^-0.5);
            subplot(nsp,1,2);
            surf(sqr(x),'EdgeAlpha',n^-0.5);
          case 'fem'
            if isempty(tb)
              subplot(nsp,1,1);
              tb = tsurf(F,[V b(:)],fphong,'EdgeAlpha',nv^-0.25);
              subplot(nsp,1,2);
              tx = tsurf(F,[V x(:)],fphong,'EdgeAlpha',nv^-0.25);
            else
              tb.Vertices = [V b(:)];
              tx.Vertices = [V x(:)];
            end
          end
          subplot(nsp,1,1);
          fs = 15;
          title('Right-hand side','FontSize',fs);
          subplot(nsp,1,2);
          title('Solution','FontSize',fs);
          append_gif();
          view(90,0);
          drawnow;
        end
     end,
  %  Compute residual on current grid
     %r=zeros(n,n);
     %tmp = b(2:n-1,2:n-1) - ( 4*x(2:n-1,2:n-1) ...
     %    - x(1:n-2,2:n-1) - x(3:n,2:n-1) ...
     %    - x(2:n-1,1:n-2) - x(2:n-1,3:n) );
     %r(2:n-1,2:n-1) = tmp;
     % r is the _pointwise_ residual
     r = b(:) - A*x(:);

  %  Restrict residual to coarse grid
     %m=(n+1)/2;
     %rhat=zeros(m,m);
     %% Coarse grid perfectly corresponds to every-other fine grid node
     %rhat(2:m-1,2:m-1)=.25*r(3:2:n-2,3:2:n-2) ...
     %   + .125*(r(2:2:n-3,3:2:n-2)+r(4:2:n-1,3:2:n-2) ...
     %   +       r(3:2:n-2,2:2:n-3)+r(3:2:n-2,4:2:n-1) ) ...
     %   + .0625*(r(2:2:n-3,2:2:n-3)+r(4:2:n-1,2:2:n-3) ...
     %   +        r(2:2:n-3,4:2:n-1)+r(4:2:n-1,4:2:n-1) );
     
     switch disc_mode
     case 'fd'
       % Using 'nearest' here, does not seem to be such a big deal. Demmel
       % originally uses 'cubic'
       rhat = interp2(sqr(M\r),-1,'nearest');
       [nhat,mhat] = size(rhat);
       hhat = 1/(nhat-1);
       Mhat = diag(sparse(repmat(hhat^2,nhat*mhat,1)));
       rhat = Mhat*rhat(:);
       Vhat = [];
       Fhat = [];
     case 'fem'
       % Better be same V used in next level
         density = nv/(max(V(:,1))*max(V(:,2)));
         m = sqrt(density) * max(V(:,1));
         n = sqrt(density) * max(V(:,2));
       if irregular
         nhat = round(0.5*n);
         mhat = round(0.5*m);
       else
         nhat = round((n-1)/2+1);
         mhat = round((m-1)/2+1);
       end
       [Vhat,Fhat] = grid_mesh(nhat,mhat,irregular);
       Mhat = massmatrix(Vhat,Fhat,'voronoi');
       % It is **Very** important for irregular grids to extract the pointwise
       % solution **before** restricting and to compute the integrated solution
       % **before** solving on the coarse level.
       rhat = Mhat*prolong(V,F,M\r(:),Vhat);
     end

  %  Solve recursively with zero initial guess
  %  Multiply residual by 4 for coarser grid
     [xhat] = mgv(zeros(nhat,mhat),rhat,jac1,jac2,Vhat,Fhat);
     % xhat is now the solution to the pointwise problem on the coarse grid.
     % Compute integrated solution before prolongating.
     xhat = xhat(:);

  %  Interpolate coarse solution to fine grid, and add to x
     %xcor=zeros(size(x));
     %% Coarse grid perfectly corresponds
     %xcor(3:2:n-2,3:2:n-2)=xhat(2:m-1,2:m-1);
     %xcor(2:2:n-1,3:2:n-2)=.5*(xhat(1:m-1,2:m-1)+xhat(2:m,2:m-1));
     %xcor(3:2:n-2,2:2:n-1)=.5*(xhat(2:m-1,1:m-1)+xhat(2:m-1,2:m));
     %xcor(2:2:n-1,2:2:n-1)=.25*(xhat(1:m-1,1:m-1)+xhat(1:m-1,2:m)+ ...
     %                           xhat(2:m,1:m-1)+xhat(2:m,2:m) );
     % 'cubic' has little/no improvement over 'linear' {Demmel's default}, but
     % some improvement over 'nearest
     switch disc_mode
     case 'fd'
       xcor = vec(interp2(reshape(xhat,nhat,mhat),1,'linear'));
     case 'fem'
       xcor = prolong(Vhat,Fhat,xhat(:),V);
     end
     %save('mgv');
     % xcor is now the integrated solution, compute pointwise solution before
     % adding to current solution
     %tsurf(Fhat,[Vhat xhat],fphong);
     %input('low')
     %tsurf(F,[V xcor],fphong);
     %input('pointwise')
     xcor = xcor(:);
     %tsurf(F,[V xcor],fphong);
     %input('intregal')
     %tsurf(F,[V x],fphong);
     %input('x')

     z=x+xcor;
     %tsurf(F,[V z],fphong);
     %input('z')
  %  Perform jac2 steps of weighted Jacobi with weight w
     tb = [];
     tx = [];
     for i=1:jac2;
        %z(2:n-1,2:n-1)=(1-w)*z(2:n-1,2:n-1) ...
        % + (w/4)*( z(1:n-2,2:n-1) + z(3:n,2:n-1) ...
        %         + z(2:n-1,1:n-2) + z(2:n-1,3:n) ...
        %         + b(2:n-1,2:n-1));
        z = z+w*(-A*z(:)+b(:));
        if vis
          switch disc_mode
          case 'fd'
            subplot(nsp,1,1);
            surf(sqr(b),'EdgeAlpha',n^-0.5);
            subplot(nsp,1,2);
            surf(sqr(z),'EdgeAlpha',n^-0.5);
          case 'fem'
            if isempty(tb)
              subplot(nsp,1,1);
              tb = tsurf(F,[V b(:)],fphong,'EdgeAlpha',nv^-0.25);
              subplot(nsp,1,2);
              tx = tsurf(F,[V z(:)],fphong,'EdgeAlpha',nv^-0.25);
            else
              tb.Vertices = [V b(:)];
              tx.Vertices = [V z(:)];
            end
          end
          subplot(nsp,1,1);
          title('Right-hand side','FontSize',fs);
          subplot(nsp,1,2);
          title('Solution','FontSize',fs);
          append_gif();
          view(90,0);
          drawnow;
        end
     end,
  end
  norm(A*z(:)-b(:))
  return
  end
