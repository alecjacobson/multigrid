function Z = vcycle_rec(B,level)
  mu = 20;
  omega = 0.25;
  [V,F] = create_regular_grid(2^level+1,2^level+1,0,0);
  [CV,CF] = create_regular_grid(2^(level-1)+1,2^(level-1)+1,0,0);
  n = size(V,1);
  A = -cotmatrix(V,F);
  M = massmatrix(V,F,'voronoi');
  if isempty(B)
    B = M*ones(n,1);
  end
  b = unique(outline(F));
  A(b,:) = sparse(1:numel(b),b,1,numel(b),n);
  B(b) = 0;

  if level == 1
    w = 1e9;
    Z = (A+sparse(b,b,w,n,n)) \ B;
    return;
  end

  % presmooth
  Richardson = @(B,X) X+omega*(B-A*X);
  Z = zeros(n,1);
  for smooth_step = 1:mu
    Z = Richardson(B,Z);
    Z(b) = 0;
    if mod(smooth_step,10)==1 && level>4
      tsurf(F,[V Z],'EdgeAlpha',1/sqrt(size(V,1)));
      drawnow;
    end
  end

  cycles = 5;
  for cycle = 1:cycles
    % Restrict
    residual = B-A*Z;
    Cresidual = prolong(V,F,residual,CV);

    CY = vcycle_rec(Cresidual,level-1);

    Z = Z + prolong(CV,CF,CY,V);
    % post smooth
    for smooth_step = 1:mu
      Z = Richardson(B,Z);
      Z(b) = 0;
      if mod(smooth_step,10)==1 && level>4
        tsurf(F,[V Z],'EdgeAlpha',1/sqrt(size(V,1)));
        drawnow;
      end
    end
  end

end
