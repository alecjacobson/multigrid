%function Z = vcycle(mu)
  mu = 1000;
  omega = 0.25;
  levels = 6;
  residual = cell(levels,1);
  % Richardson iteration
  Richard = cell(levels,1);
  Z = cell(levels,1);
  A = cell(levels,1);
  B = cell(levels,1);
  V = cell(levels,1);
  F = cell(levels,1);
  for level = 1:levels
    [V{level},F{level}] = create_regular_grid(2^level+1,2^level+1,0,0);
  end

  for level = levels:-1:1
    n = size(V{level},1);
    A{level} = -cotmatrix(V{level},F{level});
    M = massmatrix(V{level},F{level},'voronoi');
    B{level} = M*ones(n,1);
    b = unique(outline(F{level}));
    if level == 1
      Z{level} = (A{level}+sparse(b,b,w,n,n)) \ ...
        (residual{level});
    else
      Richard{level} = @(B,X) X+omega*(B-A{level}*X);
      if level == levels
        residual{level} = B{level};
      end
      Z{level} = zeros(n,1);
      for smooth_step = 1:mu
        Z{level} = Richard{level}(residual{level},Z{level});
        Z{level}(b) = 0;
        tsurf(F{level},[V{level} Z{level}],'EdgeAlpha',1/sqrt(size(V{level},1)));
        if mod(smooth_step,10)==1
          drawnow;
        end
      end
      % Restrict
      residual{level-1} = prolong(V{level},F{level}, ...
        residual{level} - A{level}*Z{level}, ...
        V{level-1});
    end
  end


  for level = 2:levels
    Z{level} = Z{level} + prolong(V{level-1},F{level-1},Z{level-1},V{level});
    b = unique(outline(F{level}));
    for smooth_step = 1:mu
      Z{level}(b) = 0;
      Z{level} = Richard{level}(residual{level},Z{level});
      Z{level}(b) = 0;
      tsurf(F{level},[V{level} Z{level}],'EdgeAlpha',1/sqrt(size(V{level},1)));
      if mod(smooth_step,10)==1
        drawnow;
      end
    end
  end

%end
