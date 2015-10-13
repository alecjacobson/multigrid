function [V,T,F] = irregular_tetrahedral_mesh(nx)
  fprintf('regular_tetrahedral_mesh\n');
  [V,T,F] = regular_tetrahedral_mesh(nx);
  in = setdiff(1:size(V,1),F);
  h = 1/(nx-1);
  d = 0.9;
  V(in,:) = V(in,:) + h*d*(2*rand(numel(in),3)-1);
  fprintf('tetgen\n');
  [V,T] = tetgen(V,F,'Flags','-Y');
  %[V,F] = create_regular_grid(nx,nx);
  %V(:,3) = 0;
  %V = bsxfun(@minus,V,[0.5 0.5 0.5]);
  %n = size(V,1);
  %V = [ ...
  %  V*axisangle2matrix([0 1 0],pi/2); ...
  %  V*axisangle2matrix([0 1 0],0); ...
  %  V*axisangle2matrix([0 1 0],pi); ...
  %  V*axisangle2matrix([0 1 0],-pi/2); ...
  %  V*axisangle2matrix([1 0 0],pi/2); ...
  %  V*axisangle2matrix([1 0 0],-pi/2); ...
  %  ];
  %F = [F;n+[F;n+[F;n+[F;n+[F;n+F]]]]];
  %V = bsxfun(@plus,V,[0.5 0.5 0.5]);
  %[V,~,IM] = remove_duplicate_vertices(V,eps);
  %F = IM(F);

  %max_vol = mean(doublearea(V,F)).^2;
  %tsurf(F,V);
  %drawnow;
  %[V,T] = tetgen(V,F,'Flags',sprintf('-q10 -a%0.17f',max_vol));
end
