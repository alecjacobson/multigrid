
iters = 1000000;
tol = 1e-7;
sam = 2;
while true
  [V,F] = create_regular_grid(sam+1,sam+1,0,0);
  if sam > 2
    Zs0 = prolong(CV,CF,CZs,V);
    Zr0 = prolong(CV,CF,CZr,V);
  else
    Zr0 = zeros(size(V,1),1);
    Zs0 = zeros(size(V,1),1);
  end
  b = unique(outline(F));
  bc = zeros(numel(b),1);
  w = 1e4;
  A = -cotmatrix(V,F) + sparse(b,b,w,size(V,1),size(V,1));
  M = massmatrix(V,F,'voronoi');
  B = M*ones(size(V,1),1) + sparse(b,1,w,size(V,1),1);


    R = B-A*Zr0;
    Yr0 = zeros(size(V,1),1);
    Yr = solve_pcg(A,R,tol/norm(R),iters,Yr0);
    Zr = Yr + Zr0;

  Zs = solve_pcg(A,B,tol/norm(B),iters,Zs0);
  %Zs = A\B;

  nsp = 4;
  subplot(1,nsp,1);
  tsurf(F,[V Zs],'EdgeAlpha',1/sam);
  subplot(1,nsp,2);
  tsurf(F,[V R],'EdgeAlpha',1/sam);
  subplot(1,nsp,3);
  tsurf(F,[V Yr],'EdgeAlpha',1/sam);
  subplot(1,nsp,4);
  tsurf(F,[V Zr],'EdgeAlpha',1/sam);
  % prepare for next iter
  CF = F;
  CV = V;
  CZs = Zs;
  CZr = Zr;
  sam = sam*2;
  input('');
end
