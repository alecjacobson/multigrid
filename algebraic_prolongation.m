function P = algebraic_prolongation(A,l)

  % "Algebraic multigrid by smoothed aggregation for second and fourth order
  % elliptic problems" [Vanek et al. 1995]
  
  % problem size
  n = size(A,1);
  % entries in A
  [I,J] = find(A);
  D = diag(A);
  % Geometric mean of diagonal: Gij = sqrt(Aii Ajj)
  G = sparse(I,J,sqrt(D(I).*D(J)),n,n);
  % strong connections
  epsilon = 0.08 * (1/2)^(l-1);
  % Neighborhood matrix
  % (can't use >= if we want to stay sparse)
  N = abs(A) > epsilon*G;

  % Step 1:
  R = 1:n;
  C = zeros(n,1);
  j = 0;
  for i = 1:n
    Ni = find(N(i,:));
    if numel(Ni)>1 && all(ismember(Ni,R))
      j = j+1;
      C(Ni) = j;
      R = setdiff(R,Ni);
    end
  end

  % step 2:
  for i = 1:max(C)
    Ci = find(C'==i);

    for j = Ci
      Nj = setdiff(find(N(j,:)),j);
      Pj = intersect(Nj,R);

      for k = Pj
        C(k) = i;
        R = setdiff(R,k);
      end
    end
  end

  % Disconnected nodes are all lumped into a single aggregate
  C(C==0) = max(C)+1;

  P = sparse(1:n,C',1,n,max(C));

  AF = A.*N;
  % diag(AF) == diag(A) since N has ones along the diagonal
  dAF = diag(A) - sum(A-AF,2);
  AF = AF + diag(dAF);

  % special Jacobi iteration
  D = diag(diag(AF));
  w = 2/3;
  P = P - w * D \(AF *P);

  %[NI,NJ] = find(N);
  %clf;
  %hold on;
  %tsurf(F,V);
  %text(V(:,1),V(:,2),V(:,1)*0,num2str(C),'BackgroundColor',[.7 .7 .7],'Fontsize',13);
  %plot_edges(V,[NI NJ],'g','LineWidth',4);
  %hold off;
  %axis equal;

end
