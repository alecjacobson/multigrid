% This demonstrates that Galerkin MG is practical on regular grids and on
% irregular grids in 2D.


nruns = 10;

S = (2.^(7:-1:2))+1;
nnzr =    zeros(numel(S),nruns);
nnzpr =   zeros(numel(S),nruns);
nnzr_P =  zeros(numel(S),nruns);
nnzpr_P = zeros(numel(S),nruns);

%create = @(s) create_regular_grid(s,s);
%create = @(s) create_irregular_grid(s,s);
%create = @(s) regular_tetrahedral_mesh(s);
create = @(s) irregular_tetrahedral_mesh(s);

for run = 1:nruns
%clf;
%hold on;
  s = S(1);
  fprintf('create\n');
  [V,Ele] = create(s);
  fprintf('cotmatrix\n');
  L1 = cotmatrix(V,Ele);
  % running operator
  L = L1;
  
  for si = 1:numel(S)
    fprintf('%d/%d of %d/%d ...\n',si,numel(S),run,nruns);
    V_prev = V;
    Ele_prev = Ele;
    s = S(si);
  
    if si ~= 1
      [V,Ele] = create(s);
    end
    %tsurf(Ele,bsxfun(@plus,V,[si*1.05 0]),'LineWidth',0.5);
  
    P = prolongation(V,Ele,V_prev);
    Li = cotmatrix(V,Ele);
    L = P'*L*P;
    nnzr   (si,run) = nnz(Li)/numel(Li);
    nnzpr  (si,run) = nnz(Li)/size(Li,1);
    nnzr_P (si,run) = nnz(L)/numel(L);
    nnzpr_P(si,run) = nnz(L)/size(L,1);
  end
%hold off
%axis equal;
%pause

  clf;
  hold on;
  %plot(nnzpr(:),'-','LineWidth',4);
  %plot(nnzpr_P(:),'--','LineWidth',5);
  plot(  nnzpr(:,1:run),'-o','LineWidth',1,'MarkerSize',13,'Color','0 0.447 0.741');
  plot(nnzpr_P(:,1:run),'-o','LineWidth',1,'MarkerSize',10,'Color','0.85 0.325 0.098');
  hold off;
  title('# non-zeros per row','FontSize',15);
end
