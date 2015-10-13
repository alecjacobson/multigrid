function Y = occurs_once(X)
  sX = sort(X,2);
  [u,m,n] = unique(sX,'rows');
  c = accumarray(n(:),1);
  sX1 = u(c==1,:);
  sJ1 = m(c==1,:);
end
