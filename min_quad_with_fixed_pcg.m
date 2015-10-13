function Z = min_quad_with_fixed_pcg(A,B,known,Y,tol,iter,Z0)
  n = size(A,1);
  in = setdiff(1:n,known)';
  Ain = A(in,in);
  L = ichol(Ain);
  Z = zeros(n,1);
  Z(known) = Y;
  if ~isempty(Z0)
    Z0 = Z0(in,:);
  end
  Z(in) = pcg(Ain,-(A(in,known)*Y)-0.5*B(in),tol,iter,L,L',Z0);
end
