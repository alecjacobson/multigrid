function Z = solve_pcg(A,B,tol,iters,Z0)
  L = ichol(A);
  %L = [];
  Z = pcg(A,B,tol,iters,L,L',Z0);
end
