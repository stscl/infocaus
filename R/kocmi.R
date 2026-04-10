kocmi = \(data, target, agent, conds, knockoff, null_knockoff = NULL, 
          type = c("cont", "disc"), nboots = 1e4, k = 3, threads = 1, 
          seed = 123456789, base = exp(1), method = "equal", contain_null = TRUE) {
  type = match.arg(type)
  mat = .convert2mat(data, type)
  return(RcppKOCMI(mat, target, agent, conds, knockoff, null_knockoff, type, nboots, 
                   k, 0, threads, seed, base, method, contain_null))
}
