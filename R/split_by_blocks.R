split_by_blocks <- function(hh, prop = 0.8, iterations = 5000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.list(hh), length(hh) >= 1)
  
  block_sizes <- vapply(hh, length, integer(1))
  all_ids     <- sort(unique(unlist(hh)))
  n_total     <- length(all_ids)
  target      <- round(prop * n_total)
  
  best_sel   <- rep(FALSE, length(hh))
  best_sum   <- 0L
  best_diff  <- Inf
  
  for (it in seq_len(iterations)) {
    ord <- sample(seq_along(hh))
    sel <- rep(FALSE, length(hh))
    s   <- 0L
    
    # 贪心：按随机顺序遍历 block，若加入能让总数更接近 target 就加入
    for (j in ord) {
      new_s <- s + block_sizes[j]
      if (abs(new_s - target) <= abs(s - target)) {
        sel[j] <- TRUE
        s <- new_s
      }
      # 小优化：已经非常接近就可以提前结束
      if (abs(s - target) <= 1L) break
    }
    
    d <- abs(s - target)
    if (d < best_diff || (d == best_diff && s > best_sum)) {
      best_diff <- d
      best_sum  <- s
      best_sel  <- sel
      if (best_diff == 0L) {
        # 已经正好命中目标，直接退出
        break
      }
    }
  }
  
  train_ids <- sort(unique(unlist(hh[best_sel])))
  test_ids  <- setdiff(all_ids, train_ids)
  
  # 生成按样本编号排序的 split 向量
  split <- rep(NA_character_, max(all_ids))
  split[train_ids] <- "train"
  split[test_ids]  <- "test"
  names(split) <- seq_along(split)
  
  list(
    train_ids = train_ids,
    test_ids  = test_ids,
    split     = split,
    sizes     = c(train = length(train_ids), test = length(test_ids), total = n_total),
    target    = c(target_train = target, achieved_train = length(train_ids), diff = length(train_ids) - target),
    block_sizes = block_sizes
  )
}
