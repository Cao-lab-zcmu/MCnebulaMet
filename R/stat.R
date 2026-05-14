run_univariate <- function(mat, group, method = c("t.test", "wilcox.test"),
                          p_adjust_method = "BH") {
  
  method <- match.arg(method)
  
  g1 <- unique(group)[1]
  g2 <- unique(group)[2]
  
  idx1 <- which(group == g1)
  idx2 <- which(group == g2)
  
  res_list <- apply(mat, 1, function(x) {
    
    x1 <- x[idx1]
    x2 <- x[idx2]
    
    # fold change（log2 scale）
    fc <- mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)
    
    # 统计检验
    pval <- tryCatch({
      if (method == "t.test") {
        t.test(x1, x2)$p.value
      } else {
        wilcox.test(x1, x2)$p.value
      }
    }, error = function(e) NA)
    
    c(log2FC = fc, pvalue = pval)
  })
  
  res <- as.data.frame(t(res_list))
  
  res$padj <- p.adjust(res$pvalue, method = p_adjust_method)
  res$.features_id <- rownames(res)
  
  res <- dplyr::as_tibble(res) |>
    dplyr::relocate(.features_id)
  
  return(res)
}
