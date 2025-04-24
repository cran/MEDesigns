#' Analysis of ME-designs in Block Set-up
#'
#' @param data Columns of dataset should be in order of block, line1,line2, cross number and response.
#'@description
#'For a given field data it will provide analysis result through ANOVA table including gca and sca effect analysis.
#'
#' @returns Returns the ANOVA table of gca and sca effect analysis.
#' @export
#'
#' @examples
#' library(MEDesigns)
#' MEBanalysis(MEdata)
MEBanalysis<-function(data){
  colnames(data)<-c("block", 'line1', 'line2', 'cross', 'yld')
  # Convert to factors
  data$block <- as.factor(data$block)
  data$cross <- as.factor(data$cross)
  
  # Fit ANOVA model
  model <- aov(yld ~ block + cross, data = data)
  anova_table <- summary(aov(model))
  print(anova_table)
  # Extract error MS for F-tests
  error_df <- anova_table[[1]]["Residuals", "Df"]
  error_ms <- anova_table[[1]]["Residuals", "Mean Sq"]
  
  # Prepare matrices for GCA and SCA analysis
  n_no <- max(data$line2)
  c_no <- max(as.numeric(data$cross))
  n_obs <- nrow(data)
  
  # Create indicator matrices
  sca <- model.matrix(~ cross - 1, data)
  block_mat <- model.matrix(~ block - 1, data)
  #col_mat <- model.matrix(~ as.factor(column) - 1, data)
  
  # Create X1 and X2 matrices
  X1 <- sca
  X2 <- cbind(1, block_mat)
  
  # Create C matrix
  C <- (t(X1) %*% X1) - (t(X1) %*% X2) %*% MASS::ginv(t(X2) %*% X2) %*% (t(X2) %*% X1)
  
  # Create Q matrix (relationship between lines and crosses)
  line_cross <- cbind(data$line1, data$line2)
  x11 <- matrix(0, nrow = n_obs, ncol = max(line_cross))
  
  for (i in 1:n_obs) {
    for (j in 1:2) {
      if (line_cross[i,j] > 0) {
        x11[i, line_cross[i,j]] <- x11[i, line_cross[i,j]] + 1
      }
    }
  }
  
  rep_mat <- t(sca) %*% sca
  q1 <- t(x11) %*% sca
  q <- q1 / rep_mat[1,1]
  qq <- q %*% t(q)
  inv_qq <- solve(qq)
  
  # Create H matrices
  jpv <- matrix(1, nrow = n_no, ncol = c_no)
  jpv1 <- jpv * (0.5 / c_no)
  h1 <- (inv_qq %*% q) - jpv1
  h2 <- diag(c_no) - (t(q) %*% inv_qq %*% (q))
  
  # Calculate C components
  c11 <- (h1) %*% C %*% t(h1)
  c22 <- (h2) %*% C %*% t(h2)
  
  # Calculate cross totals
  cross_t <- aggregate(yld ~ cross, data = data, sum)$yld
  block_t <- aggregate(yld ~ block, data = data, sum)$yld
  #col_t <- aggregate(yld ~ column, data = data, sum)$yld
  
  # Calculate adjusted cross totals
  n1 <- t(sca) %*% block_mat
  #n2 <- t(sca) %*% col_mat
  #m_mat <- t(row_mat) %*% col_mat
  k1 <- t(block_mat) %*% block_mat
  #k2 <- t(col_mat) %*% col_mat
  
  # adj_cross <- cross_t - (n1 %*% MASS::ginv(k1) %*% row_t) -
  #   (n2 - n1 %*% MASS::ginv(k1) %*% m_mat) %*%
  #   MASS::ginv(k2 - t(m_mat) %*% MASS::ginv(k1) %*% m_mat) %*%
  #   (col_t - t(m_mat) %*% MASS::ginv(k1) %*% row_t)
  
  adj_cross <- cross_t - (n1 %*% MASS::ginv(k1) %*% block_t)
  # Calculate SS for GCA and SCA
  df_gca <- n_no - 1
  df_sca <- c_no - n_no
  
  ss_gca <- t(adj_cross) %*% t(h1) %*% MASS::ginv(c11) %*% (h1) %*% adj_cross
  ss_sca <- t(adj_cross) %*% t(h2) %*% MASS::ginv(c22) %*% (h2) %*% adj_cross
  
  ms_gca <- ss_gca / df_gca
  ms_sca <- ss_sca / df_sca
  
  # Calculate F-values and p-values
  f_gca <- ms_gca / error_ms
  f_sca <- ms_sca / error_ms
  
  p_gca <- 1 - pf(f_gca, df_gca, error_df)
  p_sca <- 1 - pf(f_sca, df_sca, error_df)
  
  # Critical F-values
  f_crit_gca <- qf(0.975, df_gca, error_df)
  f_crit_sca <- qf(0.975, df_sca, error_df)
  
  # Format p-values
  format_pval <- function(p) {
    if (p < 0.0001) return("<.0001")
    return(sprintf("%.4f", p))
  }
  
  # Print results
  cat("GCA Analysis:\n")
  cat(sprintf("DF: %d, SS: %.2f, MS: %.2f, F-value: %.2f, Critical F: %.2f, p-value: %s\n",
              df_gca, ss_gca, ms_gca, f_gca, f_crit_gca, format_pval(p_gca)))
  
  cat("\nSCA Analysis:\n")
  cat(sprintf("DF: %d, SS: %.2f, MS: %.2f, F-value: %.2f, Critical F: %.2f, p-value: %s\n",
              df_sca, ss_sca, ms_sca, f_sca, f_crit_sca, format_pval(p_sca)))
}