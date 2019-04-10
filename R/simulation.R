# n50 <- function(m){
#   s <- sum(m)
#   m <- m[order(-m)]
#   c_s <- cumsum(m)
#   i <- which(c_s > 0.5*s)[1]
#
#   return(m[i])
# }
#
# plot_mol_hist <- function(m_l){
#   plot_ly(x = m_l, type = "histogram")
# }
#
# subplot(plot_ly(x = seq(0, 0.5, 0.01), y = dgamma(seq(0, 0.5, 0.01), shape = 3.76608671, rate = 54.72492228)), plot_ly(x  = x$x, type = "histogram"))
