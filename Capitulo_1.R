TCS <- function(x){
  (100/x)*4+(x/2)*8
}

TOC <- function()

resposta <- NULL

for (x in seq(1:15)) {
  i <- x
  resposta1 <- TCS(x)
  resposta <- rbind(resposta, data.frame(i,resposta1))
}




min(resposta$resposta1)


