"cor_calcu"<-function (x, y = NULL, use = "pairwise.complete.obs", 
          alternative = c("two.sided", "less", "greater"), ...) 
{
  ia = match.arg(alternative)
  cor = cor(x, y, use = use, ...)
  x = as.matrix(x)
  finMat = !is.na(x)
  if (is.null(y)) {
    np = t(finMat) %*% finMat
  } else {
    y = as.matrix(y)
    np = t(finMat) %*% (!is.na(y))
  }
  Z = 0.5 * log((1 + cor)/(1 - cor)) * sqrt(np - 2)
  if (ia == "two.sided") {
    T = sqrt(np - 2) * abs(cor)/sqrt(1 - cor^2)
    p = 2 * pt(T, np - 2, lower.tail = FALSE)
  } else if (ia == "less") {
    T = sqrt(np - 2) * cor/sqrt(1 - cor^2)
    p = pt(T, np - 2, lower.tail = TRUE)
  } else if (ia == "greater") {
    T = sqrt(np - 2) * cor/sqrt(1 - cor^2)
    p = pt(T, np - 2, lower.tail = FALSE)
  }
  list(cor = cor, p = p, Z = Z, t = T, nObs = np)
}
