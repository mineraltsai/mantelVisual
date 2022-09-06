"corrposcalcu"<-function (corr,
                          method = c("circle", "square", "ellipse", "number",
                                     "shade", "color", "pie"),
                          type = c("full", "lower", "upper"),
                          col = NULL,
                          col.lim = NULL,
                          is.corr = TRUE,
                          add = FALSE, diag = TRUE,
                          p.mat = NULL)
{
  method = match.arg(method)
  type = match.arg(type)



  if (is.null(col.lim)) {
    if (is.corr) {
      col.lim = c(-1, 1)
    }
    else {
      if (!diag) {
        diag(corr) = NA
      }
      col.lim = c(min(corr, na.rm = TRUE), max(corr, na.rm = TRUE))
    }
  }
  SpecialCorr = 0
  if (is.corr) {
    if (min(corr, na.rm = TRUE) < -1 - .Machine$double.eps^0.75 ||
        max(corr, na.rm = TRUE) > 1 + .Machine$double.eps^0.75) {
      stop("The matrix is not in [-1, 1]!")
    }
    SpecialCorr = 1

  }
  intercept = 0
  zoom = 1
  if (!is.corr) {
    c_max = max(corr, na.rm = TRUE)
    c_min = min(corr, na.rm = TRUE)
    if ((col.lim[1] > c_min) | (col.lim[2] < c_max)) {
      stop("Wrong color: matrix should be in col.lim interval!")
    }
    if (diff(col.lim)/(c_max - c_min) > 2) {
      warning("col.lim interval too wide, please set a suitable value")
    }
    if (c_max <= 0 | c_min >= 0) {
      intercept = -col.lim[1]
      zoom = 1/(diff(col.lim))
      if (col.lim[1] * col.lim[2] < 0) {
        warning("col.lim interval not suitable to the matrix")
      }
    }
    else {
      stopifnot(c_max * c_min < 0)
      stopifnot(c_min < 0 && c_max > 0)
      intercept = 0
      zoom = 1/max(abs(col.lim))
      SpecialCorr = 1
    }
    corr = (intercept + corr) * zoom
  }
  col.lim2 = (intercept + col.lim) * zoom
  int = intercept * zoom

  n = nrow(corr)
  m = ncol(corr)
  min.nm = min(n, m)
  ord = 1:min.nm
  if (is.null(rownames(corr))) {
    rownames(corr) = 1:n
  }
  if (is.null(colnames(corr))) {
    colnames(corr) = 1:m
  }
  apply_mat_filter = function(mat) {
    x = matrix(1:n * m, nrow = n, ncol = m)
    switch(type, upper = mat[row(x) > col(x)] <- Inf, lower = mat[row(x) <
                                                                    col(x)] <- Inf)
    if (!diag) {
      diag(mat) = Inf
    }
    return(mat)
  }
  getPos.Dat = function(mat) {
    tmp = apply_mat_filter(mat)
    Dat = tmp[is.finite(tmp)]
    ind = which(is.finite(tmp), arr.ind = TRUE)
    Pos = ind
    Pos[, 1] = ind[, 2]
    Pos[, 2] = -ind[, 1] + 1 + n
    PosName = ind
    PosName[, 1] = colnames(mat)[ind[, 2]]
    PosName[, 2] = rownames(mat)[ind[, 1]]
    return(list(Pos, Dat, PosName))
  }


  Pos = getPos.Dat(corr)[[1]]
  PosName = getPos.Dat(corr)[[3]]
  DAT = getPos.Dat(corr)[[2]]
  pNew = getPos.Dat(p.mat)[[2]]

  corrPos = data.frame(PosName, Pos, DAT)
  colnames(corrPos) = c("xName", "yName", "x", "y", "corr")
  if (!is.null(p.mat)) {
    corrPos = cbind(corrPos, pNew)
    colnames(corrPos)[6] = c("p.value")
  }
  corrPos = corrPos[order(corrPos[, 3], -corrPos[, 4]), ]
  rownames(corrPos) = NULL
  res = list(corrPos = corrPos)
}
