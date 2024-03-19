

batchParam <- function(FUN, ...) {
  .mapply(
    FUN = FUN, 
    dots = as.list.data.frame(data.frame(...)), # use ?base::data.frame to recycle the length
    MoreArgs = NULL
  )
}



as_geomfunctionargs <- function(x) {
  # `x` is return of [batchParam]
  ans <- lapply(x, FUN = as.list)
  names(ans) <- vapply(x, FUN = function(i) {
    nm <- names(i)
    nm[nm == 'alpha'] <- '\u03b1'
    nm[nm == 'nu'] <- '\u03bd'
    nm[nm == 'xi'] <- '\u03be'
    nm[nm == 'omega'] <- '\u03c9'
    paste0(sprintf(fmt = '%s = %.3g', nm, i), collapse = '\n')
  }, FUN.VALUE = NA_character_)
  return(ans)
}








