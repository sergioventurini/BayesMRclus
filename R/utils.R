expit <- function(x) {
	return(1/(1 + exp(-x)))
}

list2matrix <- function(D) {
	return(t(sapply(D, as.numeric)))
}

list2array <- function(D) {
	S <- length(D)
	n <- nrow(as.matrix(D[[1]]))
	p <- ncol(as.matrix(D[[1]]))
	
	out <- array(NA, dim = c(n, p, S))
	for (s in 1:S) {
		out[, , s] <- as.matrix(D[[s]])
	}
	
	return(out)
}

list.sum <- function(x) {
	m <- length(x)
	z <- x[[1]]
	if (m == 1) {
		return(z)
	}
	for (j in 2:m) {
		z <- z + x[[j]]
	}
	
  return(z)
}

#' Auxiliary function to recursively check NAs in a list.
#'
#' \code{check_list_na()} compares two lists and fills in the missing
#'   elements in the first with those included in the second. The
#'   comparison is recursive in the sense that the process is repeated for
#'   all lists included in those given.
#'
#' @param orig A list whose content must be checked.
#' @param des A list to use as a reference with which compare the first one.
#'
#' @return A list with all elements added.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @examples
#' # [nothing for now]
#'
#' @export
check_list_na <- function(orig, des) {
  check_it <- function(o, d) {
    d.nm <- names(d)
    d.na <- is.na(match(d.nm, names(o)))
    if (any(d.na))
      o <- c(o, d[d.nm[which(d.na)]])

    return(o)
  }

  if (!is.list(orig))
    stop("the 'orig' argument must be a list")
  if (!is.list(des))
    stop("the 'des' argument must be a list")

  orig_new <- check_it(orig_new <- orig, des)

  for (el in 1:length(orig_new)) {
    if (is.list(orig_new[[el]]))
      orig_new[[el]] <- check_list_na(orig_new[[el]], des[[el]])
  }

  return(orig_new)
}

colMedians <- function(x, na.rm = TRUE, dims = 1L) {
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.array(x) || length(dn <- dim(x)) < 2L) 
    stop("'x' must be an array of at least two dimensions")
  if (dims < 1L || dims > length(dn) - 1L) 
    stop("invalid 'dims'")
  n <- prod(dn[id <- seq_len(dims)])
  dn <- dn[-id]
  z <- apply(x, 2, median, nar.rm = na.rm)
  if (length(dn) > 1L) {
    dim(z) <- dn
    dimnames(z) <- dimnames(x)[-id]
  }
  else names(z) <- dimnames(x)[[dims + 1L]]
  return(z)
}

bayesmr_pb <- function(min = 0, max = 1, initial = 0, char = "=", width = 49, skip = 5) {
  e <- new.env(parent = emptyenv())
  e$.val <- initial
  e$.killed <- FALSE
  e$.nb <- 0L
  nw <- nchar(char, "w")
  if (is.na(width)) {
    width <- getOption("width")
    width <- width - 10L
    width <- trunc(width/nw)
  }
  if (max <= min)
    stop("must have 'max' > 'min'")
  empty_string <- strrep(" ", skip)
  message(empty_string, "0%   10   20   30   40   50   60   70   80   90   100%")
  message(empty_string, "[----|----|----|----|----|----|----|----|----|----]")
  utils::flush.console()
  up3 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    e$.val <- value
    nb <- round(width * (value - min)/(max - min))
    if (nb == e$.nb) 
      return()
    message(paste0("\r     |", strrep(" ", nw * width + 6)), appendLF = FALSE)
    message(paste(c("\r     |", rep.int(char, nb), rep.int(" ", nw * (width - nb)), "|"), collapse = ""),
      appendLF = FALSE)
    utils::flush.console()
    e$.nb <- nb
  }
  getVal <- function() e$.val
  kill <- function() if (!e$.killed) {
    message(" ")
    utils::flush.console()
    e$.killed <- TRUE
  }
  up3(initial)
  structure(list(getVal = getVal, up = up3, kill = kill), class = "txtProgressBar")
}

bayesmr_setpb <- function(pb, value) {
    oldval <- pb$getVal()
    pb$up(value)
    invisible(oldval)
}

#' Check for suggested package (requireNamespace) and throw error if necessary
#'
#' @noRd
#' @param pkg Package name as a string.
#' @param min_version Optionally, a minimum version number as a string.
#' @return TRUE, invisibly, if no error is thrown.
#'
suggested_package <- function(pkg, min_version = NULL) {
  stopifnot(length(pkg) == 1, is.character(pkg))
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Please install the ",
      pkg, " package to use this function.",
      call. = FALSE
    )
  }

  if (!is.null(min_version)) {
    stopifnot(is.character(min_version))
    if (utils::packageVersion(pkg) < package_version(min_version)) {
      stop(
        "Version >=", min_version, " of the ",
        pkg, " package is required to use this function.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Explicit and/or regex parameter selection
#'
#' @noRd
#' @param explicit Character vector of selected parameter names.
#' @param patterns Character vector of regular expressions.
#' @param complete Character vector of all possible parameter names.
#' @return Characeter vector of combined explicit and matched (via regex)
#'   parameter names, unless an error is thrown.
#'
select_pars <- function(explicit = character(), patterns = character(), complete = character()) {
  stopifnot(is.character(explicit),
            is.character(patterns),
            is.character(complete))

  if (!length(explicit) && !length(patterns))
    return(complete)

  if (length(explicit)) {
    if (!all(explicit %in% complete)) {
      not_found <- which(!explicit %in% complete)
      stop(
        "Some 'pars' don't match parameter names: ",
        paste(explicit[not_found], collapse = ", ")
      )
    }
  }

  if (!length(patterns)) {
    return(unique(explicit))
  } else {
    regex_pars <-
      unlist(lapply(seq_along(patterns), function(j) {
        grep(patterns[j], complete, value = TRUE)
      }))
    if (!length(regex_pars))
      stop("no matches for 'regex_pars'.", call. = FALSE)
  }

  unique(c(explicit, regex_pars))
}

choose_colors <- function(n) {
  all_clrs <- unlist(bayesplot::color_scheme_get())
  clrs <- switch(
    as.character(n),
    "1" = get_color("m"),
    "2" = get_color(c("l", "d")),
    "3" = get_color(c("l", "m", "d")),
    "4" = all_clrs[-c(2, 4)],
    "5" = all_clrs[-3],
    "6" = all_clrs,
    rep_len(all_clrs, n)
  )
  unname(rev(clrs))
}

# Access a subset of the scheme colors
#
# @param level A character vector of level names (see scheme_level_names()). The
#   abbreviations "l", "lh", "m", "mh", "d", and "dh" can also be used instead
#   of the full names.
# @return A character vector of color values.
#
# [Source: bayesplot]
get_color <- function(levels) {
  sel <- which(!levels %in% scheme_level_names())
  if (length(sel)) {
    levels[sel] <- sapply(levels[sel], full_level_name)
  }
  stopifnot(all(levels %in% scheme_level_names()))
  color_vals <- bayesplot::color_scheme_get()[levels]
  unlist(color_vals, use.names = FALSE)
}

full_level_name <- function(x) {
  switch(x,
         l = "light", lh = "light_highlight",
         m = "mid", mh = "mid_highlight",
         d = "dark", dh = "dark_highlight")
}

# Color scheme level names
#
# [Source: bayesplot]
scheme_level_names <- function() {
  c("light",
    "light_highlight",
    "mid",
    "mid_highlight",
    "dark",
    "dark_highlight")
}

# Print a matrix in a pretty way
print_matrix <- function(mat, rownm = NULL, colnm = NULL, colwidth = 10, between_cols = 2, ndigits = 2, shift = 0,
  isint = FALSE) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  if (is.null(rownm)) {
    rownm <- paste0("Row ", 1:nr)
  }
  if (is.null(colnm)) {
    colnm <- paste0("Col ", 1:nc)
  }
  stopifnot(length(rownm) == nr, length(colnm) == nc)
  maxwidth <- max(nchar(format(round(mat, digits = ndigits), nsmall = ifelse(isint, 0, 2))))
  if (colwidth < maxwidth) colwidth <- maxwidth
  mat_str <- format(round(mat, digits = ndigits), nsmall = ifelse(isint, 0, 2), width = colwidth)
  if (any(is.na(mat))) mat_str <- gsub("NA", " -", mat_str)

  # if (any(nchar(rownm) > colwidth)) {
  #   rownm <- abbreviate(rownm, minlength = colwidth, strict = TRUE, named = FALSE)
  # }
  if (any(nchar(colnm) > colwidth)) {
    colnm <- abbreviate(colnm, minlength = colwidth, strict = TRUE, named = FALSE)
  }
  firstcolwidth <- max(nchar(rownm))
  empty_string_shift <- strrep(" ", shift)
  empty_string_firstcol <- strrep(" ", firstcolwidth)
  empty_string_between_cols <- strrep(" ", between_cols)
  empty_string_rows <- empty_string_cols <- character(nc)

  cat(empty_string_shift, sep = "")
  cat(empty_string_firstcol, sep = "")
  cat(empty_string_between_cols, sep = "")
  for (j in 1:nc) {
    empty_string_cols[j] <- strrep(" ", colwidth - nchar(colnm[j]))
    cat(empty_string_cols[j], colnm[j], sep = "")
    cat(empty_string_between_cols, sep = "")
  }
  cat("\n")
  for (i in 1:nr) {
    empty_string_rows[i] <- strrep(" ", firstcolwidth - nchar(rownm[i]))
    cat(empty_string_shift, sep = "")
    cat(rownm[i], empty_string_rows[i], sep = "")
    cat(empty_string_between_cols, sep = "")
    for (j in 1:nc) {
      cat(mat_str[i, j], sep = "")
      cat(empty_string_between_cols, sep = "")
      if (j == nc) cat("\n")
    }
  }
}

#' Auxiliary function to compute the I-squared heterogeneity index.
#'
#' \code{check_list_na()} compares two lists and fills in the missing
#'   elements in the first with those included in the second. The
#'   comparison is recursive in the sense that the process is repeated for
#'   all lists included in those given.
#'
#' @param y A list whose content must be checked.
#' @param s A list to use as a reference with which compare the first one.
#'
#' @return A length-one numeric vector.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @examples
#' # [nothing for now]
#'
#' @export
I2 <- function(y, s) {
  k <- length(y)
  w <- 1/s^2
  sum.w <- sum(w)
  mu.hat <- sum(y * w)/sum.w
  Q <- sum(w * (y - mu.hat)^2)
  Isq <- (Q - (k - 1))/Q
  Isq <- max(0, Isq)
  return(Isq)
}
