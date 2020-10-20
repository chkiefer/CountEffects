# Desired new syntax ------------------------------------------------

# model0 <- '
# outcome:
# y |~ poisson
#
# treatment:
# x |~ categorical
#
# covariates:
# gender |~ categorical
# eta + z1 |~ mvnorm
# z2 ~ eta + z1
# z2 |~ poisson
#
# ## Bla bla bla
# latent variables:
# eta =~ 1*ind1 + la2*ind2 + la3*ind3
# ind1 ~ 0*1
# eta |~ normal
# ind1 + ind2 + ind3 |~ skewnormal
#'

ParseModelString <- function(model.syntax){
  if (length(model.syntax) == 0) {
    stop("CountEffects ERROR: empty model syntax")
  }
  model.syntax <- gsub("[#!].*(?=\n)", "", model.syntax,
                       perl = TRUE)
  model.syntax <- gsub(";", "\n", model.syntax,
                       fixed = TRUE)
  model.syntax <- gsub("[ \t]+", "", model.syntax,
                       perl = TRUE)
  model.syntax <- gsub("\n{2,}", "\n", model.syntax,
                       perl = TRUE)
  model.syntax <- gsub(pattern = "˜", replacement = "~",
                       model.syntax)
  model <- unlist(strsplit(model.syntax, "\n"))
  model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL",
                       model)
  operators <- c("=~", "~~",
                 "~", "==", ":=",
                 ":",  "%", "\\|~")
  start.idx <- grep(paste(operators, collapse = "|"),
                    model.simple)
  if (length(start.idx) == 0L) {
    stop("CountEffects ERROR: model does not contain model syntax (no operators found)")
  }
  if (start.idx[1] > 1L) {
    for (el in 1:(start.idx[1] - 1L)) {
      if (nchar(model.simple[el]) > 0L) {
        warning("CountEffects WARNING: no operator found in this syntax line: ",
                model.simple[el], "\n", "                  This syntax line will be ignored!")
      }
    }
  }
  end.idx <- c(start.idx[-1] - 1, length(model))
  model.orig <- model
  model <- character(length(start.idx))
  for (i in 1:length(start.idx)) {
    model[i] <- paste(model.orig[start.idx[i]:end.idx[i]],
                      collapse = "")
  }
  model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL",
                       model)
  idx.wrong <- which(!grepl(paste(operators, collapse = "|"),
                            model.simple))
  if (length(idx.wrong) > 0) {
    cat("CountEffects: missing operator in formula(s):\n")
    print(model[idx.wrong])
    stop("CountEffects ERROR: syntax error in model syntax")
  }
  idx.wrong <- which(grepl("^\\+", model))
  if (length(idx.wrong) > 0) {
    cat("CountEffects: some formula(s) start with a plus (+) sign:\n")
    print(model[idx.wrong])
    stop("CountEffects ERROR: syntax error in lavaan model syntax")
  }

  FLAT.lhs <- character(0)
  FLAT.op <- character(0)
  FLAT.rhs <- character(0)
  FLAT.rhs.mod.idx <- integer(0)
  FLAT.block <- integer(0)
  FLAT.fixed <- character(0)
  FLAT.start <- character(0)
  FLAT.lower <- character(0)
  FLAT.upper <- character(0)
  FLAT.label <- character(0)
  FLAT.prior <- character(0)
  FLAT.efa <- character(0)
  FLAT.idx <- 0L
  MOD.idx <- 0L
  CON.idx <- 0L
  MOD <- vector("list", length = 0L)
  CON <- vector("list", length = 0L)
  BLOCK <- character(0)
  BLOCK_OP <- FALSE
  for (i in 1:length(model)) {
    x <- model[i]
    if (debug) {
      cat("formula to parse:\n")
      print(x)
      cat("\n")
    }
    line.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL",
                        x)
    if (grepl("=~", line.simple, fixed = TRUE)) {
      op <- "=~"
    } else if (grepl("~~", line.simple, fixed = TRUE)) {
      op <- "~~"
    } else if (grepl("|~", line.simple, fixed = TRUE)) {
      op <- "\\|~"
    } else if (grepl("~", line.simple, fixed = TRUE)) {
      op <- "~"
    } else if (grepl("==", line.simple, fixed = TRUE)) {
      op <- "=="
    } else if (grepl(":=", line.simple, fixed = TRUE)) {
      op <- ":="
    } else if (grepl(":", line.simple, fixed = TRUE)) {
      op <- ":"
    } else if (grepl("%", line.simple, fixed = TRUE)) {
      op <- "%"
    } else {
      stop("unknown operator in ", model[i])
    }
    if (substr(x, 1, 6) == "label(") stop("label modifier can not be used on the left-hand side of the operator")
    op.idx <- regexpr(op, x)

    lhs <- substr(x, 1L, op.idx - 1L)
    rhs <- substr(x, op.idx + attr(op.idx, "match.length"),
                  nchar(x))
    if (substr(rhs, 1, 1) == "+") {
      rhs <- substr(rhs, 2, nchar(rhs))
    }
    if (op == "==" || op == ":=") {
      lhs <- gsub("\\\"", "", lhs)
      rhs <- gsub("\\\"", "", rhs)
      CON.idx <- CON.idx + 1L
      CON[[CON.idx]] <- list(op = op, lhs = lhs, rhs = rhs,
                             user = 1L)
      next
    }

    if (op == "\\|~" ) {
      lhs <- gsub("\\\"", "", lhs)
      rhs <- gsub("\\\"", "", rhs)
      CON.idx <- CON.idx + 1L
      CON[[CON.idx]] <- list(op = op, lhs = lhs, rhs = rhs,
                             user = 1L)
      next
    }
    if (op == ":") {
      BLOCK_OP <- TRUE
      if (nchar(rhs) != 0L) {
        stop("CountEffects ERROR: Opposed to lavaan multilevel syntax, no syntax is allowed after :")
      }
      lhs.orig <- lhs
      lhs <- tolower(lhs)
      if (!lhs %in% c("outcome", "treatment", "covariates",
                      "latentvariables")) {
        stop("CountEffects ERROR: unknown block identifier")
      }
      FLAT.idx <- FLAT.idx + 1L
      FLAT.lhs[FLAT.idx] <- lhs
      FLAT.op[FLAT.idx] <- op
      FLAT.rhs[FLAT.idx] <- rhs
      FLAT.fixed[FLAT.idx] <- ""
      FLAT.start[FLAT.idx] <- ""
      FLAT.lower[FLAT.idx] <- ""
      FLAT.upper[FLAT.idx] <- ""
      FLAT.label[FLAT.idx] <- ""
      FLAT.prior[FLAT.idx] <- ""
      FLAT.efa[FLAT.idx] <- ""
      FLAT.rhs.mod.idx[FLAT.idx] <- 0L

      BLOCK <- lhs
      FLAT.block[FLAT.idx] <- BLOCK
      next
    }

    LHS <- strsplit(lhs, split = "+", fixed = TRUE)[[1]]
    LHS <- gsub("^\\S*\\*", "", LHS)
    if (!all(make.names(LHS) == LHS)) {
      stop("CountEffects ERROR: left hand side (lhs) of this formula:\n    ",
           lhs, " ", op, " ", rhs, "\n    contains either a reserved word (in R) or an illegal character: ",
           dQuote(LHS[!make.names(LHS) == LHS]), "\n    See ?reserved for a list of reserved words in R",
           "\n    Please use a variable name that is not a reserved word in R",
           "\n    and use only characters, digits, or the dot symbol.")
    }
    lhs.formula <- as.formula(paste("~", lhs))
    lhs.out <- syntax_parse_rhs(rhs = lhs.formula[[2L]],
                                    op = op)
    lhs.names <- names(lhs.out)
    rhs <- gsub("\\(?([-]?[0-9]*\\.?[0-9]*)\\)?\\?",
                "start(\\1)\\*", rhs)
    RHS <- strsplit(rhs, split = "+", fixed = TRUE)[[1]]
    RHS.names <- gsub("^\\S*\\*", "", RHS)
    BAD <- c("if", "else", "repeat", "while",
             "function", "for", "in")
    if (any(RHS.names %in% BAD)) {
      stop("lavaan ERROR: right hand side (rhs) of this formula:\n    ",
           lhs, " ", op, " ", rhs, "\n    contains either a reserved word (in R) or an illegal character: ",
           dQuote(RHS.names[which(RHS.names %in% BAD)[1]]),
           "\n    See ?reserved for a list of reserved words in R",
           "\n    Please use a variable name that is not a reserved word in R",
           "\n    and use only characters, digits, or the dot symbol.")
    }
    RHS <- strsplit(rhs, split = "+", fixed = TRUE)[[1]]
    RHS.labels <- gsub("\\*\\S*$", "", RHS)
    if (any(RHS.labels %in% BAD)) {
      stop("lavaan ERROR: right hand side (rhs) of this formula:\n    ",
           lhs, " ", op, " ", rhs, "\n    contains either a reserved word (in R) or an illegal character: ",
           dQuote(RHS.names[which(RHS.labels %in% BAD)[1]]),
           "\n    See ?reserved for a list of reserved words in R",
           "\n    Please use a variable name that is not a reserved word in R",
           "\n    and use only characters, digits, or the dot symbol.")
    }
    rhs.formula <- as.formula(paste("~", rhs))
    out <- syntax_parse_rhs(rhs = rhs.formula[[2L]],
                                op = op)
    if (debug)
      print(out)
    for (l in 1:length(lhs.names)) {
      for (j in 1:length(out)) {
        if (names(out)[j] == "intercept") {
          if (op == "~") {
            rhs.name <- ""
          }
          else {
            stop("lavaan ERROR: right-hand side of formula contains an invalid variable name:\n    ",
                 x)
          }
        }
        else if (names(out)[j] == "..zero.." &&
                 op == "~") {
          rhs.name <- ""
        }
        else if (names(out)[j] == "..constant.." &&
                 op == "~") {
          rhs.name <- ""
        }
        else {
          rhs.name <- names(out)[j]
        }
        if (op == "=~" && lhs.names[l] == names(out)[j]) {
          stop("lavaan ERROR: latent variable `",
               lhs.names[l], "' can not be measured by itself")
        }
        if (op != "~~") {
          idx <- which(FLAT.lhs == lhs.names[l] & FLAT.op ==
                         op & FLAT.block == BLOCK & FLAT.rhs == rhs.name)
          if (length(idx) > 0L) {
            stop("lavaan ERROR: duplicate model element in: ",
                 model[i])
          }
        }
        else {
          idx <- which(FLAT.lhs == rhs.name & FLAT.op ==
                         "~~" & FLAT.block == BLOCK & FLAT.rhs ==
                         lhs.names[l])
          if (length(idx) > 0L) {
            stop("lavaan ERROR: duplicate model element in: ",
                 model[i])
          }
        }
        FLAT.idx <- FLAT.idx + 1L
        FLAT.lhs[FLAT.idx] <- lhs.names[l]
        FLAT.op[FLAT.idx] <- op
        FLAT.rhs[FLAT.idx] <- rhs.name
        FLAT.block[FLAT.idx] <- BLOCK
        FLAT.fixed[FLAT.idx] <- ""
        FLAT.start[FLAT.idx] <- ""
        FLAT.label[FLAT.idx] <- ""
        FLAT.lower[FLAT.idx] <- ""
        FLAT.upper[FLAT.idx] <- ""
        FLAT.prior[FLAT.idx] <- ""
        FLAT.efa[FLAT.idx] <- ""
        mod <- list()
        rhs.mod <- 0L
        if (length(lhs.out[[l]]$efa) > 0L) {
          mod$efa <- lhs.out[[l]]$efa
          FLAT.efa[FLAT.idx] <- paste(mod$efa, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$fixed) > 0L) {
          mod$fixed <- out[[j]]$fixed
          FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$start) > 0L) {
          mod$start <- out[[j]]$start
          FLAT.start[FLAT.idx] <- paste(mod$start, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$lower) > 0L) {
          mod$lower <- out[[j]]$lower
          FLAT.lower[FLAT.idx] <- paste(mod$lower, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$upper) > 0L) {
          mod$upper <- out[[j]]$upper
          FLAT.upper[FLAT.idx] <- paste(mod$upper, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$label) > 0L) {
          mod$label <- out[[j]]$label
          FLAT.label[FLAT.idx] <- paste(mod$label, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$prior) > 0L) {
          mod$prior <- out[[j]]$prior
          FLAT.prior[FLAT.idx] <- paste(mod$prior, collapse = ";")
          rhs.mod <- 1L
        }
        if (op == "=~" && rhs == "0") {
          mod$fixed <- 0
          FLAT.rhs[FLAT.idx] <- FLAT.lhs[FLAT.idx]
          FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse = ";")
          rhs.mod <- 1L
        }
        FLAT.rhs.mod.idx[FLAT.idx] <- rhs.mod
        if (rhs.mod > 0L) {
          MOD.idx <- MOD.idx + 1L
          MOD[[MOD.idx]] <- mod
        }
      }
    }
  }
}





syntax_parse_rhs <- function (rhs, op = "")
{
  out <- list()
  repeat {
    if (length(rhs) == 1L) {
      out <- c(vector("list", 1L), out)
      NAME <- all.vars(rhs)
      if (length(NAME) > 0L) {
        names(out)[1L] <- NAME
      }
      else {
        if (as.character(rhs) == "1") {
          names(out)[1L] <- "intercept"
        }
        else if (as.character(rhs) == "0") {
          names(out)[1L] <- "..zero.."
          out[[1L]]$fixed <- 0
        }
        else {
          names(out)[1L] <- "..constant.."
          out[[1L]]$fixed <- 0
        }
      }
      break
    }
    else if (rhs[[1L]] == "*") {
      out <- c(vector("list", 1L), out)
      NAME <- all.vars(rhs[[3L]])
      if (length(NAME) > 0L) {
        rhs3.names <- all.names(rhs[[3L]])
        if (rhs3.names[1L] == ":") {
          NAME <- paste(NAME[1L], ":", NAME[2L],
                        sep = "")
        }
        names(out)[1L] <- NAME
      }
      else {
        names(out)[1L] <- "intercept"
      }
      i.var <- all.vars(rhs[[2L]], unique = FALSE)
      if (length(i.var) > 0L) {
        out[[1L]]$label <- i.var
      }
      else {
        out[[1L]] <- syntax_get_modifier(rhs[[2L]])
      }
      break
    }
    else if (rhs[[1L]] == ":") {
      out <- c(vector("list", 1L), out)
      NAME <- all.vars(rhs)
      NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
      names(out)[1L] <- NAME
      break
    }
    else if (rhs[[1L]] == "+") {
      out <- c(vector("list", 1L), out)
      if (length(rhs[[3L]]) == 3L && rhs[[3L]][[1]] ==
          "*") {
        NAME <- all.vars(rhs[[3L]][[3]])
        if (length(NAME) > 0L) {
          rhs3.names <- all.names(rhs[[3L]][[3]])
          if (rhs3.names[1L] == ":") {
            NAME <- paste(NAME[1L], ":", NAME[2L],
                          sep = "")
          }
          names(out)[1L] <- NAME
        }
        else {
          names(out)[1L] <- "intercept"
        }
        i.var <- all.vars(rhs[[3]][[2L]], unique = FALSE)
        if (length(i.var) > 0L) {
          out[[1L]]$label <- i.var
        }
        else {
          out[[1L]] <- syntax_get_modifier(rhs[[3]][[2L]])
        }
      }
      else if (length(rhs[[3L]]) == 3L && rhs[[3L]][[1]] ==
               ":") {
        NAME <- all.vars(rhs[[3L]])
        NAME <- paste(NAME[1L], ":", NAME[2L],
                      sep = "")
        names(out)[1L] <- NAME
      }
      else {
        NAME <- all.vars(rhs[[3]])
        if (length(NAME) > 0L) {
          names(out)[1L] <- NAME
        }
        else {
          if (as.character(rhs[[3]]) == "1") {
            names(out)[1L] <- "intercept"
          }
          else if (as.character(rhs[[3]]) == "0") {
            names(out)[1L] <- "..zero.."
            out[[1L]]$fixed <- 0
          }
          else {
            names(out)[1L] <- "..constant.."
            out[[1L]]$fixed <- 0
          }
        }
      }
      rhs <- rhs[[2L]]
    }
    else {
      stop("lavaan ERROR: I'm confused parsing this line: ",
           rhs, "\n")
    }
  }
  if (length(out) > 1L) {
    rhs.names <- names(out)
    while (!is.na(idx <- which(duplicated(rhs.names))[1L])) {
      dup.name <- rhs.names[idx]
      orig.idx <- match(dup.name, rhs.names)
      merged <- c(out[[orig.idx]], out[[idx]])
      if (!is.null(merged))
        out[[orig.idx]] <- merged
      out <- out[-idx]
      rhs.names <- names(out)
    }
  }
  out
}


syntax_get_modifier <- function (mod)
{
  if (length(mod) == 1L) {
    if (is.numeric(mod))
      return(list(fixed = mod))
    if (is.na(mod))
      return(list(fixed = as.numeric(NA)))
    if (is.character(mod))
      return(list(label = mod))
  }
  else if (mod[[1L]] == "start") {
    cof <- unlist(lapply(as.list(mod)[-1], eval, envir = NULL,
                         enclos = NULL))
    return(list(start = cof))
  }
  else if (mod[[1L]] == "lower") {
    cof <- unlist(lapply(as.list(mod)[-1], eval, envir = NULL,
                         enclos = NULL))
    return(list(lower = cof))
  }
  else if (mod[[1L]] == "upper") {
    cof <- unlist(lapply(as.list(mod)[-1], eval, envir = NULL,
                         enclos = NULL))
    return(list(upper = cof))
  }
  else if (mod[[1L]] == "equal") {
    label <- unlist(lapply(as.list(mod)[-1], eval, envir = NULL,
                           enclos = NULL))
    return(list(label = label))
  }
  else if (mod[[1L]] == "label") {
    label <- unlist(lapply(as.list(mod)[-1], eval, envir = NULL,
                           enclos = NULL))
    label[is.na(label)] <- ""
    return(list(label = label))
  }
  else if (mod[[1L]] == "prior") {
    prior <- unlist(lapply(as.list(mod)[-1], eval, envir = NULL,
                           enclos = NULL))
    return(list(prior = prior))
  }
  else if (mod[[1L]] == "efa") {
    efa <- unlist(lapply(as.list(mod)[-1], eval, envir = NULL,
                         enclos = NULL))
    return(list(efa = efa))
  }
  else if (mod[[1L]] == "c") {
    cof <- unlist(lapply(as.list(mod)[-1], eval, envir = NULL,
                         enclos = NULL))
    if (all(is.na(cof))) {
      return(list(fixed = rep(as.numeric(NA), length(cof))))
    }
    else if (is.numeric(cof))
      return(list(fixed = cof))
    else if (is.character(cof)) {
      cof[is.na(cof)] <- ""
      return(list(label = cof))
    }
    else {
      stop("lavaan ERROR: can not parse modifier:",
           mod, "\n")
    }
  }
  else {
    cof <- try(eval(mod, envir = NULL, enclos = NULL), silent = TRUE)
    if (inherits(cof, "try-error")) {
      stop("lavaan ERROR: evaluating modifier failed: ",
           paste(as.character(mod)[[1]], "()", sep = ""),
           "\n")
    }
    else if (is.numeric(cof)) {
      return(list(fixed = cof))
    }
    else if (is.character(cof)) {
      return(list(label = cof))
    }
    else {
      stop("lavaan ERROR: can not parse modifier: ",
           paste(as.character(mod)[[1]], "()", sep = ""),
           "\n")
    }
  }
}


