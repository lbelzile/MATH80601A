
#' Update for Markov chain Monte Carlo
#' 
#' This function computes a single update for a Metropolis-Hastings 
#' algorithm for univariate proposal based on random walk, adjusted Langevin
#' dynamics with gradient or a quadratic Laplace approximation.
#' 
#' @inheritParams univariate_laplace_update
#' @return a list with components
#' \itemize{
#' \item{\code{ll}}: value of the log-likelihood, with potential additional parameters passed as attributes;
#' \item{\code{cur}}: values of the parameters after update;
#' \item{\code{accept}}: logical indicating the proposal has been accepted (\code{TRUE}) or rejected (\code{FALSE}).
#' }
#' @export
mcmc_update <- function(
    type = c("rw","mala","laplace"),
    par_curr,
    par_name,
    loglik,
    loglik_grad = NULL,
    loglik_hessian = NULL,
    logprior,
    logprior_grad = NULL,
    logprior_hessian = NULL,
    lb = -Inf,
    ub = Inf,
    transform = FALSE,
    damping = 1,
    sd_prop,
    ...){
  type <- match.arg(type)
  if(type == "rw"){
    univariate_rw_update(par_curr = par_curr, par_name = par_name, loglik = loglik, logprior = logprior,
                         sd_prop = sd_prop, lb = lb, ub = ub, transform = transform,...)
  } else if(type == "mala"){
    univariate_rw_update(par_curr = par_curr, par_name = par_name, loglik = loglik, logprior = logprior,
                         loglik_grad = loglik_grad, logprior_grad = logprior_grad, damping = damping, transform = transform,
                         sd_prop = sd_prop, lb = lb, ub = ub, ...)
  } else if(type == "laplace"){
    univariate_laplace_update(par_curr = par_curr, par_name = par_name, loglik = loglik, logprior = logprior,
                         loglik_grad = loglik_grad, logprior_grad = logprior_grad, damping = damping,
                         logprior_hessian = logprior_hessian, loglik_hessian = loglik_hessian, lb = lb, ub = ub, ...)
  }
}

#' Univariate updates based on Laplace approximation
#'
#' Following Rue and Held (2004), we perform
#' a random walk based on a normal approximation
#' at the \code{cur} parameter value. This requires both
#' the gradient and the hessian of the log-likelihood and log-posterior;
#' if the latter are left undefined, they are obtained using numerical
#' differentiation via the \code{numDeriv} package.
#' In such cases, users cannot have a formal argument of the functions matching those of
#' \code{grad} and \code{hessian} (typically \code{x} or \code{method})
#'
#' The algorithm uses a random walk proposal
#' and a Metropolis acceptance step.
#'
#' @param cur [double] curr value of the parameter
#' @param par_name [string] the name of the argument for the function
#' @param loglik [function] log likelihood function
#' @param loglik_gradient [function] derivative of the log likelihood with respect to parameter of interest
#' @param loglik_hessian [function] second derivative of the log likelihood with respect to the parameter of interest
#' @param logprior [function] log prior of parameter
#' @param logprior_gradient [function] derivative of the log prior with respect to parameter of interest
#' @param logprior_hessian [function] second derivative of the log prior with respect to the parameter of interest
#' @param lb [scalar] lower bound for the parameter
#' @param ub [scalar] upper bound for the parameter
#' @param damping [scalar] contraction factor for the Newton update
#' @param ... additional arguments passed to the log likelihood function and its derivatives
#' @return a new value for the parameter of interest
#' @export
#' @keywords internal
#' @references Rue, H. and L. Held (2005), Gaussian Markov random fields, CRC press, section 4.4.1
univariate_laplace_update <-
  function(par_curr,
           par_name,
           loglik,
           loglik_grad = NULL,
           loglik_hessian = NULL,
           logprior,
           logprior_grad = NULL,
           logprior_hessian = NULL,
           lb = -Inf,
           ub = Inf,
           damping = 1,
           ...){
    ellipsis <- list(...)
    args <- formalArgs(loglik)
    stopifnot(length(par_curr) == 1L,
              is.finite(par_curr),
              length(lb) == 1,
              length(ub) == 1,
              lb < ub,
              length(damping) == 1L,
              damping > 0,
              damping <= 1,
              par_name %in% args
    )
    if(is.null(logprior_grad)){
       if(isTRUE(any(c("x","method","method.args","side") %in% args))){
          stop("Some arguments in \"logprior\" match those of \"numDeriv\" functions.")
       }
       logprior_grad <- function(...){
          logprior_alt <- function(x, ...){
             args <- list(...)
             args[[par_name]] <- x
             do.call(logprior, args = args)
          }
          numDeriv::grad(func = logprior_alt, x = args[[par_name]], ...)
       }
    } else{
       if(!isTRUE(all(args == formalArgs(logprior_grad)))){
          stop("Invalid signature for \"logprior_grad\".")
       }
    }
    if(is.null(loglik_grad)){
       if(isTRUE(any(c("x","method","method.args","side") %in% args))){
          stop("Some arguments in \"loglik\" match those of \"numDeriv\" functions.")
       }
       loglik_grad <- function(...){
          loglik_alt <- function(x, ...){
             args <- list(...)
             args[[par_name]] <- x
             do.call(loglik, args = args)
          }
          numDeriv::grad(func = loglik_alt, x = args[[par_name]], ...)
       }
    } else{
       if(!isTRUE(all(args == formalArgs(loglik_grad)))){
          stop("Invalid signature for \"loglik_grad\".")
       }
    }
    # Check for Hessians
    if(is.null(logprior_hessian)){
       if(isTRUE(any(c("x","method","method.args","side") %in% args))){
          stop("Some arguments in \"logprior\" match those of \"numDeriv\" functions.")
       }
       logprior_hessian <- function(...){
          logprior_alt <- function(x, ...){
             args <- list(...)
             args[[par_name]] <- x
             do.call(logprior, args = args)
          }
          numDeriv::grad(func = logprior_alt, x = args[[par_name]], ...)
       }
    } else{
       if(!isTRUE(all(args == formalArgs(logprior_hessian)))){
          stop("Invalid signature for \"logprior_hessian\".")
       }
    }
    if(is.null(loglik_hessian)){
       if(isTRUE(any(c("x","method","method.args","side") %in% args))){
          stop("Some arguments in \"loglik\" match those of \"numDeriv\" functions.")
       }
       loglik_hessian <- function(...){
          loglik_alt <- function(x, ...){
             args <- list(...)
             args[[par_name]] <- x
             do.call(loglik, args = args)
          }
          numDeriv::hessian(func = loglik_alt, x = args[[par_name]], ...)
       }
    } else{
       if(!isTRUE(all(args == formalArgs(loglik_hessian)))){
          stop("Invalid signature for \"loglik_hessian\".")
       }
    }
    # Parameter value to be returned
    par_new <- par_curr # if all things fail
    # Copy arguments into a list to call function
    loglik_args <- ellipsis[names(ellipsis) %in% args]
    # Override value of parameter
    loglik_args[[par_name]] <- par_curr
    logpost_grad_curr <-
      do.call(what = loglik_grad,
              args = loglik_args) +
      do.call(what = logprior_grad,
              args = loglik_args)
    logpost_hessian_curr <-
      do.call(what = loglik_hessian,
              args = loglik_args) +
      do.call(what = logprior_hessian,
              args = loglik_args)
    mean_curr <- par_curr -
      damping * logpost_grad_curr/logpost_hessian_curr
    precision_curr <- -logpost_hessian_curr
    if(!isTRUE(precision_curr > 0)){
      return(par_new)
    }
    par_prop <- rtnorm(n = 1,
                       a = lb,
                       b = ub,
                       mean = mean_curr,
                       sd = sqrt(1/precision_curr))
    logpost_curr <-
      do.call(what = loglik,
              args = loglik_args) +
      do.call(what = logprior,
              args = loglik_args)
    # Reverse move
    loglik_args[[par_name]] <- par_prop
    # TODO determine whether we allow users to pass
    # other arguments to logprior
    logpost_grad_prop <-
      do.call(what = loglik_grad,
              args = loglik_args) +
      do.call(what = logprior_grad,
              args = loglik_args)
    logpost_hessian_prop <-
      do.call(what = loglik_hessian,
              args = loglik_args) +
      do.call(what = logprior_hessian,
              args = loglik_args)
    mean_prop <- par_prop -
      damping * logpost_grad_prop/logpost_hessian_prop
    precision_prop <- -logpost_hessian_prop
    if(!isTRUE(precision_prop > 0)){
      return(par_new)
    }
    logpost_prop <-
      do.call(what = loglik,
              args = loglik_args) +
      do.call(what = logprior,
              args = loglik_args)
    log_MH_ratio <-
      logpost_prop - logpost_curr +
      dtnorm(par_curr,
             a = lb,
             b = ub,
             mean = mean_prop,
             sd = sqrt(1/precision_prop),
             log = TRUE) -
      dtnorm(par_prop,
             a = lb,
             b = ub,
             mean = mean_curr,
             sd = sqrt(1/precision_curr),
             log = TRUE)
    if(log_MH_ratio > log(runif(1))){
      par_new <- par_prop
    }
    return(par_new)
  }


#' Multivariate updates based on Laplace approximation
#'
#' Following Rue and Held (2004), we perform
#' a random walk based on a multivariate normal approximation
#' at the \code{cur} parameter value.
#'
#' @param cur [double] current value of the vector parameter
#' @param par_name [string] the name of the argument for the function
#' @param loglik [function] log likelihood function
#' @param loglik_gradient [function] derivative of the log likelihood with respect to parameter of interest
#' @param loglik_hessian [function] hessian matrix, the second derivative of the log likelihood with respect to the parameter of interest
#' @param logprior [function] log prior of parameter
#' @param logprior_gradient [function] derivative of the log prior with respect to parameter of interest
#' @param logprior_hessian [function] second derivative of the log prior with respect to the parameter of interest
#' @param lb [double] vector of lower bounds for the parameters
#' @param ub [double] vector of upper bounds for the parameters
#' @param damping [double] double or vector of contraction factor for the Newton update
#' @param ... additional arguments passed to the log likelihood function and its derivatives
#' @return a new vector for the parameter of interest
#' @export
#' @references Rue, H. and L. Held (2005), Gaussian Markov random fields, CRC press, section 4.4.2
multivariate_laplace_update <-
  function(par_curr,
           par_name,
           loglik,
           loglik_grad = NULL,
           loglik_hessian = NULL,
           logprior,
           logprior_grad = NULL,
           logprior_hessian = NULL,
           lb = -Inf,
           ub = Inf,
           damping = 1,
           ...){
    ellipsis <- list(...)
    args <- formalArgs(loglik)
    stopifnot(length(par_curr) == 1L,
              is.finite(par_curr),
              length(lb) == 1,
              length(ub) == 1,
              lb < ub,
              length(damping) == 1L,
              damping > 0,
              damping <= 1,
              par_name %in% args
    )
    # Parameter value to be returned
    par_new <- par_curr # if all things fail
    # Copy arguments into a list to call function
    loglik_args <- ellipsis[names(ellipsis) %in% args]
    # Check if loglik_grad and logprior_grad are provided
    if(is.null(logprior_grad)){
       if(isTRUE(any(c("x","method","method.args","side") %in% args))){
          stop("Some arguments in \"logprior\" match those of \"numDeriv\" functions.")
       }
       logprior_grad <- function(...){
          logprior_alt <- function(x, ...){
             args <- list(...)
             args[[par_name]] <- x
             do.call(logprior, args = args)
          }
          numDeriv::grad(func = logprior_alt, x = args[[par_name]], ...)
       }
    } else{
       if(!isTRUE(all(args == formalArgs(logprior_grad)))){
          stop("Invalid signature for \"logprior_grad\".")
       }
    }
    if(is.null(loglik_grad)){
       if(isTRUE(any(c("x","method","method.args","side") %in% args))){
          stop("Some arguments in \"loglik\" match those of \"numDeriv\" functions.")
       }
       loglik_grad <- function(...){
          loglik_alt <- function(x, ...){
             args <- list(...)
             args[[par_name]] <- x
             do.call(loglik, args = args)
          }
          numDeriv::grad(func = loglik_alt, x = args[[par_name]], ...)
       }
    } else{
       if(!isTRUE(all(args == formalArgs(loglik_grad)))){
          stop("Invalid signature for \"loglik_grad\".")
       }
    }
    # Check for Hessians
    if(is.null(logprior_hessian)){
       if(isTRUE(any(c("x","method","method.args","side") %in% args))){
          stop("Some arguments in \"logprior\" match those of \"numDeriv\" functions.")
       }
       logprior_hessian <- function(...){
          logprior_alt <- function(x, ...){
             args <- list(...)
             args[[par_name]] <- x
             do.call(logprior, args = args)
          }
          numDeriv::grad(func = logprior_alt, x = args[[par_name]], ...)
       }
    } else{
       if(!isTRUE(all(args == formalArgs(logprior_hessian)))){
          stop("Invalid signature for \"logprior_hessian\".")
       }
    }
    if(is.null(loglik_hessian)){
       if(isTRUE(any(c("x","method","method.args","side") %in% args))){
          stop("Some arguments in \"loglik\" match those of \"numDeriv\" functions.")
       }
       loglik_hessian <- function(...){
          loglik_alt <- function(x, ...){
             args <- list(...)
             args[[par_name]] <- x
             do.call(loglik, args = args)
          }
          numDeriv::hessian(func = loglik_alt, x = args[[par_name]], ...)
       }
    } else{
       if(!isTRUE(all(args == formalArgs(loglik_hessian)))){
          stop("Invalid signature for \"loglik_hessian\".")
       }
    }
    # Override value of parameter
    loglik_args[[par_name]] <- par_curr
    logpost_grad_curr <-
      do.call(what = loglik_grad,
              args = loglik_args) +
      do.call(what = logprior_grad,
              args = loglik_args)
    logpost_hessian_curr <-
      do.call(what = loglik_hessian,
              args = loglik_args) +
      do.call(what = logprior_hessian,
             args = loglik_args)
    eigen_precision <- eigen(-logpost_hessian_curr)
    if(!isTRUE(all(eigen_precision$values > 0))){
      return(par_new)
    }
    covar_curr <- tcrossprod(eigen_precision$vectors %*% diag(1/sqrt(eigen_precision$values)))
    mean_curr <- par_curr -
      damping * logpost_grad_curr %*% covar_curr
    precision_curr <- -logpost_hessian_curr
    par_prop <- TruncatedNormal::rtmvnorm(
      n = 1,
      lb = lb,
      ub = ub,
      mu = mean_curr,
      sigma = covar_curr)
    logpost_curr <-
      do.call(what = loglik,
              args = loglik_args) +
      do.call(what = logprior,
              args = loglik_args)
    # Reverse move
    loglik_args[[par_name]] <- par_prop
    logpost_grad_prop <-
      do.call(what = loglik_grad,
              args = loglik_args) +
      do.call(what = logprior_grad,
              args = loglik_args)
    logpost_hessian_prop <-
      do.call(what = loglik_hessian,
              args = loglik_args) +
      do.call(what = logprior_hessian,
              args = loglik_args)
    eigen_precision <- eigen(-logpost_hessian_prop)
    if(!isTRUE(all(eigen_precision$values > 0))){
      return(par_new)
    }
    covar_prop <- tcrossprod(eigen_precision$vectors %*% diag(1/sqrt(eigen_precision$values)))
    mean_prop <- as.numeric(par_prop -
      damping * logpost_grad_prop %*% covar_prop)
    logpost_prop <-
      do.call(what = loglik,
              args = loglik_args) +
      do.call(what = logprior,
              args = loglik_args)
    log_MH_ratio <-
      logpost_prop - logpost_curr +
      TruncatedNormal::dtmvnorm(par_curr,
             lb = lb,
             ub = ub,
             mu = mean_prop,
             sigma = covar_prop,
             log = TRUE) -
      TruncatedNormal::dtmvnorm(par_prop,
             lb = lb,
             ub = ub,
             mu = mean_curr,
             sigma = covar_curr,
             log = TRUE)
    if(log_MH_ratio > log(runif(1))){
      par_new <- par_prop
    }
    return(par_new)
  }


#' Univariate updates based on MALA
#'
#' @keywords internal
#' @export
univariate_mala_update <-
  function(par_curr,
           par_name,
           loglik,
           loglik_grad = NULL,
           logprior,
           logprior_grad = NULL,
           sd_prop,
           lb = -Inf,
           ub = Inf,
           damping = 1,
           transform = FALSE,
           ...){
    ellipsis <- list(...)
    args <- formalArgs(loglik)
    stopifnot(length(par_curr) == 1L,
              is.finite(par_curr),
              length(lb) == 1,
              length(ub) == 1,
              lb < ub,
              length(damping) == 1L,
              damping > 0,
              damping <= 1,
              par_name %in% args,
              is.logical(transform),
              length(transform) == 1L
              )
    if(lb == -Inf & ub == Inf){
      transform <- FALSE
    }
# Parameter value to be returned
par_new <- par_curr # if all things fail
# Copy arguments into a list to call function
loglik_args <- ellipsis[names(ellipsis) %in% args]

# Check if loglik_grad and logprior_grad are provided
if(is.null(logprior_grad)){
   if(isTRUE(any(c("x","method","method.args","side") %in% args))){
      stop("Some arguments in \"logprior\" match those of \"numDeriv\" functions.")
   }
   logprior_grad <- function(...){
      logprior_alt <- function(x, ...){
         args <- list(...)
         args[[par_name]] <- x
         do.call(logprior, args = args)
      }
      numDeriv::grad(func = logprior_alt, x = args[[par_name]], ...)
   }
} else{
   if(!isTRUE(all(args == formalArgs(logprior_grad)))){
      stop("Invalid signature for \"logprior_grad\".")
   }
}
if(is.null(loglik_grad)){
   if(isTRUE(any(c("x","method","method.args","side") %in% args))){
      stop("Some arguments in \"loglik\" match those of \"numDeriv\" functions.")
   }
   loglik_grad <- function(...){
      loglik_alt <- function(x, ...){
         args <- list(...)
         args[[par_name]] <- x
         do.call(loglik, args = args)
      }
      numDeriv::grad(func = loglik_alt, x = args[[par_name]], ...)
   }
} else{
   if(!isTRUE(all(args == formalArgs(loglik_grad)))){
      stop("Invalid signature for \"loglik_grad\".")
   }
}

# Override value of parameter
loglik_args[[par_name]] <- par_curr
logpost_curr <-
  do.call(what = loglik,
          args = loglik_args) +
  do.call(what = logprior,
          args = loglik_args)

if(!transform){
  logpost_grad_curr <-
    do.call(what = loglik_grad,
            args = loglik_args) +
    do.call(what = logprior_grad,
            args = loglik_args)
  mean_curr <- par_curr +
    damping * 0.5 * logpost_grad_curr
  par_prop <- rtnorm(n = 1,
                     a = lb,
                     b = ub,
                     mean = mean_curr,
                     sd = sd_prop)
  # Reverse move
  loglik_args[[par_name]] <- par_prop
  logpost_grad_prop <-
    do.call(what = loglik_grad,
            args = loglik_args) +
    do.call(what = logprior_grad,
            args = loglik_args)
  mean_prop <- par_prop +
    damping * 0.5 * logpost_grad_prop
  logpost_prop <-
    do.call(what = loglik,
            args = loglik_args) +
    do.call(what = logprior,
            args = loglik_args)
  log_MH_ratio <-
    logpost_prop - logpost_curr +
    dtnorm(par_curr,
           a = lb,
           b = ub,
           mean = mean_prop,
           sd = sd_prop,
           log = TRUE) -
    dtnorm(par_prop,
           a = lb,
           b = ub,
           mean = mean_curr,
           sd = sd_prop,
           log = TRUE)
} else {# TRANSFORM
  tpar_curr <- transfo(par_curr, lb = lb, ub = ub)
  jac_curr <- jac_inv_transfo(tpar_curr, lb = lb, ub = ub, log = FALSE)
  logpost_grad_tcurr <-
    (do.call(what = loglik_grad,
            args = loglik_args) +
    do.call(what = logprior_grad,
            args = loglik_args)) * jac_curr +
    dlogjac_inv_transfo(tpar = tpar_curr, lb = lb, ub = ub)
  tmean_curr <- tpar_curr + damping * 0.5 * logpost_grad_tcurr
  tpar_prop <- rnorm(n = 1,
                    mean = tmean_curr,
                    sd = sd_prop)
  # Reverse move
  par_prop <- inv_transfo(tpar = tpar_prop, lb = lb, ub = ub)
  loglik_args[[par_name]] <- par_prop
  logpost_prop <-
    do.call(what = loglik,
            args = loglik_args) +
    do.call(what = logprior,
            args = loglik_args)
  jac_prop <- jac_inv_transfo(tpar = tpar_prop, lb = lb, ub = ub)
  logpost_grad_tprop <-
    (do.call(what = loglik_grad,
            args = loglik_args) +
    do.call(what = logprior_grad,
            args = loglik_args)) * jac_prop +
      dlogjac_inv_transfo(tpar = tpar_prop, lb = lb, ub = ub)
  tmean_prop <- tpar_prop + damping * 0.5 * logpost_grad_tprop
  log_MH_ratio <-
    logpost_prop - logpost_curr +
    log(jac_prop) - log(jac_curr) +
    dnorm(tpar_curr,
           mean = tmean_prop,
           sd = sd_prop,
           log = TRUE) -
    dnorm(tpar_prop,
           mean = tmean_curr,
           sd = sd_prop,
           log = TRUE)
}
if(log_MH_ratio > log(runif(1))){
  par_new <- par_prop
  attr(par_new, "accept") <- TRUE
} else{
  attr(par_new, "accept") <- FALSE
}

return(par_new)
}

#' Transform parameter to unconstrained scale
#'
#' @param par scalar parameter
#' @param lb lower bound
#' @param ub upper bound
#' @reference Section 56 of the Stan Reference manual version 2.9 at \url{https://github.com/stan-dev/stan/releases/download/v2.9.0/stan-reference-2.9.0.pdf}
transfo <- function(par, lb, ub){
  stopifnot(length(par) == 1L,
            length(lb) == 1L,
            length(ub) == 1L,
            isTRUE(lb < ub))
  if(lb == -Inf & ub == Inf){
    return(par)
  } else if(lb > -Inf & ub == Inf){
    return(log(par - lb))
  } else if(lb == -Inf & ub < Inf){
    return(log(ub - par))
  } else if(lb > -Inf & ub < Inf){
    return(qlogis((par - lb) / (ub - lb)))
  }
}

jac_inv_transfo <- function(tpar, lb, ub, log = FALSE){
  if(lb == -Inf & ub == Inf){
    ljac <- log(tpar)
  }
  if(lb > -Inf & ub == Inf){
    ljac <- tpar
  } else if(lb == -Inf & ub < Inf){
    ljac <- tpar
  } else if(lb > -Inf & ub < Inf){
    ljac <- log(ub - lb) + plogis(tpar, log.p = TRUE) + plogis(tpar, log.p = TRUE, lower.tail = FALSE)
  }
  if(log){
    return(ljac)
  } else{
    return(exp(ljac))
  }
}

dlogjac_inv_transfo <- function(tpar, lb, ub){
  if(lb == -Inf & ub == Inf){
    return(0)
  }
  if(lb > -Inf & ub == Inf){
    return(1)
  } else if(lb == -Inf & ub < Inf){
    return(1)
  } else{
   -1 + 2*plogis(-tpar)
  }
}

inv_transfo <- function(tpar, lb, ub){
  stopifnot(length(tpar) == 1L,
            length(lb) == 1L,
            length(ub) == 1L,
            isTRUE(lb < ub))
  if(lb == -Inf & ub == Inf){
    return(tpar)
  } else if(lb > -Inf & ub == Inf){
    return(exp(tpar) + lb)
  } else if(lb == -Inf & ub < Inf){
    return(ub - exp(tpar))
  } else{
    return(lb + (ub - lb)*plogis(tpar))
  }
}

#' Univariate updates based on random walk
#'
#' The algorithm uses a random walk proposal.
#'
#' @inheritParams univariate_laplace_update
#' @return a new value for the parameter of interest
#' @export
#' @keywords internal
univariate_rw_update <-
  function(par_curr,
           par_name,
           loglik,
           logprior,
           sd_prop,
           lb = -Inf,
           ub = Inf,
           transform = FALSE,
           ...){
    ellipsis <- list(...)
    args <- formalArgs(loglik)
    stopifnot(length(par_curr) == 1L,
              is.finite(par_curr),
              length(lb) == 1,
              length(ub) == 1,
              lb < ub)
    if(lb == -Inf & ub == Inf){
      transform <- FALSE
    }
    # Parameter value to be returned
    par_new <- par_curr # if all things fail
    # Copy arguments into a list to call function
    loglik_args <- ellipsis[names(ellipsis) %in% args]
    # Override value of parameter
    loglik_args[[par_name]] <- par_curr
    logpost_curr <-
      do.call(what = loglik,
              args = loglik_args) +
      do.call(what = logprior,
              args = loglik_args)
    if(!transform){
    par_prop <- rtnorm(n = 1,
                       a = lb,
                       b = ub,
                       mean = par_curr,
                       sd = sd_prop)
    adj <- dtnorm(par_curr,
                  a = lb,
                  b = ub,
                  mean = par_prop,
                  sd = sd_prop,
                  log = TRUE) -
      dtnorm(par_prop,
             a = lb,
             b = ub,
             mean = par_curr,
             sd = sd_prop,
             log = TRUE)
    } else{ # Transformation is TRUE
      # Y = f(X)
      # Transform parameter to the real line
      trpar_curr <- transfo(par = par_curr, lb = lb, ub = ub)
      # Sample a proposal from the normal
      trpar_prop <- rnorm(n = 1, mean = trpar_curr, sd = sd_prop)
      # Compute the difference of log-jacobians
      adj <- jac_inv_transfo(tpar = trpar_prop, lb = lb, ub = ub, log = TRUE) -
        jac_inv_transfo(tpar = trpar_curr, lb = lb, ub = ub, log = TRUE)
      par_prop <- inv_transfo(tpar = trpar_prop, lb = lb, ub = ub)
    }
    # Reverse move
    loglik_args[[par_name]] <- par_prop
    logpost_prop <-
      do.call(what = loglik,
              args = loglik_args) +
      do.call(what = logprior,
              args = loglik_args)
    log_MH_ratio <-
      logpost_prop - logpost_curr + adj
    if(log_MH_ratio > log(runif(1))){
      par_new <- par_prop
    }
    return(par_new)
  }

#' Update function based on Metropolis-Hasting algorithm
#'
#' Proposals are made based on a (conditional) truncated Gaussian distribution at indices \code{ind}
#' given other components.
#' @param cur current value of the vector of parameters of which one component is to be updated.
#' @param lb lower bounds for the parameter of interest. Default to \code{-Inf} if argument is missing.
#' @param ub upper bounds for the parameters of interest. Default to \code{Inf} if argument is missing.
#' @param ind indices of \code{cur} to update.
#' @param prior.fun log prior function, a function of the parameter to update.
#' @param lik.fun log-likelihood function, a function of the parameters to update. The user who wants to keep
#' results of intermediate can pass them via attributes.
#' @param ll value of the log-likelihood at \code{cur}.
#' @param pmu proposal mean; if missing, default to random walk.
#' @param pcov covariance matrix of the proposal for \code{cur}.
#' @param cond logical; should updates be made conditional on other values in \code{cur}. Default to \code{TRUE}.
#' @param transform logical; should individual updates for each parameter be performed on the unconstrained space? If \code{TRUE},
#' argument \code{cond} is ignored.
#' @param ... additional arguments passed to the function, currently ignored.
#' @return a list with components
#' \itemize{
#' \item{\code{ll}}: value of the log-likelihood, with potential additional parameters passed as attributes;
#' \item{\code{cur}}: values of the parameters after update;
#' \item{\code{accept}}: logical indicating the proposal has been accepted (\code{TRUE}) or rejected (\code{FALSE}).
#' }
#' @export
mh.fun <- function(cur, lb, ub, prior.fun, lik.fun, ll, ind, pmu, pcov, cond = TRUE, transform = FALSE, ...) {
  # overwrite missing values with default setting
  if (missing(ind)) {
    lc <- length(cur)
    ind <- 1:lc
  } else {
    lc <- length(ind)
  }
  if (missing(lb)) {
    lb <- rep(-Inf, length.out = lc)
  }
  if (missing(ub)) {
    ub <- rep(Inf, length.out = lc)
  }
  # Transform to different scale if transfo is TRUE, see STAN manual for transformations
  # This help with skewed variables
  # Sanity checks for length of arguments - only lengths
  pcov <- as.matrix(pcov) # handle scalar case
  stopifnot(length(lb) == lc, lc == length(ub), ncol(pcov) == nrow(pcov), ncol(pcov) == length(cur))

  if (transform) {
    propRW <- rprop.rwmh(cur = cur[ind], cov = as.matrix(pcov[ind, ind]), lbound = lb, ubound = ub)
    prop <- cur
    prop[ind] <- propRW$prop
    jac <- propRW$logjac
  } else {
    # Copy value
    prop <- cur
    if (missing(pmu)) {
      pmu <- cur
      RW <- TRUE
    } else {
      stopifnot(length(pmu) == length(cur))
    }

    if (!cond) {
      sig <- as.matrix(pcov[ind, ind])
      # Sample new proposal from truncated Normal centered at current value
      prop[ind] <- TruncatedNormal::mvrandn(l = lb, mu = pmu[ind], u = ub, Sig = sig, n = 1)
      if (RW) {
        # mean is at previous iteration, so density term drops out because of symmetry,
        # but not normalizing constant of truncated dist (since means differ)
        jac <- suppressWarnings(-log(TruncatedNormal::mvNcdf(l = lb - prop[ind], u = ub - prop[ind], Sig = sig, n = 1000)$prob) +
                                  log(TruncatedNormal::mvNcdf(l = lb - cur[ind], u = ub - cur[ind], Sig = sig, n = 1000)$prob))
      } else { # fixed mean, so normalizing constants for truncated components cancels out of the ratio of transition kernel densities
        jac <- mgp::dmvnorm(x = cur[ind], mean = pmu[ind], sigma = sig, logd = TRUE) -
          mgp::dmvnorm(x = prop[ind], mean = pmu[ind], sigma = sig, logd = TRUE)
      }
    } else { # cond == TRUE
      prop[ind] <- rcondmvtnorm(n = 1L, ind = ind, x = cur[-ind], lbound = lb, ubound = ub, mu = pmu, sigma = pcov, model = "norm")
      jac <- dcondmvtnorm(n = 1L, ind = ind, x = cur, lbound = lb, ubound = ub, mu = pmu, sigma = pcov, model = "norm", log = TRUE) -
        dcondmvtnorm(n = 1L, ind = ind, x = prop, lbound = lb, ubound = ub, mu = pmu, sigma = pcov, model = "norm", log = TRUE)
    }
  }
  ll.p <- lik.fun(prop)
  if (missing(ll)) {
    ll.c <- lik.fun(cur)
  } else {
    ll.c <- ll
  }
  prior.p <- prior.fun(prop)
  prior.c <- prior.fun(cur)
  acpt <- as.vector(ll.p - ll.c + prior.p - prior.c + jac)
  if (isTRUE(log(runif(1)) < acpt)) {
    return(list(ll = ll.p, cur = prop, accept = TRUE))
  } else {
    return(list(ll = ll.c, cur = cur, accept = FALSE))
  }
}


#' Variance adaptation
#'
#' Adapt standard deviation of a proposal for a Markov chain Monte Carlo algorithm, based on the acceptance rate.
#' The function is targeting an acceptance probability of \code{target} for univariate Metropolis--Hastings proposals, which is
#' optimal for certain Gaussian regimes. If the number of attempts is large enough and the acceptance rate is not close
#' to the target, the standard deviation of the proposal will be increased or reduced and the number acceptance and attempts
#' reinitialized. If no change is made, the components \code{acc} will be the same as acceptance and similarly for attempts.
#'
#' @param attempts integer indicating number of attempts
#' @param acceptance integer giving the number of moves accepted by the algorithm
#' @param sd.p standard deviation of proposal
#' @param target target acceptance rate
#' @return a list with components \code{sd}, \code{acc}, \code{att}
#' @export
adaptive <- function(attempts, acceptance, sd.p, target = 0.234) {
  stopifnot(sd.p > 0)
  if (attempts < acceptance) {
    stop("Invalid input: the number of attempts must be larger than the number of acceptance")
  }
  att <- attempts
  acc <- acceptance
  newsd <- sd.p
  dist_target <- qlogis(acc / att) - qlogis(target)
  dist_target <- pmax(pmin(3, dist_target), -3)
  multfact <- plogis(dist_target)/0.5
  if (att > 20 & abs(dist_target) > 2) {
      newsd <- sd.p * multfact
      att <- 0L
      acc <- 0L
  }
  if (att > 30 & abs(dist_target) > 1.25) {
      newsd <- sd.p * multfact
      att <- 0L
      acc <- 0L
  }
  if (att > 50) {
    newsd <- sd.p * multfact
    att <- 0L
    acc <- 0L
  }
  return(list(sd = newsd, acc = acc, att = att))
}

