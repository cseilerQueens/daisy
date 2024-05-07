################################################################################
#' Normalize parameter values
#' @description This function allows the user to interrupt and resume the optimization
#' 
#' @param object a GA object
#' @return a GA object and a pop.rds file
#'
#' @examples
#' 
#' GA <- ga(type = "real-valued",
#'          fitness =  function(x) -Rastrigin(x[1], x[2]),
#'          lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
#'          popSize = 10, maxiter = 3, seed = 1,
#'          postFitness = postfit)
#' 
#' @export
postfit <- function(object)
{
  pop <- object@population
  # update info
  if(!file.exists("pop.rds")) assign(".pop", NULL, envir = globalenv())
  .pop <- get(".pop", envir = globalenv())
  saveRDS(.pop, "pop.rds")
  .pop <- readRDS("pop.rds")
  assign(".pop", append(.pop, list(pop)), envir = globalenv()) 
  # output the input ga object (this is needed!!)
  saveRDS(object, "object.rds")
  object <- readRDS("object.rds")
  object 
}

