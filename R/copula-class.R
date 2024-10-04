library("R6")
#' Classe Copulas GEV
#' @title Classe Copulas GEV
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @description The class provides tools to define the type and dependency parameter attributes of a copula, and the associated methods such as its value, conditional probability and failure probability
#' @field copula.type An extreme-value copula class within gumbel, galambos, huslerReiss, tawn,
#' @field rho The value of Spearman's dependence index
#' @field theta.u1 Parameter Vector
#' @field theta.u2 Parameter Vector
#' @field name Name for test
#' @section Methods:
#' \describe{
#'  \item{\code{new(cop.type, Spearman.rho)}}{This method is used to create object of this class with \code{cop.type} as copula class and \code{Spearman.rho} the Spearmon's dependency index value.}
#'  \item{\code{initialize(cop.type, Spearman.rho)}}{This method is used to create object of this class with \code{cop.type} as copula class and \code{Spearman.rho} the Spearmon's dependency index value.}
#'  \item{\code{set.class.cop(copula.name)}}{This method is used to set da copula features using \code{copula.name}.}
#'  \item{\code{get.class.cop()}}{This method is used to create object of this class with  Spearmon's dependency index value.}}
copulagev <- R6Class("copulagev",
                    public = list(
                        #' @description
                        #' Create an object of copulagev class
                        #' @param cop.type Name of the copula
                        #' @param Spearman.rho Coeficient of Spearman correlation
                        initialize = function(cop.type=NULL, Spearman.rho=NULL){
                            # cat(paste0("Hello, my name is ", self$name, ".\n"))
                            self$set.class.cop(cop.type,Spearman.rho)
                        },
                        name = "Yuri",
                        rho = NULL,
                        copula.type = NULL,
                        theta.u1 = list(xi=NULL,mu=NULL,sigma=NULL,delta=NULL),
                        theta.u2 = list(xi=NULL,mu=NULL,sigma=NULL,delta=NULL),
                        #' @description
                        #' Set the main features of a copula
                        #' @param cop.type Name of the copula
                        #' @param Spearman.rho Coeficient of Spearman correlation
                        set.class.cop = function(copula.name,rho.coef){
                            if(!is.null(copula.name) && !is.null(rho.coef)){
                                df=3
                                type.cop = ifelse(tolower(copula.name) %in% c("gumbel","galambos","huslerReiss","tawn","tev"),tolower(copula.name),"gumbel")
                                ev.Copula = switch(
                                    type.cop,
                                    "gumbel"= copula::gumbelCopula(copula::iRho(copula::gumbelCopula(), rho=rho.coef), dim=2),
                                    "galambos"= copula::galambosCopula(copula::iRho(copula::galambosCopula(), rho=rho.coef)),
                                    "huslerReiss"= copula::huslerReissCopula(copula::iRho(copula::huslerReissCopula(), rho=rho.coef)),
                                    "tawn"= copula::tawnCopula(copula::iRho(copula::tawnCopula(), rho=rho.coef)),
                                    "tev"= copula::tevCopula(copula::iRho(copula::tevCopula(), rho=rho.coef), df = df)
                                )
                                private$class.copula = ev.Copula
                            }
                        },
                        #' @description
                        #' Get the main features of a copula
                        get.class.cop = function(){
                            if(is.null(private$class.copula)){
                                self$set.class.cop()
                            }
                            return(private$class.copula)
                        },
                        #' @description
                        #' Assembly all copula data including the type of correlation, the value of dependence coefficient, marginals parameters
                        assembly = function(){
                            param = list(theta1=self$theta.u1,theta2=self$theta.u2,class.copula=private$class.copula)
                            return(param)
                        }
                    ),
                    private = list(class.copula = NULL)
                    )


