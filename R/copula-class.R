library("R6")
#' Classe copulagev
#' @title Classe da Copulas GEV
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @description A classe representa as copulas e fornece um ferramental para copula
#' @field copula.type An extreme-value copula class within gumbel, galambos, huslerReiss, tawn,
#' @field rho The value of Spearman's dependence index
#' @field theta.u1 Vetor de parâmetro da marginal 1
#' @field theta.u2 Vetor de parâmetro da marginal 2
#' @field name Name for test
#' @section Methods:
#' \describe{
#'  \item{\code{new(cop.type, Spearman.rho)}}{This method is used to create object of this class with \code{cop.type} as copula class and \code{Spearman.rho} the Spearmon's dependency index value.}
#'  \item{\code{initialize(cop.type, Spearman.rho)}}{This method is used to create object of this class with \code{cop.type} as copula class and \code{Spearman.rho} the Spearmon's dependency index value.}
#'  \item{\code{set.class.cop(copula.name)}}{This method is used to set da copula features using \code{copula.name}.}
#'  \item{\code{get.class.cop()}}{This method is used to create object of this class with  Spearmon's dependency index value.}}
copulagev <- R6Class("copulagev",
                    public = list(
                        #' @description Construtor da classe copulagev
                        #' @param cop.type Nome da copula
                        #' @param Spearman.rho Coeficiente de correlação de Spearman
                        initialize = function(cop.type=NULL, Spearman.rho=NULL){
                            self$set.class.cop(cop.type,Spearman.rho)
                        },
                        rho = NULL,
                        copula.type = NULL,
                        theta.u1 = list(mu=NULL,sigma=NULL,xi=NULL,delta=NULL),
                        theta.u2 = list(mu=NULL,sigma=NULL,xi=NULL,delta=NULL),
                        #' @description Retorna o nome da classe
                        name = function(){
                            return(private$class.name)
                        },
                        #' @description Define o tipo de copula extrema
                        #' @param cop.type Nome da copula
                        #' @param Spearman.rho Coeficiente de correlação de Spearman
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
                        #' @description Get the main features of a copula
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
                    private = list(
                        class.copula = NULL,
                        class.name = "copulagev"
                        )
                    )


