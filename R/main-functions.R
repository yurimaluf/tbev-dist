library(R6)
# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   To indicate that it needs to be exposed. run devtools::document()


#' The transformation function that gives the name transformed BEV distribution
#'
#' @param x input of random variable to transformed.
#' @param mu parameter location
#' @param delta distortion parameter
#'
#' @return The value of the probability failure
#' @export
transf.BEV = function(x,mu,delta){
    transf = (x-mu)*abs(x-mu)^delta
    return(transf)
}


#' The function computes the probability of failure given by P(X<Y), where X~F and Y~F are random variables with GEV distribution. It has a special importance in reliability engineering where X and Y represent stress and strength, respectively. The vector parameter (xi,mu,sigma) is given by location mu, scale sigma and shape xi parameter.
#'
#' @param theta1 Vector of parameter (xi,mu,sigma,delta) of the distribuicao 1 F_X
#' @param theta2 Vector of parameter (xi,mu,sigma,delta) of the distribuicao 2 F_Y
#' @param copula.type An extreme-value copula class within gumbel, galambos, huslerReiss, tawn, and tev
#' @param rho The value of Spearman's dependence index
#'
#' @return The value of probability failure
#' @export
failure.tbev = function(copula.object){
    # Integrating data
    param = copula.object$assembly()

    # Figureing out integral
    result=integrate(DC.integrand,lower=0.001,upper=0.999,param)
    return(result$value)
}

#' @export
DC.integrand = function(t,param){
    # Computing the arguments of DC(,)
    x=bgev::qbgev(t, mu=param$theta1$mu, sigma=param$theta1$sigma, xi=param$theta1$xi,delta=param$theta1$delta)
    u=bgev::pbgev(x, mu=param$theta2$mu, sigma=param$theta2$sigma, xi=param$theta2$xi,delta=param$theta2$delta)
    return(tbevdist::D.copula(t,u,param$class.copula))
}

#' @export
D.copula = function(u1,u2,class.copula){
    u2=u2[1]
    return((copula::cCopula(cbind(u2, u1),class.copula))[,2])
}

#' @export
ev.copula = function(u1,u2,class.copula){
    return(copula::pCopula(c(u1,u2),class.copula))
}

#' @export
ev.class.copula = function(type,rho,df=3){
    type.cop = ifelse(tolower(type) %in% c("gumbel","galambos","huslerReiss","tawn","tev"),tolower(type),"gumbel")
    ev.Copula = switch(
        type.cop,
        "gumbel"= copula::gumbelCopula(copula::iRho(copula::gumbelCopula(), rho=rho), dim=2),
        "galambos"= copula::galambosCopula(copula::iRho(copula::galambosCopula(), rho=rho)),
        "huslerReiss"= copula::huslerReissCopula(copula::iRho(copula::huslerReissCopula(), rho=rho)),
        "tawn"= copula::tawnCopula(copula::iRho(copula::tawnCopula(), rho=rho)),
        "tev"= copula::tevCopula(copula::iRho(copula::tevCopula(), rho=rho), df = df)
    )
    return(ev.Copula)
}

#' @export
mount.theta = function(xi,mu,sigma,delta){
    return(list(xi=xi,mu=mu,sigma=sigma,delta=delta))
}









