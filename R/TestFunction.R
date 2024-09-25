#' The function computes the probability of failure given by P(X<Y). It has a special importance in reliability engineering where the random variables X and Y represent stress and strength, respectively.
#'
#' @param r represents the parameter of dependence function
#' @param theta1 Vector of parameter (xi,mu,sigma) of the distribuicao 1 F_X
#' @param theta2 Vector of parameter (xi,mu,sigma) of the distribuicao 2 F_Y
#'
#' @return The value of the probability failure
#'
failure.tbev.2 = function(r,theta1,theta2){
    # Def. do Intervalo
    intervalo = tbev.interval(theta1, theta2)
    # Calculo da Integral
    result=integrate(failure.integrand.2,lower=intervalo[1],upper=intervalo[2],r=r,theta1=theta1,theta2=theta2)
    return(result$value)
}


#' Distribution function of transformed GEV distribution
#'
#' @param x vector of quantiles
#' @param theta parameter vector (xi,mu,sigma,delta) where xi, mu, sigma and delta are respectively the shape, location and scale parameter and delta is the parameter of transformed function Tx
#'
#' @return vector of probabilities
#'
pbev = function(x,theta){
    y=transform.BEV(x,mu=theta[2],delta=theta[4])
    return(evd::pgev(y,loc=0, scale=theta[3], shape=theta[1]))
}


tbev.domain = function(theta){
    if(theta[1] == 0){
        x.min=0
        xi=0
    }else if(theta[1] > 0){
        x.min=theta[2]-theta[3]/theta[1]
        xi=1
    }else{
        x.min=theta[2]-theta[3]/theta[1]
        xi=-1
    }
    return(list(x.min=x.min, xi.signal=xi))
}

#' Verifica o intervalo do suporte da distribuição GEV para cada componente
#'
#' @param theta1 Vetor de parametros (xi,mu,sigma2) da distribuicao 1
#' @param theta2 Vetor de parametros (xi,mu,sigma2) da distribuicao 2
#'
#' @return gives the support of distribution according to parameters
#' @export
tbev.interval = function(theta1, theta2) {
    domain.1 = tbev.domain(theta1)
    domain.2 = tbev.domain(theta2)
    # Verifica a intersecção entre os intervalos
    if (domain.1$xi.signal * domain.2$xi.signal < 0) {
        #sinais de xi diferentes
        ind.positive = which(max(c(
            domain.1$xi.signal, domain.2$xi.signal
        )) == c(domain.1$xi.signal, domain.2$xi.signal),
        arr.ind = T)
        ind.negative = which(min(c(
            domain.1$xi.signal, domain.2$xi.signal
        )) == c(domain.1$xi.signal, domain.2$xi.signal),
        arr.ind = T)
        start.point.positive = c(domain.1$x.min, domain.2$x.min)[ind.positive]
        start.point.negative = c(domain.1$x.min, domain.2$x.min)[ind.negative]
        if (start.point.positive - start.point.negative > 0) {
            intervalo = c(NULL, NULL)
        } else{
            intervalo = c(start.point.positive, start.point.negative)
        }
    } else {
        #sinais de xi iguais
        intervalo = sort(domain.1$xi.signal * c(max(domain.1$x.min * domain.1$xi.signal,domain.2$x.min * domain.2$xi.signal),Inf))
    }
    return(intervalo)
}


hptbev.aux1 = function(u1,u2,r){
    return(Cr(u1,u2,r)*((-log(u1))^(r-1))/u1)
}

failure.integrand.2 = function(x,fbg,theta1,theta2,r){
    u1=fbg(x,theta1)
    u2=fbg(x,theta2)
    #hptbev.aux1(u1,u2,r)
    return(Cr(u1,u2,r)*((-log(u1))^(r-1))/u1)
}


