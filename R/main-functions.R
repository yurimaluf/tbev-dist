Gr = function(t,s,r){
     ##return(exp((-log(t))^r+(-log(s))^r))
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

Cr = function(t,s,r){
    return(((-log(t))^r+(-log(s))^r)^(1/r))
}

hptbev.aux1 = function(u1,u2,r){
   return(Cr(u1,u2,r)*((-log(u1))^(r-1))/u1)
}


failure.integrand = function(u,theta1,theta2,r){
    u1=evd::qgev(u,loc=theta1[2], scale=theta1[3], shape=theta1[1])
    u2=evd::pgev(u1,loc=theta2[2], scale=theta2[3], shape=theta2[1])
    return(Cr(u,u2,r)*exp(-u2)*((-log(u))^(r-1))/u)
}

failure.integrand.2 = function(x,fbg,theta1,theta2,r){
    u1=fbg(x,theta1)
    u2=fbg(x,theta2)
    #hptbev.aux1(u1,u2,r)
    return(Cr(u1,u2,r)*((-log(u1))^(r-1))/u1)
}



#' The function computes the probability of failure given by P(X<Y), where X~F and Y~F are random variables with GEV distribution. It has a special importance in reliability engineering where X and Y represent stress and strength, respectively. The vector parameter (xi,mu,sigma) is given by location mu, scale sigma and shape xi parameter.
#'
#' @param r represents the parameter of dependence function.
#' @param theta1 Vector of parameter (xi,mu,sigma) of the distribuicao 1 F_X
#' @param theta2 Vector of parameter (xi,mu,sigma) of the distribuicao 2 F_Y
#'
#' @return The value of the probability failure
#' @export
failure.tbev = function(r,theta1,theta2){
    # Calculo da Integral
    result=integrate(failure.integrand,lower=0,upper=1,r=r,theta1=theta1,theta2=theta2)
    return(result$value)
}

#' The function computes the probability of failure given by P(X<Y). It has a special importance in reliability engineering where the random variables X and Y represent stress and strength, respectively.
#'
#' @param r represents the parameter of dependence function
#' @param theta1 Vector of parameter (xi,mu,sigma) of the distribuicao 1 F_X
#' @param theta2 Vector of parameter (xi,mu,sigma) of the distribuicao 2 F_Y
#'
#' @return The value of the probability failure
#' @export
#'
failure.tbev.2 = function(r,theta1,theta2){
    # Def. do Intervalo
    intervalo = tbev.interval(theta1, theta2)
    # Calculo da Integral
    result=integrate(failure.integrand.2,lower=intervalo[1],upper=intervalo[2],r=r,theta1=theta1,theta2=theta2)
    return(result$value)
}




