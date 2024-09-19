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



#' Title Verifica o intervalo do suporte da distribuição GEV para cada componente
#'
#' @param theta1 Vetor de parametros (xi,mu,sigma2) da distribuicao 1
#' @param theta2 Vetor de parametros (xi,mu,sigma2) da distribuicao 2
#'
#' @return
#' @export
#'
#' @examples
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
        intervalo = sort(domain.1$xi.signal * c(
            max(
                domain.1$x.min * domain.1$xi.signal,
                domain.2$x.min * domain.2$xi.signal
            ),
            Inf
        ))
    }
}

Cr = function(t,s,r){
    return((-log(t))^r+(-log(s))^r)
}

hptbev.aux1 = function(u1,u2,r){
   return(Cr(u1,u2,r)*((-log(u1))^(r-1))/u1)
}

hptbev.aux2 = function(x,fbg,theta1,theta2,r){
    u1=fbg(x,theta1)
    u2=fbg(x,theta2)
    #hptbev.aux1(u1,u2,r)
    return(Cr(u1,u2,r)*((-log(u1))^(r-1))/u1)
}


#' Compute the probability of failure given by P(X<Y). It has a special importance in reliability engineering where the random variables X and Y represent stress and strength, respectively.
#'
#' @param r
#' @param theta1
#' @param theta2
#'
#' @return
#' @export
#'
#' @examples
h.failure.tbev = function(r,theta1,theta2){
    #
    xi1=theta1[1]
    mu1=theta1[2]
    sigma1=theta1[3]
    #
    xi2=theta2[1]
    mu2=theta2[2]
    sigma2=theta2[3]
    #

    #
    integrate(hptbev.aux2,r,theta1,theta2)

}






