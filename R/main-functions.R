Gr = function(t,s,r){
     ##return(exp((-log(t))^r+(-log(s))^r))
}


Cr = function(t,s,r){
    return(((-log(t))^r+(-log(s))^r)^(1/r))
}

transform.BEV = function(x,mu,delta){
    transf = (x-mu)*abs(x-mu)^delta
    return(transf)
}

failure.integrand = function(u,theta1,theta2,r){
    x=bgev::qbgev(u,mu=theta1[2], sigma=theta1[3], xi=theta1[1],delta=theta1[4])
    u2=bgev::pgev(x,mu=theta2[2], sigma=theta2[3], xi=theta2[1],delta=theta2[4])
    return(Cr(u,u2,r)*exp(-u2)*((-log(u))^(r-1))/u)
}


#' The function computes the probability of failure given by P(X<Y), where X~F and Y~F are random variables with GEV distribution. It has a special importance in reliability engineering where X and Y represent stress and strength, respectively. The vector parameter (xi,mu,sigma) is given by location mu, scale sigma and shape xi parameter.
#'
#' @param r represents the parameter of dependence function.
#' @param theta1 Vector of parameter (xi,mu,sigma,delta) of the distribuicao 1 F_X
#' @param theta2 Vector of parameter (xi,mu,sigma,delta) of the distribuicao 2 F_Y
#'
#' @return The value of the probability failure
#' @export
failure.tbev = function(r,theta1,theta2){
    # Calculo da Integral
    result=integrate(failure.integrand,lower=0,upper=1,r=r,theta1=theta1,theta2=theta2)
    return(result$value)
}





