#' @title Quasi-Cauchy quantile regression
#' @name qcreg
#'
#' @description Returns an object of class \code{rq()} that represents a Quasi-Cauchy quantile regression fit. Quasi-Cauchy quantile regression is useful when you want to perform quantile regression analysis on data limited to the unit range.
#'
#'
#' @usage qcreg(formula, data, tau=0.5, npi=100, criterion="bic", tau_i=0.05, tau_f=0.95) 
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms, separated by + operators, on the right.
#' @param data a \code{data.frame()} composed of the variables that will be used in the model.
#' @param tau the quantile to be estimated, this is a number strictly between 0 and 1. The default value is 0.5.
#' @param npi (optional) the number of Pi's that will be considered for choosing the Pi that best fits the model. The default value is 100.
#' @param criterion (optional) criterion to decide the Pi that fits the model. Choose "aic" for AIC, "bic" for BIC and "R2" for pseudo-R2. Or, indicate a numerical value between 0<Pi<pi to use a particular Pi. The default is the automatic choice of Pi following the BIC criterion.
#' @param tau_i (optional) if you want to estimate several quantiles simultaneously, enter the lower limit of the range of coatis you want to estimate here. The default value is 0.05.
#' @param tau_f (optional) if you want to estimate several quantiles simultaneously, enter the upper limit of the range of coatis you want to estimate here. The default value is 0.95.
#'
#'
#'
#' @details The Quasi-Cauchy quantile regression model is based on the traditional quantile model, proposed by Koenker (2005) (<doi:10.1017/CBO9780511754098>), to which the Quasi-Cauchy link function is added, allowing the estimation of quantile regression when modeling a variable of nature limited to the ranges \[0,1\], (0,1\], \[0,1) or (0,1). For more details on Quasi-Cauchy quantile regression, see de Oliveira, Ospina, Leiva, Figueroa-Zuniga and Castro (2023) (<doi:10.3390/fractalfract7090667>).
#'
#'
#' @return \code{qcreg()} returns an object of class \code{rq()}, hence all outputs of an \code{rq()} object are accessible.
#' @return \code{index} returns the Pi value used in estimating the model and 4 goodness-of-fit criteria, namely: AIC, BIC, pseudo-R2, adjusted pseudo-R2. 
#' @return \code{effects} returns the marginal effect on the average.
#' @return \code{quantregplot} returns argument for graphical visualization of estimates (and confidence intervals) considering a range of values for tau instead of a single value.
#' @return \code{pis} returns the values of Pi considered in the procedure for choosing the ideal Pi, as well as the corresponding goodness-of-fit criterion values. Available only when Pi is chosen via goodness-of-fit criteria.
#' @author Jose Sergio Case de Oliveira
#' 
#' 
#' @references \[1\] Koenker, R. W. (2005). Quantile Regression, Cambridge U. Press. <doi:10.1017/CBO9780511754098>
#' @references \[2\] de Oliveira, J.S.C.; Ospina, R.; Leiva, V.; Figueroa-Zuniga, J.; Castro, C. (2023). Quasi-Cauchy Regression Modeling for Fractiles Based on Data Supported in the Unit Interval. Fractal Fract. 7, 667. <doi:10.3390/fractalfract7090667>
#' 
#' 
#' 
#' @examples
#' 
#' data("Democratization", package = "qcauchyreg")
#'
#' fit <- qcreg(democratization ~ schooling + press_freedom, data = Democratization, criterion=1)
#' summary(fit)
#' fit$effects
#' fit$index
#'
#' 
#' 
#' 
#' @examples
#' 
#'
#' data("Poverty", package = "qcauchyreg")
#'
#' fit2 <- qcreg(poverty ~ population + illiteracy + pc_income, 
#' data = Poverty, npi=50, criterion="bic")
#' summary(fit2)
#' fit2$effects
#' fit2$index
#'
#' plot(fit2$pis, type="l")
#'
#' plot(fit2$quantregplot) 
#' 
#' 
#' @export 
utils::globalVariables(c("stats", "AIC", "as.formula", "model.extract", "model.frame"))
qcreg=function(formula, data, tau=0.5, npi=100, criterion="bic", tau_i=0.05, tau_f=0.95){

if(inherits(formula, "formula")==0) stop("The argument formula must be a formula. See the documentation!")
if(inherits(data, "data.frame")==0) stop("The argument data must be a data.frame. See the documentation!")
if(inherits(tau, "numeric")==0) stop("The argument tau must be a numeric. See the documentation!")
if(inherits(npi, "numeric")==0) stop("The argument npi must be a numeric. See the documentation!")
if(inherits(criterion, "character")==1){if(sum(criterion==c("aic","bic","R2"))==0){stop("Invalid character argument. See the documentation!")}}
if(inherits(criterion, "numeric")==1){if(criterion>=pi){stop("Invalid numeric argument. 0<PI<pi. See the documentation!")}}
if(inherits(criterion, "numeric")==1){if(criterion<=0){stop("Invalid numeric argument. 0<PI<pi. See the documentation!")}}
if(inherits(tau_i, "numeric")==0) stop("The argument tau_i must be a numeric. See the documentation!")
if(inherits(tau_f, "numeric")==0) stop("The argument tau_f must be a numeric. See the documentation!")
if(tau_i>tau_f) stop("The argument tau_i must be less than tau_f. See the documentation!")
if(tau_i>0.999) stop("The argument tau_i must be less than 1 and greater than 0. See the documentation!")
if(tau_i<0.001) stop("The argument tau_i must be less than 1 and greater than 0. See the documentation!")
if(tau_f>0.999) stop("The argument tau_f must be less than 1 and greater than 0. See the documentation!")
if(tau_f<0.001) stop("The argument tau_i must be less than 1 and greater than 0. See the documentation!")
HH=model.frame(formula,data=data)
Y=as.matrix(model.extract(HH, "response"))
if(max(Y)>1) stop("The response variable must be in the unit range. See the documentation!")
if(max(Y)<0) stop("The response variable must be in the unit range. See the documentation!")
tau=tau
if(inherits(criterion, "numeric")==1){PI_ext=criterion}
nrep=npi
PI=seq(0.001,I(pi-0.001),length=nrep)
n=dim(HH)[1]
if(inherits(criterion, "numeric")==0){
A=matrix(0,nrep,4)
for(i in 1:nrep){
G = function(x){I(tan(PI[i]*((x)-0.5)))}
YY = G(Y)
HHH=cbind(HH,YY)
form2 = as.formula(paste("YY~",paste0(formula)[3])) 
fe = quantreg::rq(form2, tau=tau, data=HHH)
fe0 = quantreg::rq(YY~1, tau=tau, data=HHH)
R2 = 1 - fe$rho/fe0$rho
A[i,1]=PI[i]
A[i,2]=R2
A[i,3]=AIC(fe)[1]
A[i,4]=AIC(fe, k=log(n))[1]}
colnames(A) = c("PI","R2","AIC","BIC")
A=data.frame(A)
PI_o = A$PI[which.max(A$R2)]
R_o = A[which.max(A[,2]),2]
PI_o2 = A$PI[which.min(A$AIC)]
AIC_o = A[which.min(A[,3]),3]
PI_o3 = A$PI[which.min(A$BIC)]
BIC_o = A[which.min(A[,4]),4]
Ind = rbind(c(round(R_o,6), PI_o),c(round(AIC_o,6), PI_o2), c(round(BIC_o,6), PI_o3))
colnames(Ind) = c("Index", "PI")
rownames(Ind) = c("Pseudo-R2","AIC","BIC")}
if(criterion=="R2"){ 
G = function(x){I(tan(PI_o*((x)-0.5)))} 
Pp=cbind(A$PI, A$R2)
colnames(Pp) = c("PI", "Pseudo-R2")
YY = G(Y)
HHH=cbind(HH,YY)
form3 = as.formula(paste("YY~",paste0(formula)[3]))
gy_o <- quantreg::rq(form3, tau=tau, data=HHH)
fe0 = quantreg::rq(YY~1, tau=tau, data=HHH)
R2 = 1 - gy_o$rho/fe0$rho
AA=matrix(0,1,5)
colnames(AA) = c("PI","Pseudo-R2","Adjust. pseudo-R2","AIC","BIC")
rownames(AA) = c("Values")
AA[1,1]=PI_o
AA[1,2]=R2
AA[1,3]=1-(1-R2)*(n-1)/(n-dim(HH)[2])
AA[1,4]=AIC(gy_o)[1]
AA[1,5]=AIC(gy_o, k=log(n))[1]}
if(criterion=="aic"){ 
G = function(x){I(tan(PI_o2*((x)-0.5)))} 
Pp=cbind(A$PI, A$AIC)
colnames(Pp) = c("PI", "AIC")
YY = G(Y)
HHH=cbind(HH,YY)
form3 = as.formula(paste("YY~",paste0(formula)[3]))
gy_o <- quantreg::rq(form3, tau=tau, data=HHH)
fe0 = quantreg::rq(YY~1, tau=tau, data=HHH)
R2 = 1 - gy_o$rho/fe0$rho
AA=matrix(0,1,5)
colnames(AA) = c("PI","Pseudo-R2","Adjust. pseudo-R2","AIC","BIC")
rownames(AA) = c("Values")
AA[1,1]=PI_o2
AA[1,2]=R2
AA[1,3]=1-(1-R2)*(n-1)/(n-dim(HH)[2])
AA[1,4]=AIC(gy_o)[1]
AA[1,5]=AIC(gy_o, k=log(n))[1]}
if(criterion=="bic"){ 
G = function(x){I(tan(PI_o3*((x)-0.5)))} 
Pp=cbind(A$PI, A$BIC)
colnames(Pp) = c("PI", "BIC")
YY = G(Y)
HHH=cbind(HH,YY)
form3 = as.formula(paste("YY~",paste0(formula)[3]))
gy_o <- quantreg::rq(form3, tau=tau, data=HHH)
fe0 = quantreg::rq(YY~1, tau=tau, data=HHH)
R2 = 1 - gy_o$rho/fe0$rho
AA=matrix(0,1,5)
colnames(AA) = c("PI","Pseudo-R2","Adjust. pseudo-R2","AIC","BIC")
rownames(AA) = c("Values")
AA[1,1]=PI_o3
AA[1,2]=R2
AA[1,3]=1-(1-R2)*(n-1)/(n-dim(HH)[2])
AA[1,4]=AIC(gy_o)[1]
AA[1,5]=AIC(gy_o, k=log(n))[1]}
if(inherits(criterion, "numeric")==1){ 
G = function(x){I(tan(PI_ext*((x)-0.5)))} 
YY = G(Y)
HHH=cbind(HH,YY)
form3 = as.formula(paste("YY~",paste0(formula)[3]))
gy_o <- quantreg::rq(form3, tau=tau, data=HHH)
fe0 = quantreg::rq(YY~1, tau=tau, data=HHH)
R2 = 1 - gy_o$rho/fe0$rho
A=matrix(0,1,5)
colnames(A) = c("PI","Pseudo-R2","Adjust. pseudo-R2","AIC","BIC")
rownames(A) = c("Values")
A[1,1]=PI_ext
A[1,2]=R2
A[1,3]=1-(1-R2)*(n-1)/(n-dim(HH)[2])
A[1,4]=AIC(gy_o)[1]
A[1,5]=AIC(gy_o, k=log(n))[1]}
betas=as.matrix(gy_o$coef)
Xhat = apply(gy_o$x,2,mean)
if(criterion=="R2"){
F_E2 = (PI_o*(1+(sum(Xhat*betas))^2))^(-1)}
if(criterion=="aic"){
F_E2 = (PI_o2*(1+(sum(Xhat*betas))^2))^(-1)}
if(criterion=="bic"){
F_E2 = (PI_o3*(1+(sum(Xhat*betas))^2))^(-1)}
if(inherits(criterion, "numeric")==1){
F_E2 = (PI_ext*(1+(sum(Xhat*betas))^2))^(-1)}
E2 = F_E2*betas
Ef=cbind(E2) 
colnames(Ef) = c("Em")
quantreg.plot=summary(quantreg::rq(form3, data = HHH, tau = seq(tau_i, tau_f, by = 0.05)))
if(inherits(criterion, "numeric")==0){gy_o$index=t(AA)}
if(inherits(criterion, "numeric")==1){gy_o$index=t(A)}
gy_o$effects=Ef
gy_o$quantregplot=quantreg.plot
if(inherits(criterion, "numeric")==0){gy_o$pis=Pp}
return(gy_o)}
