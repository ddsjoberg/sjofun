#' Model-based test for treatment efficacy
#'
#' The model-based sinlge-arm comparison addresses the fundamental question of whether the current treatment provides a clinically significant
#' improvement over prior treatments in the population. The proposed test statistic computes the difference between the
#' observed outcome from the current treatment and the covariate-specific predicted outcome based on a model of the
#' historical data. Thus, the difference between the observed and predicted quantities is attributed to the current
#' treatment.
#'
#' Heller, Glenn, Michael W. Kattan, and Howard I. Scher. "Improving the decision to pursue a phase 3 clinical trial
#' by adjusting for patient-specific factors in evaluating phase 2 treatment efficacy data." Medical Decision Making
#' 27.4 (2007): 380-386.
#'
#' @author Daniel D Sjoberg \email{sjobergd@@mskcc.org}
#'
#' @param data a data frame containing the outcome and the outcome predictions.
#' @param model glm object of the predictive model estimated on historical cohort.  If specified, the outcome, covars, cov, and coef objects will be extracted from object.
#' @param outcome the outcome, or response variable name. Must be a variable contained within the data frame specified in data=.
#' @param covars vector of covariate/predictor variable(s) names. Must be a variable(s) contained within the data frame specified in data=.  If model includes an intercept, user must include the column of ones in the covars vector
#' @param cov Variance-covaraince matrix of the beta coefficients from predictive model.
#' @param coef Vector of the beta coefficients from predictive model.  If model includes an intercept, a vector of ones must appear in data.
#' @param type Type of predictive model used.  logistic is currently the only valid input.
#' @param output.details Save additional information.  Default is FALSE.
#'
#' @return Returns a list of results from analysis.
#'
#' @examples
#' set.seed(23432)
#' #simulating historic dataset and creating prediction model.
#' marker=rnorm(500, sd = 2)
#' respond=runif(500)<plogis(marker)
#' historic.data=data.frame(respond,marker)
#' model.fit=glm(data=historic.data, formula = respond ~ marker, family = binomial(logit))
#'
#' #simulating new data, with higher response rate
#' new.data = marker=rnorm(50, sd = 2)
#' respond=runif(50)<plogis(marker + 1)
#' new.data=data.frame(respond,marker)
#'
#' #comparing outcomes in new data to those predicted in historic data
#' # z-statistic = 2.412611 indicates signficant difference
#' model.diff(data = new.data, model = model.fit)
#'
#' #comparing model based difference with binomial test
#' #p-value of 0.3222 indicates we fail to reject null hypothesis
#' binom.test(x=sum(new.data$respond), n=nrow(new.data), p = 0.5, alternative = c("two.sided"))
#'
#' @export
#'
model.diff=function(data, outcome, covars, model, cov, coef, type="logistic", output.details = FALSE) {
  #initializing empty list to store results
  results=list()

  #checking to assess whether prior model was provided.
  if(missing(model) & (missing(cov) | missing(coef)))  {
    stop("Must provide model object (model), or both covariance matrix (cov) and coefficient vector (coef).")
  } else if(!missing(model) & (!missing(cov) | !missing(coef)) ) {
    stop("Provide model object (model), or covariance matrix (cov) and coefficient vector (coef).  But not both.")
  }


  #extracting covariance and coefficient if model object provided, calculating linear predictor
  if(!missing(model)){
    #extracting linear predictor
    pr.xb=predict(model, newdata = data)
    #extracting variance-covariance matrix
    cov=vcov(model)
    #extracting vector of coefficients
    coef=model$coefficients
    #extracting covaraite names
    covars=names(model$coefficients)
    #if model includes intercept, adding vector of ones to data
    if(covars[1]=="(Intercept)") data["(Intercept)"]=1
    #extracting outcome name
    outcome=all.vars(attributes(model$terms)$variables)[1]
  }

  #subsetting on complete data
  data=data[complete.cases(data[,c(outcome, covars)]),]

  #if the model object not specified, using use supplied info to calculate the XB
  if(missing(model)) pr.xb=data[,covars] %*% coef


  #calculating difference in reposnse and predicted outcome
  pr.outcome=exp(pr.xb)/(1+exp(pr.xb))
  delta=data[,outcome] - pr.outcome

  #saving observation level data
  if(output.details==TRUE){
    results[[1]]=cbind(data[,c(outcome,covars)],pr.outcome,pr.xb,delta)
    names(results[[1]])[1:length(c(outcome,covars))]=c(outcome,covars)
    names(results)=c("data")
  }
  #calculating variance of Y- pr(y)
  V1=mean(delta^2)
  V2=cov
  d=colMeans(apply(data[,covars], 2, function(x) x*pr.outcome*(1-pr.outcome)))
  V=as.numeric(nrow(data)*(V1 + t(d) %*% V2 %*% d))
  var.list=list(V, V1, V2, d)
  names(var.list)=c("var.tot","V1","V2","d")
  results=c(results,var.list)

  #calculating test statistic
  S.star=as.vector(sum(delta)/sqrt(V))
  results=c(results,S.star)
  names(results)[length(results)]="z"

  #returning final result
  return(results)
}
