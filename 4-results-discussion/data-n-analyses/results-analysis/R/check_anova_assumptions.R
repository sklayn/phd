check_anova_assumptions <- function(variable, factor, data = NULL) {
  ## Checks whether the assumptions for a one-way ANOVA are met. 
  ## 
  ## Arguments: variable - string containing name of variable to be tested 
  ##              (the dependent variable)
  ##            factor - string containing name of factor variable 
  ##              (the independent variable)
  ##            data - (optional) data frame containing the variables
  
  # construct the model formula, of the form y ~ x, then fit the linear model
  model <- formula(paste(variable, factor, sep = "~"))
  lm.fit <- aov(model, data = data)
  
  # normality: Shapiro-Wilk test on the residuals of the linear model fit
  shapiro.out <- shapiro.test(residuals(lm.fit))
  
  # homogeneity of variances: Bartlett test
  bartlett.out <- bartlett.test(model, data = data)
  
  # messages for each outcome
  messages.outcomes <- list(shap.sign = "The residuals are not normally distributed. Consider transforming your data!",
							              shap.insign = "The residuals are normally distributed.",
							              bartlett.sign = "The variances are not homogeneous. Consider transforming your data!", 
							              bartlett.insign = "The variances are homogeneous.")
  
  # check the outcome of each test for significance and print the p-value and the corresponding message.
  if(shapiro.out$p.value < 0.05) {
    print(paste("Shapiro-Wilk test: p =", signif(shapiro.out$p.value, 3), 
                messages.outcomes["shap.sign"], sep = " "), quote = FALSE)

  } else {
    print(paste("Shapiro-Wilk test: p =", signif(shapiro.out$p.value, 3), 
                messages.outcomes["shap.insign"], sep = " "), quote = FALSE)  
  }

   
  if(bartlett.out$p.value < 0.05) {
    print(paste("Bartlett test: p =", signif(bartlett.out$p.value, 3), 
                messages.outcomes["bartlett.sign"], sep = " "), quote = FALSE)    
  
  } else {
    print(paste("Bartlett test: p =", signif(bartlett.out$p.value, 3), 
                messages.outcomes["bartlett.insign"], sep = " "), quote = FALSE)
  }

}