# AMS394 PROJECT by Xiyi Lin

partial.paired.test <- function(x, ...) UseMethod("partial.paired.test")
partial.paired.test <- function(x, y, alternative = c("two.sided", "less", "greater"),
                        mu = 0, conf.level = 0.95, method = c('modified','corrected'), conf.int = TRUE,
                        ...) {
  
  # check if any argument is missing
  if (is.null(x)) 
    stop("Error: input x is missing for performing partial paired test")
  
  if (is.null(y)) 
    stop("Error: input y is missing for performing partial paired test")
  
  if(!missing(conf.level) &&
     (length(conf.level) != 1 || !is.finite(conf.level) ||
      conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  
  if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  
  # sample size check
  if(length(x) < 3) 
    stop("not enough 'x' observations")
  
  if(length(y) < 3) 
    stop("not enough 'y' observations")
  
  # specifying the method to use: by default, modified t-test
  ifelse(is.null(method)=='TRUE',met <- "modified", met <- match.arg(method))
  
  alpha <- 1 - conf.level
  alt <- match.arg(alternative)
  
  # create a data frame and empty vectors
  DF <- data.frame (x,y)
  N <- length(x)
  n1_x <- c()
  n1_y <- c()
  n2_x <- c()
  n2_y <- c()
  n3_x <- c()
  n3_y <- c()
  n4_x <- c()
  n4_y <- c()
  
  # sorting data into four sets based on NAs
  for (i in 1:N) {
    if(is.na(DF[i,][1]) == 'FALSE' && is.na(DF[i,][2]) == 'FALSE')
    {n1_x <- c(n1_x, DF[i,][[1]]) 
    n1_y <- c(n1_y, DF[i,][[2]]) }
    if(is.na(DF[i,][1]) == 'FALSE' && is.na(DF[i,][2]) == 'TRUE')
    {n2_x <- c(n2_x, DF[i,][[1]]) 
    n2_y <- c(n2_y, DF[i,][[2]])}
    if(is.na(DF[i,][1]) == 'TRUE' && is.na(DF[i,][2]) == 'FALSE')
    {n3_x <- c(n3_x, DF[i,][[1]]) 
    n3_y <- c(n3_y, DF[i,][[2]])}
    if(is.na(DF[i,][1]) == 'TRUE' && is.na(DF[i,][2]) == 'TRUE')
    {n4_x <- c(n4_x, DF[i,][[1]]) 
    n4_y <- c(n4_y, DF[i,][[2]])}
  }
  
  #n1 <- matrix(data = c(n1_x,n1_y), byrow = TRUE, nrow = 2) # completely matched
  #n2 <- matrix(data = c(n2_x,n2_y), byrow = TRUE, nrow = 2) # sample 2 missing
  #n3 <- matrix(data = c(n3_x,n3_y), byrow = TRUE, nrow = 2) # sample 1 missing
  #n4 <- matrix(data = c(n4_x,n4_y), byrow = TRUE, nrow = 2) # all NAs
  
  # data set size ni
  n1size <- length(n1_x)
  n2size <- length(n2_x)
  n3size <- length(n3_x)
  n4size <- length(n4_x)
  
  
  # non-normal assumption
  if(shapiro.test(x)$p.value < (1-conf.level) | shapiro.test(y)$p.value < (1-conf.level)){
    cat("Warning in assumption of normality: at least one sample fails the test at the significance level, use Wilcoxon test")
    # no degree of freedom for Wilcoxon test
    
    # 9 cases
    #1. n1>0, n2>0, n3>0, n4=0: 
    if(length(n1)>0 && length(n2)>0 && length(n3)>0)
    {
      # Two Sample paired Wilcoxon test: use n1 set only
      two_wtest <- wilcox.test(n1[1,], n1[2,], paired = TRUE, mu, alternative = alt, conf.int = TRUE, conf.level)
      Ws <- two_wtest$statistic
      cint <- two_wtest$conf.int
      pval <- two_wtest$p.value
    }
    
    #2. n1>0, n2=0, n3>0, n4=0: warning message; paired t-test/two independent samples
    if(length(n1)>0 && length(n2)==0 && length(n3)>0)
    {
      # Two Sample paired Wilcoxon test: use n1 set only
      two_wtest <- wilcox.test(n1[1,], n1[2,], paired = TRUE, mu, alternative = alt, conf.int = TRUE, conf.level)
      Ws <- two_wtest$statistic
      cint <- two_wtest$conf.int
      pval <- two_wtest$p.value
    }
    
    #3. n1>0, n2>0, n3=0, n4=0: warning message; paired t-test/two independent samples
    if(length(n1)>0 && length(n2)>0 && length(n3)==0)
    { 
      # Two Sample paired Wilcoxon test: use n1 set only
      two_wtest <- wilcox.test(n1[1,], n1[2,], paired = TRUE, mu, alternative = alt, conf.int = TRUE, conf.level)
      Ws <- two_wtest$statistic
      cint <- two_wtest$conf.int
      pval <- two_ttest$p.value
    }
    
    #4. n1=0, n2>0, n3>0, n4=0: two independent samples
    if(length(n1)==0 && length(n2)>0 && length(n3)>0)
    { 
      # Two Sample non-paired Wilcoxon test
      two_wtest <- wilcox.test(n2[1,], n3[2,], paired = FALSE, mu, alternative = alt, conf.int = TRUE, conf.level)
      Ws <- two_wtest$statistic
      cint <- two_wtest$conf.int
      pval <- two_wtest$p.value
    }
    
    #5. n1>0, n2=0, n3=0, n4=0: completely paired Wilcoxon test
    if(length(n1)>0 && length(n2)==0 && length(n3)==0)
    {
      # Two Sample paired Wilcoxon test: use n1 set only
      two_wtest <- wilcox.test(n1[1,], n1[2,], paired = TRUE, mu, alternative = alt, conf.int = TRUE, conf.level)
      Ws <- two_wtest$statistic
      cint <- two_wtest$conf.int
      pval <- two_ttest$p.value
    }
    
    #6. n1=0, n2>0, n3=0, n4=0: warning; one sample Wilcoxon test
    if(length(n1)==0 && length(n2)>0 && length(n3)==0)
    {
      cat("Warning: may have some inaccuracy due to using one sample Wilcoxon test for partially paired data")
      # One Sample Wilcoxon test: use n2 set only
      cat("Exact confidence interval for the median in the one-sample case.")
      one_wtest <- wilcox.test(n2[1,], paired = FALSE, mu, alternative = alt, conf.int = TRUE, conf.level)
      Ws <- one_wtest$statistic
      cint <- one_wtest$conf.int
      pval <- one_ttest$p.value
    }
    
    #7. n1=0, n2=0, n3>0, n4=0: warning; one sample t-test
    if(length(n1)==0 && length(n2)==0 && length(n3)>0)
    {
      cat("Warning: may have some inaccuracy due to using one sample t-test for partially paired data")
      # One Sample Wilcoxon test: use n3 set only
      cat("Exact confidence interval for the median in the one-sample case.")
      one_wtest <- wilcox.test(n3[2,], paired = FALSE, mu, alternative = alt, conf.int = TRUE, conf.level)
      Ws <- one_wtest$statistic
      cint <- one_wtest$conf.int
      pval <- one_ttest$p.value
    }
    
    #8. n1=0, n2=0, n3=0, n4=0: Error
    if(length(n1)==0 && length(n2)==0 && length(n3)==0)
      stop("Error: cannot perform Wilcoxon test since all data are missing")
    
    #9. n1=0, n2=0, n3=0, n4!=0: Error
    if(length(n1)==0 && length(n2)==0 && length(n3)==0 && length(n4) != 0)
      stop("Error: cannot perform Wilcoxon test based on all missing values")
    
    # output layout
    
      names(Ts) <- "W"
      names(mu) <- "difference in means"
      attr(cint, "conf.level") <- conf.level
      rval <- list(statistic=Ws, p.value=pval, 
                   conf.int=cint, null.value = mu, alternative=alt,
                   method="Partially paired Wilcoxon rank test")
      class(rval) <- "htest"
      rval
      
    }
  
  # normal assumption
  if((shapiro.test(x)$p.value > (1-conf.level) && shapiro.test(y)$p.value > (1-conf.level)) | length(x)>30) {
  # Warning: Using Central Limit Theorem to approximately assume that the population is normal based on large sample sizes
    
    # 9 cases
  # perform var.test when two independent samples have different sizes 
  # except when n1,n2,n3 all=0, we automatically discard n4 from our analysis
  #1. n1>0, n2>0, n3>0, n4=0: modified-t or corrected-z
  if(length(n1)>0 && length(n2)>0 && length(n3)>0)
  {
    if (met == "modified") {    ###  modified t-test 
    mean_D <- mean(n1[1,]-n1[2,])
    std_D <- sd(n1[1,]-n1[2,])
    mean_T <- mean(n2[1,])
    std_T <- sd(n2[1,])
    mean_N <- mean(n3[2,])
    std_N <- sd(n3[2,])
    n_H <- 2 * n2size * n3size / (n2size + n3size)
    modified_t <- (n1size*mean_D+n_H*(mean_T-mean_N)-mu)/sqrt(n1size*std_D^2+n_H^2*(std_T^2/n2size+std_N^2/n3size))
    Ts <- modified_t
    # confidence interval
    dfd <- 3 # degree of freedom = 3
    cint <- switch(alt,
                   "less" = c(-Inf, modified_t + qt(conf.level, dfd)),
                   "greater" = c(modified_t - qt(conf.level, dfd), Inf),
                   "two.sided" = {
                     ci <- qt(1 - alpha/2, dfd)
                     modified_t + c(-ci, ci) }
                   )
    se <- sqrt(n1size*std_D^2+n_H^2*(std_T^2/n2size+std_N^2/n3size))
    cint <- mu + cint * se
    # p-value    
    pval <- switch(alt,
                   "less" = pt(modified_t, dfd),
                   "greater" = pt(modified_t, dfd, lower.tail=FALSE),
                   "two.sided" = 2 * pt(-abs(modified_t), dfd)
                   )
    }
    
    if (met == "corrected") {   ### corrected z-test
    mean_Tstar <- mean(c(n1[1,],n2[1,]))
    std_Tstar <- sd(c(n1[1,],n2[1,]))
    mean_Nstar <- mean(c(n1[2,],n3[2,]))
    std_Nstar <- sd(c(n1[2,],n3[2,]))
    cov_TN <- cov(n1[1,],n1[2,])
    corrected_z <- (mean_Tstar-mean_Nstar-mu)/sqrt(std_Tstar^2/(n1size+n2size)+std_Nstar^2/(n1size+n3size)-2*n1size*cov_TN/((n1size+n2size)*(n1size+n3size)))
    Ts <- corrected_z
    dfd <- 0 # no degree of freedom for z
    stderr <- sqrt(std_Tstar^2/(n1size+n2size)+std_Nstar^2/(n1size+n3size)-2*n1size*cov_TN/((n1size+n2size)*(n1size+n3size)))
    # confidence interval
    cint <- switch(alt,
                   "less" = c(-Inf, corrected_z * stderr + qnorm(conf.level) * stderr),
                   "greater" = c(corrected_z * stderr - qnorm(conf.level) * stderr, Inf),
                   "two.sided" = {
                     c(corrected_z * stderr - qnorm((1 - alpha/2)) * stderr, corrected_z *
                         stderr + qnorm((1 - alpha/2)) * stderr) } 
                   )
    cint <- cint + mu
    # p-value
    pval <- switch(alt,
                   "less" = pnorm(corrected_z),
                   "greater" = 1 - pnorm(corrected_z),
                   "two.sided" = 2 * pnorm( - abs(corrected_z))
                   )
    }
  }
  
  #2. n1>0, n2=0, n3>0, n4=0: warning message; paired t-test/two independent samples
  if(length(n1)>0 && length(n2)==0 && length(n3)>0)
  {
    # variance test
    ifelse(var.test(n1[1,],c(n1[2,],n3[2,]))$p.value > 1-conf.level,var <- "TRUE",var <- "FALSE") 
    
    # Two Sample t-test 
    two_ttest <- t.test(n1[1,],c(n1[2,],n3[2,]), alternative = alt, mu, paired = FALSE, var.equal = var, conf.level) 
    cat("Warning: may have some inaccuracy due to using Two Sample t-test for partially paired data")
    Ts <- two_ttest$statistic
    dfd <- two_ttest$parameter
    cint <- two_ttest$conf.int
    pval <- two_ttest$p.value
    
    # Alternatively, paired t-test
    paired_ttest <- t.test(n1[1,],n1[2,], alternative = alt, mu, paired = TRUE, var.equal = var, conf.level)
    cat("Warning: may have some inaccuracy due to using Completely Paired t-test for partially paired data")
    Ts1 <- paired_ttest$statistic
    dfd1 <- paired_ttest$parameter
    cint1 <- paired_ttest$conf.int
    pval1 <- paired_ttest$p.value
  }
  
  #3. n1>0, n2>0, n3=0, n4=0: warning message; paired t-test/two independent samples
  if(length(n1)>0 && length(n2)>0 && length(n3)==0)
  { 
    # variance test
    ifelse(var.test(n1[2,],c(n1[1,],n2[1,]))$p.value > 1-conf.level,var <- "TRUE",var <- "FALSE") 
    
    # Two Sample t-test 
    two_ttest <- t.test(c(n1[1,],n2[1,]), n1[2,], alternative = alt, mu, paired = FALSE, 
                        var.equal = var, conf.level) 
    cat("Warning: may have some inaccuracy due to using Two Sample t-test for partially paired data")
    Ts <- two_ttest$statistic
    dfd <- two_ttest$parameter
    cint <- two_ttest$conf.int
    pval <- two_ttest$p.value
    
    # Alternatively, paired t-test
    paired_ttest <- t.test(n1[1,],n1[2,], alternative = alt, mu, paired = TRUE, 
                           var.equal = var, conf.level)
    cat("Warning: may have some inaccuracy due to using Completely Paired t-test for partially paired data")
    Ts1 <- paired_ttest$statistic
    dfd1 <- paired_ttest$parameter
    cint1 <- paired_ttest$conf.int
    pval1 <- paired_ttest$p.value
  }
  
  #4. n1=0, n2>0, n3>0, n4=0: two independent samples
  if(length(n1)==0 && length(n2)>0 && length(n3)>0)
  { 
    # variance test
    ifelse(var.test(n2[1,],n3[2,])$p.value > 1-conf.level,var <- "TRUE",var <- "FALSE") 
    
    # subtract sample 1(x) from n2 and sample 2(y) from n3 to form an independent sample
    two_ttest <- t.test(n2[1,], n3[2,], alternative = alt, mu, paired = FALSE, 
                        var.equal = var, conf.level) # Two Sample t-test
    Ts <- two_ttest$statistic
    dfd <- two_ttest$parameter
    cint <- two_ttest$conf.int
    pval <- two_ttest$p.value
  }
  
  #5. n1>0, n2=0, n3=0, n4=0: completely paired t test
  if(length(n1)>0 && length(n2)==0 && length(n3)==0)
  {
    paired_ttest <- t.test(n1[1,],n1[2,], alternative = alt, mu, paired = TRUE, 
                           var.equal = FALSE, conf.level)
    Ts <- paired_ttest$statistic
    dfd <- paired_ttest$parameter
    cint <- paired_ttest$conf.int
    pval <- paired_ttest$p.value
  }
  
  #6. n1=0, n2>0, n3=0, n4=0: warning; one sample t-test
  if(length(n1)==0 && length(n2)>0 && length(n3)==0)
  {
    cat("Warning: may have some inaccuracy due to using one sample t-test for partially paired data")
    one_ttest <- t.test(n2[1,], alternative = alt, mu, paired = FALSE, 
                        var.equal = FALSE, conf.level=0.95) # One Sample t-test
    Ts <- one_ttest$statistic
    dfd <- one_ttest$parameter
    cint <- one_ttest$conf.int
    pval <- one_ttest$p.value
  }
  
  #7. n1=0, n2=0, n3>0, n4=0: warning; one sample t-test
  if(length(n1)==0 && length(n2)==0 && length(n3)>0)
  {
    cat("Warning: may have some inaccuracy due to using one sample t-test for partially paired data")
    one_ttest <- t.test(n3[2,], alternative = alt, mu, paired = FALSE, 
                        var.equal = FALSE, conf.level=0.95) # One Sample t-test
    Ts <- one_ttest$statistic
    dfd <- one_ttest$parameter
    cint <- one_ttest$conf.int
    pval <- one_ttest$p.value
  }
  
  #8. n1=0, n2=0, n3=0, n4=0: Error
  if(length(n1)==0 && length(n2)==0 && length(n3)==0)
    stop("Error: cannot perform t-test since all data are missing")
  
  #9. n1=0, n2=0, n3=0, n4!=0: Error
  if(length(n1)==0 && length(n2)==0 && length(n3)==0 && length(n4) != 0)
    stop("Error: cannot perform t-test based on all missing values")
  
  # output layout
  
  if((length(n1)>0 && length(n2)==0 && length(n3)>0) || (length(n1)>0 && length(n2)>0 && length(n3)==0)) {
    names(c(Ts,Ts1)) <- "t"
    names(c(dfd,dfd1)) <- "df"
    names(mu) <- "difference in means"
    attr(c(cint,cint1), "conf.level") <- conf.level
    rval <- list(statistic=c('two sample t-test',Ts,'paired t-test',Ts1), parameter = c(dfd,dfd1), p.value=c(pval,pval1), 
                 conf.int=c(cint,cint1), null.value = mu, alternative=alt,
                 method="Partially paired t-test")
    class(rval) <- "htest"
    rval
  }
  else {
    names(Ts) <- "t"
    names(dfd) <- "df"
    names(mu) <- "difference in means"
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic=Ts, parameter = dfd, p.value=pval, 
                 conf.int=cint, null.value = mu, alternative=alt,
                 method="Partially paired t-test")
    class(rval) <- "htest"
    rval
  }
    }
  return(rval)
}
