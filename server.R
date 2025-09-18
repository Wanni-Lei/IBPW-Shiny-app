## functions
library(tidyverse)
library(shiny)
library(devtools)
library(clinfun)
library(Rcpp)
library(RcppArmadillo)
library(shinycssloaders)
library(readxl)

## posterior probability by integral method
# url <- "https://github.com/Wanni-Lei/data/raw/refs/heads/main/bayesian_matrix_integral.xlsx"
# destfile <- "bayesian_matrix_integral.xlsx"
# 
# #Download the file
# download.file(url, destfile, mode = "wb")
# library(readxl)
# posterior_matrix <- read_excel(destfile)

## posterior probability by 10K simulation
url <- "https://github.com/Wanni-Lei/data/raw/refs/heads/main/bayesian_matrix_10kSimulation.xlsx"
destfile <- "bayesian_matrix_10kSimulation.xlsx"

#Download the file
download.file(url, destfile, mode = "wb")
library(readxl)
posterior_matrix <- read_excel(destfile)

## use rcpp generate each scenarios
library(Rcpp)

# Define the Rcpp function 
cppFunction('
  DataFrame generateResultsRcpp(int nmax) {
    std::vector<int> n_values;
    std::vector<int> n1_values;
    std::vector<int> r_values;
    std::vector<int> r1_values;
    int n1max;

    // Fix n
    for (int n = 6; n <= nmax; ++n) {
      // Fix n1
      if (n < 10) {
        n1max = n - 2;
      } else {
       n1max = n - 3;
      }

      for (int n1 = 3; n1 <= n1max; ++n1) {

        // Fix r1
        for (int r1 = (n1 - 1); r1 >= 0; --r1) {
          // Fix r
          for (int r = (r1 + (n-n1) -1); r >= (r1 + 1); --r) {
            n_values.push_back(n);
            n1_values.push_back(n1);
            r_values.push_back(r);
            r1_values.push_back(r1);
          }
        }
      }
    }

    // Create a DataFrame from the vectors
    DataFrame df = DataFrame::create(_["n"] = n_values, _["n1"] = n1_values, _["r"] = r_values, _["r1"] = r1_values);

    return df;
  }
')

# Call the Rcpp function, the maximum sample size for single arm is set to 100
results_rcpp <- generateResultsRcpp(100)
rownames(results_rcpp) <- NULL

library(tidyverse)
results_rcpp <- as.data.frame(results_rcpp)
scenarios <- results_rcpp %>% 
  arrange(n,n1,r1)




optimal_minimax <- function(pa1,pb1, alpha_input, beta_input, pa0, pb0, delta, scenarios, posterior_matrix){
  
  
  if( pa1>pb1 | pa0>pb0){ 
    
    stop("Error: pa1 must be less than pb1, and pa0 must be less or equal than pb0") 
    
  } 
  
  
  
  simon <-ph2simon(pu= max(pa1, pa0, pb0), pa= pb1, alpha_input, beta_input) 
  df <- as.data.frame(simon$out) 
  df_optimal<- df[which.min(df$`EN(p0)`), ] 
  df_minimax <- df[which.min(df$n), ]
  n_simon <- as.numeric(df_optimal[4]) ## n of Simon's optimal design 
  n_simon_minimax <- as.numeric(df_minimax[4]) ## n of Simon's minimax design
  n1_simon <- max(10, n_simon) ## restrict upper bound of n 
  
  ## chose maximum sample size for single arm 
  
  if((pb1-pa1)>= 0.3){ 
    choice<- scenarios %>% 
      filter(0.7*n_simon <= n & n <= 1.4*n1_simon ## restrict n 
             #& n1/(n-n1) >= 0.55 & n1/(n-n1) <= 1 ## restrict n1/n2 
             & pb0 >= (r1-2)/n1 & r1>= n1*(pa0-0.3) ## restrict r1 
             & pa1 <= (r-r1)/(n-n1) & (r-r1)/(n-n1) <= (pb1+0.4) ) ## restrict r2 
    
  }else{ 
    choice<- scenarios %>% 
      filter(0.7*n_simon <= n & n <= 1.2*n1_simon ## restrict n 
             & n1/(n-n1) >= 0.65 & n1/(n-n1) <= 5 ## restrict n1/n2 
             & pb0 >= (r1-3)/n1 & r1>= n1*(pa0-0.1) ## restrict r1 
             & pa1 <= (r-r1)/(n-n1) & (r-r1)/(n-n1) <= (pb1+0.1)) 
    
  } 
  
  
  
  
  ## apply the simonBayesian function to the dataframe 
  
  results <- choice %>% 
    mutate(EN_p0 = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4], 
                                                               pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0, delta=delta, posterior_matrix=posterior_matrix )$ess), 
           power = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4], 
                                                               pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0, delta=delta,posterior_matrix=posterior_matrix )$power), 
           alpha = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4], 
                                                               pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0,delta=delta, posterior_matrix=posterior_matrix )$alpha)) 
  
  
  
  
  
  valid_scenario <- results %>% 
    filter(EN_p0 <= 2*n) %>% 
    filter((1 - power) <= beta_input & alpha <= alpha_input) 
  
  
  if(nrow(valid_scenario)==0){
    n_index <- n_simon_minimax -10
    while(nrow(valid_scenario) == 0 && n_index < 100){
      
      choice<- scenarios %>%
        filter(n_index <= n & n <= (n_index +10) ## restrict n
               #& n1/(n-n1) >= 0.55 & n1/(n-n1) <= 1 ## restrict n1/n2
               & pb0 >= (r1-3)/n1 & r1>= n1*(pa0-0.3) ## restrict r1 
               & pa1 <= (r-r1)/(n-n1) & (r-r1)/(n-n1) <= (pb1+0.4))  ## restrict r2
      
      
      ## apply the simonBayesian function to the dataframe
      # results <- choice %>%
      #   mutate(EN_p0 = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4],
      #                                                              pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0, delta=delta, posterior_matrix=posterior_matrix )$ess),
      #          power = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4],
      #                                                              pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0, delta=delta,posterior_matrix=posterior_matrix )$power),
      #          alpha = apply(choice, 1, function(x) simon_bayesian(x[1], x[2], x[3], x[4],
      #                                                              pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0,delta=delta, posterior_matrix=posterior_matrix )$alpha))
      # 
      
      results0 <- apply(choice, 1, function(x) {
        y <- simon_bayesian(x[1], x[2], x[3], x[4], pa1=pa1, pb1=pb1, pa0=pa0, pb0=pb0,
                            delta=delta, posterior_matrix=posterior_matrix)
        return(c(EN_p0=y$ess, power=y$power, alpha=y$alpha))
      })
      
      results <- data.frame(choice, t(as.matrix(results0)))
         
      valid_scenario <- results %>%
        filter(EN_p0 <= 2*n) %>%
        filter((1 - power) <= beta_input & alpha <= alpha_input)
      
      n_index <- n_index +10 ##increase by 10 every time
    }
    optimal <- valid_scenario %>%
      filter(EN_p0 == min(EN_p0))%>% ## filter minimum EN_p0
      ## filter power closest to power_input
      filter(abs(power - (1-beta_input)) == min(abs(power - (1-beta_input))))%>%
      filter(abs(alpha - alpha_input) == min(abs(alpha - alpha_input))) ## alpha power closest to alpha_input
    
    
    
    minimax <- valid_scenario %>%
      filter(n == min(n))%>% ## filter minimum EN_p0
      ## filter power+alpha has minimum deviation
      filter((abs(power - (1-beta_input))+abs(alpha - alpha_input)) ==
               min(abs(power - (1-beta_input))+abs(alpha - alpha_input)))%>%
      filter(EN_p0==min(EN_p0)) ## filter minimum expected sample size
  
  }  else{ ## if feasible solution is found
    
    optimal <- valid_scenario %>%
      filter(EN_p0 == min(EN_p0))%>% ## filter minimum EN_p0
      ## filter power closest to power_input
      filter(abs(power - (1-beta_input)) == min(abs(power - (1-beta_input))))%>%
      filter(abs(alpha - alpha_input) == min(abs(alpha - alpha_input))) ## alpha power closest to alpha_input
    
    
    
    minimax <- valid_scenario %>%
      filter(n == min(n))%>% ## filter minimum EN_p0
      ## filter power+alpha has minimum deviation
      filter((abs(power - (1-beta_input))+abs(alpha - alpha_input)) ==
               min(abs(power - (1-beta_input))+abs(alpha - alpha_input)))%>%
      filter(EN_p0==min(EN_p0)) ## filter minimum expected sample size
    
  }
    
    
    df <- as.data.frame(rbind(optimal,minimax))
    n_optimal <- nrow(optimal)
    n_minimax <- nrow(minimax)
    
    # Generate unique row names
    optimal_names <- paste("optimal", seq_len(n_optimal), sep = "_")
    minimax_names <- paste("minimax", seq_len(n_minimax), sep = "_")
    
    row.names(df) <- c(optimal_names, minimax_names)
  
  
  if(nrow(valid_scenario)==0){
    stop("No feasible solution found")
  }
  
  
  return(df)
  
}



simon_bayesian <- function(n, n1, r, r1, pa1, pb1, pa0, pb0, delta, posterior_matrix){

  if( pa1 > pb1 | pa0 > pb0){
    stop("Error: pa1 must be less than pb1, and pa0 must be less or equal than pb0")
  }

  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1),
                    x2=rep(0:n2, n1+1))
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1",
                            x1>r1&(x1+x2)<=r ~ "F2",
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1),
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)),
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )

  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()

  row.names(x12.smry) <-x12.smry$result


  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1)

  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)


  ## when both 2 doses pase stage 1
  nr.x12 <- dim(x12)[1]
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))


  ab.pf <- ab %>%
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>%
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))

  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>%
    filter(bH.aPbP==1)  %>%
    filter((a1+a2)<(b1+b2))

  ## Indicator matrix:
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)

  ## Sort the data frame by ascending order of a_total and b_total,
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id())


  ## load posterior matrix
  ## changed simulation in p.rbH to 10000 times
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)

  ## add the "posterior" values from the "posterior_matrix_filter" to
  ## the "indicator" based on matching values of
  ## "a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];

  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];

  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }

  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,
                           _["group"] = indicator["group"]);
}

')
  # indicator <- indicator %>%
  #   arrange(a_total, b_total) %>%
  #   group_by(a_total, b_total) %>%
  #   mutate(sample_size = n)  %>%
  #   mutate(group = cur_group_id())

  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)

  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)


  ## add joint probability
  indicator_p <- indicator %>%
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))

  ## Pr(B>A) > delta
  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior>delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)


  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])



  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]

  ## when response rate of 2 doses are not equal
  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0

  #my_list <- round(list(n = n,n1 = n1, r = r, r1 = r1, ess=ess, power=power,alpha=alpha),3)
  vector <- round(c(ess, power, alpha), 3)
  vector_list <- as.list(vector)
  names(vector_list) <- c("ess", "power", "alpha")
  #vector_list
  return(vector_list)
}


## decision rule table
#designs <- optimal_minimax(pa1=0.05, pb1=0.4, alpha_input = 0.1, beta_input = 0.2, pa0 = 0.05, pb0 = 0.05, delta = 0.8, scenarios = scenarios, posterior_matrix = posterior_matrix)

decision_table <- function(label, designs, pa1, pb1, pa0, pb0, delta, posterior_matrix){
  
  
  if(label =="Optimal"){
    n <- designs[1, "n"]
    n1 <- designs[1, "n1"]
    r <- designs[1, "r"]
    r1 <- designs[1, "r1"] 
  }else{
    n <- designs[2, "n"]
    n1 <- designs[2, "n1"]
    r <- designs[2, "r"]
    r1 <- designs[2, "r1"] 
  }
  
  if( pa1 > pb1 | pa0 > pb0){
    stop("Error: pa1 must be less than pb1, and pa0 must be less or equal than pb0")
  }
  
  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1),
                    x2=rep(0:n2, n1+1))
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1",
                            x1>r1&(x1+x2)<=r ~ "F2",
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1),
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)),
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )
  
  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()
  
  row.names(x12.smry) <-x12.smry$result
  
  
  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1)
  
  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)
  
  
  ## when both 2 doses pase stage 1
  nr.x12 <- dim(x12)[1]
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))
  
  
  ab.pf <- ab %>%
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>%
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))
  
  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>%
    filter(bH.aPbP==1)  %>%
    filter((a1+a2)<(b1+b2))
  
  ## Indicator matrix:
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)
  
  ## Sort the data frame by ascending order of a_total and b_total,
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id())
  
  
  ## load posterior matrix
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  ## add the "posterior" values from the "posterior_matrix_filter" to
  ## the "indicator" based on matching values of
  ## "a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];

  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];

  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }

  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,
                           _["group"] = indicator["group"]);
}

')
  # indicator <- indicator %>%
  #   arrange(a_total, b_total) %>%
  #   group_by(a_total, b_total) %>%
  #   mutate(sample_size = n)  %>%
  #   mutate(group = cur_group_id())
  
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)
  
  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)
  
  
  ## add joint probability
  indicator_p <- indicator %>%
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))
  
  ## Pr(B>A) > delta
  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior>delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)
  
  
  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])
  
  
  
  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]
  
  ## when response rate of 2 doses are not equal
  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0
  
  #my_list <- round(list(n = n,n1 = n1, r = r, r1 = r1, ess=ess, power=power,alpha=alpha),3)
  vector <- round(c(ess, power, alpha), 3)
  vector_list <- as.list(vector)
  names(vector_list) <- c("ess", "power", "alpha")
  #vector_list
  
  ## this part is for decision table
  
  # 1. Get the posterior probability data for both arms passing
  response_both_pass <- indicator_p %>%
    select(a_total, b_total, posterior) %>%
    distinct() %>%
    filter(a_total > r & b_total > r)  # Both pass efficacy threshold
  
  # 2. Generate all possible response combinations  
  max_responses <- n
  
  all_combinations <- expand.grid(
    a_total = 0:max_responses,
    b_total = 0:max_responses
  ) %>%
    filter(a_total <= n & b_total <= n)
  
  # 3. Apply decision logic considering both symmetric cases
  decision_data <- all_combinations %>%
    mutate(
      a_passes = a_total > r,  # r = 2
      b_passes = b_total > r,  # r = 2
      
      label = case_when(
        # Only A passes
        a_passes & !b_passes ~ "A wins",
        # Only B passes  
        !a_passes & b_passes ~ "B wins",
        # Neither passes
        !a_passes & !b_passes ~ "no winner",
        # Both pass - check posterior probability
        a_passes & b_passes ~ "check_posterior"
      )
    )
  
  # 4. Update cases where both pass with symmetric posterior probability logic
  decision_data <- decision_data %>%
    left_join(
      response_both_pass %>% select(a_total, b_total, posterior),
      by = c("a_total", "b_total")
    ) %>%
    mutate(
      label = case_when(
        # When both pass, check posterior in both directions
        label == "check_posterior" & !is.na(posterior) & posterior > delta ~ "B wins",
        label == "check_posterior" & !is.na(posterior) & posterior < (1-delta) ~ "A wins",
        label == "check_posterior" ~ "no winner",
        TRUE ~ label
      )
    ) %>%
    select(a_total, b_total, label)
  
  # 5. Add the symmetric cases by swapping a_total and b_total for posterior decisions
  # This ensures we consider both directions when both arms pass
  symmetric_posterior <- response_both_pass %>%
    select(a_total, b_total, posterior) %>%
    # Create symmetric version by swapping columns and inverting posterior
    mutate(
      a_temp = b_total,
      b_temp = a_total, 
      posterior_inv = 1 - posterior
    ) %>%
    select(a_total = a_temp, b_total = b_temp, posterior = posterior_inv) %>%
    # Remove duplicates (cases where a_total = b_total)
    filter(a_total != b_total)
  
  # Combine original and symmetric posterior data
  all_posterior_data <- bind_rows(
    response_both_pass %>% select(a_total, b_total, posterior),
    symmetric_posterior
  ) %>%
    distinct()
  
  # 6. Re-apply decision logic with complete symmetric posterior data
  decision_data_complete <- all_combinations %>%
    mutate(
      a_passes = a_total > r,
      b_passes = b_total > r,
      
      label = case_when(
        a_passes & !b_passes ~ "A wins",
        !a_passes & b_passes ~ "B wins", 
        !a_passes & !b_passes ~ "no winner",
        a_passes & b_passes ~ "check_posterior"
      )
    ) %>%
    left_join(
      all_posterior_data,
      by = c("a_total", "b_total")
    ) %>%
    mutate(
      label = case_when(
        label == "check_posterior" & !is.na(posterior) & posterior > delta ~ "B wins",
        label == "check_posterior" & !is.na(posterior) & posterior < (1-delta) ~ "A wins", 
        label == "check_posterior" ~ "no winner",
        TRUE ~ label
      )
    ) %>%
    select(a_total, b_total, label)
  
  # 7. Create the final decision table
  ra_values <- sort(unique(decision_data_complete$a_total))
  
  final_table <- data.frame(Ya = character(), 
                            `Arm A wins` = character(),
                            `Arm B wins` = character(), 
                            `No winner` = character(),
                            stringsAsFactors = FALSE)
  
  for(ra in ra_values) {
    ra_data <- decision_data_complete %>% filter(a_total == ra) %>% arrange(b_total)
    
    # Get Yb values for each outcome
    a_wins_rb <- ra_data %>% filter(label == "A wins") %>% pull(b_total)
    b_wins_rb <- ra_data %>% filter(label == "B wins") %>% pull(b_total)  
    no_winner_rb <- ra_data %>% filter(label == "no winner") %>% pull(b_total)
    
    # Format rule strings
    format_rb_rule <- function(rb_vals) {
      if(length(rb_vals) == 0) return("_")
      
      rb_vals <- sort(unique(rb_vals))
      
      if(length(rb_vals) == 1) return(paste("Yb =", rb_vals))
      
      # Find consecutive ranges
      ranges <- list()
      i <- 1
      while(i <= length(rb_vals)) {
        start_val <- rb_vals[i]
        end_val <- start_val
        
        # Find end of consecutive sequence
        while(i < length(rb_vals) && rb_vals[i+1] == rb_vals[i] + 1) {
          i <- i + 1
          end_val <- rb_vals[i]
        }
        
        if(start_val == end_val) {
          ranges <- append(ranges, as.character(start_val))
        } else {
          ranges <- append(ranges, paste(start_val, end_val, sep = "-"))
        }
        i <- i + 1
      }
      
      # Format final output
      if(length(ranges) == 1) {
        if(grepl("-", ranges[[1]])) {
          parts <- strsplit(ranges[[1]], "-")[[1]]
          return(paste0(parts[1], " ≤ Yb ≤ ", parts[2]))
        } else {
          return(paste("Yb =", ranges[[1]]))
        }
      } else {
        return(paste0("Yb ∈ {", paste(ranges, collapse = ", "), "}"))
      }
    }
    
    final_table <- rbind(final_table, data.frame(
      Ya = paste("Ya =", ra),
      `Arm.A.wins` = format_rb_rule(a_wins_rb),
      `No.winner` = format_rb_rule(no_winner_rb),
      `Arm.B.wins` = format_rb_rule(b_wins_rb),
      #`No.winner` = format_rb_rule(no_winner_rb),
      stringsAsFactors = FALSE
    ))
  }
  
  names(final_table) <- c("Ya", "Arm A wins", "No winner", "Arm B wins")
  
  # Print table
  
 library(knitr)
return(final_table)
  
  
}








## operating characteristics function
oc <- function(n, n1, r, r1, pa1, pb1, pa0, pb0, delta, posterior_matrix){

  if( pa1>pb1 | pa0>pb0){
    stop("Error: pa1 must be less than pb1, and pa0 must be less or equal than pb0")
  }

  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1),
                    x2=rep(0:n2, n1+1))
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1",
                            x1>r1&(x1+x2)<=r ~ "F2",
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1),
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)),
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )

  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()

  row.names(x12.smry) <-x12.smry$result


  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1)

  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)


  ## when both 2 doses pase stage 1
  nr.x12 <- dim(x12)[1]
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))


  ab.pf <- ab %>%
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>%
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))

  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>%
    filter(bH.aPbP==1)  %>%
    filter((a1+a2)<(b1+b2))

  ## Indicator matrix:
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)

  ## Sort the data frame by ascending order of a_total and b_total,
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id())


  ## load posterior matrix
  ## changed simulation in p.rbH to 10000 times
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)

  ## add the "posterior" values from the "posterior_matrix_filter" to
  ## the "indicator" based on matching values of
  ##"a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];

  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];

  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }

  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,
                           _["group"] = indicator["group"]);
}

')
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id())

  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)

  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)



  indicator_p <- indicator %>%
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))


  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior> delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)


  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])



  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]

  ## when response rate of 2 doses are not equal
  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0

  A.ES.H0 <- sum(table0[1,1:3]) ## arm A failing at stage 1
  B.ES.H0 <- colSums(table0)[1] ## arm B failing at stage 1
  AB.ES.H0 <- A.ES.H0 + B.ES.H0 - table0[1,1] ## Arm A/B failing at stage 1
  A.P1.H1 <- 1- sum(table1[1,1:3]) ## arm A passing stage 1
  B.P1.H1 <- 1- colSums(table1)[1] ## arm B passing stage 1
  AB.P1.H1 <- table1[2,2]+table1[2,3]+table1[3,2]+table1[3,3]


  vector <- list("n" = n,"n1" = n1, "r" = r, "r1" = r1, "ess"=ess, "power"=power,"alpha"=alpha, "A.ES.H0"=A.ES.H0, "B.ES.H0"=B.ES.H0,"AB.ES.H0"=AB.ES.H0, "A.P1.H1"=A.P1.H1,"B.P1.H1"=B.P1.H1, "AB.P1.H1"=AB.P1.H1)
  vector <- as.data.frame(vector)
  colnames(vector) <- c("n","n1","r","r1","ess", "power", "alpha","A.ES.H0", "B.ES.H0","AB.ES.H0", "A.P1.H1","B.P1.H1", "AB.P1.H1")
  rownames(vector) <- "Results"

  return(round(vector,3))
}


## return table1
octable1 <- function(n, n1, r, r1, pa1, pb1, pa0, pb0, delta, posterior_matrix){

  if( pa1>pb1 | pa0>pb0){
    stop("Error: pa1 must be less than pb1, and pa0 must be less or equal than pb0")
  }

  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1),
                    x2=rep(0:n2, n1+1))
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1",
                            x1>r1&(x1+x2)<=r ~ "F2",
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1),
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)),
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )

  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()

  row.names(x12.smry) <-x12.smry$result


  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1)

  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)


  #colnames(table1) <- c("B.fail.stage1", "B.fail.stage2", "B.pass")
  #rownames(table1) <- c("A.fail.stage1", "A.fail.stage2", "A.pass")

  ## when both 2 doses pase stage 1
  nr.x12 <- dim(x12)[1]
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))


  ab.pf <- ab %>%
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>%
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))

  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>%
    filter(bH.aPbP==1)  %>%
    filter((a1+a2)<(b1+b2))

  ## Indicator matrix:
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)

  ## Sort the data frame by ascending order of a_total and b_total,
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id())


  ## load posterior matrix
  ## changed simulation in p.rbH to 10000 times
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)

  ## add the "posterior" values from the "posterior_matrix_filter" to
  ## the "indicator" based on matching values of
  ##"a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];

  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];

  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }

  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,
                           _["group"] = indicator["group"]);
}

')
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id())

  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)

  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)



  indicator_p <- indicator %>%
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))


  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior> delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)


  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])



  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]

  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0

  table1<- as.data.frame(table1)
  table1 <- round(table1,3)

  table1[3,3] <- paste0(table1[3,3], " (B claims as winner: ", round(bpower_1,3), ")")
  return(table1)
}




## return table 0
octable0 <- function(n, n1, r, r1, pa1, pb1, pa0, pb0, delta, posterior_matrix){

  if( pa1>pb1 | pa0>pb0){
    stop("Error: pa1 must be less than pb1, and pa0 must be less or equal than pb0")
  }

  n2 <- n-n1
  x12 <- data.frame(x1=rep(0:n1, each=n2+1),
                    x2=rep(0:n2, n1+1))
  x12 <- x12 %>% filter(!(x1<=r1&x2>0)) %>%
    mutate(result=case_when(x1<=r1 ~ "F1",
                            x1>r1&(x1+x2)<=r ~ "F2",
                            x1>r1 & (x1+x2)>r ~ "P"),
           pa1=ifelse(x1<=r1, dbinom(x1, n1, p=pa1),
                      dbinom(x1,n1, p=pa1)*dbinom(x2,n2, p=pa1)),
           pb1=ifelse(x1<=r1, dbinom(x1, n1, p=pb1),
                      dbinom(x1,n1, p=pb1)*dbinom(x2,n2, p=pb1)),
           pa0=ifelse(x1<=r1, dbinom(x1, n1, p=pa0),
                      dbinom(x1,n1, p=pa0)*dbinom(x2,n2, p=pa0)),
           pb0=ifelse(x1<=r1, dbinom(x1, n1, p=pb0),
                      dbinom(x1,n1, p=pb0)*dbinom(x2,n2, p=pb0))
    )

  x12.smry <- x12[,-(1:2)] %>%
    group_by(result) %>%
    summarize(p.pa1=sum(pa1),
              p.pb1=sum(pb1),
              p.pa0=sum(pa0),
              p.pb0=sum(pb0)) %>% data.frame()

  row.names(x12.smry) <-x12.smry$result


  table1 <- outer(x12.smry$p.pa1,x12.smry$p.pb1)

  table0 <- outer(x12.smry$p.pa0,x12.smry$p.pb0)


  #colnames(table0) <- c("B.fail.stage1", "B.fail.stage2", "B.pass")
  #rownames(table0) <- c("A.fail.stage1", "A.fail.stage2", "A.pass")
  ## when both 2 doses pase stage 1
  nr.x12 <- dim(x12)[1]
  ab <- data.frame(a1=rep(x12$x1,nr.x12),
                   a2=rep(x12$x2,nr.x12),
                   b1=rep(x12$x1,each=nr.x12),
                   b2=rep(x12$x2,each=nr.x12))


  ab.pf <- ab %>%
    mutate(a.out=case_when(a1<=r1 ~ "aF1",
                           a1>r1&(a1+a2)<=r ~ "aF2",
                           a1>r1 & (a1+a2)>r ~ "aP"),
           b.out=case_when(b1<=r1 ~ "bF1",
                           b1>r1&(b1+b2)<=r ~ "bF2",
                           b1>r1 & (b1+b2)>r ~ "bP")) %>%
    mutate(bH.aPbP=ifelse(a.out=="aP"& b.out=="bP",1,0))

  ## filter only (a1+a2)<(b1+b2)
  bW.aPbP <- ab.pf %>%
    filter(bH.aPbP==1)  %>%
    filter((a1+a2)<(b1+b2))

  ## Indicator matrix:
  indicator <- data.frame(a1 = bW.aPbP$a1,
                          a2 = bW.aPbP$a2,
                          b1 = bW.aPbP$b1,
                          b2 = bW.aPbP$b2,
                          a_total = bW.aPbP$a1+bW.aPbP$a2,
                          b_total = bW.aPbP$b1+bW.aPbP$b2,
                          posterior = 0)

  ## Sort the data frame by ascending order of a_total and b_total,
  ## create groups for unique a_total, b_total
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id())


  ## load posterior matrix
  ## changed simulation in p.rbH to 10000 times
  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)

  ## add the "posterior" values from the "posterior_matrix_filter" to
  ## the "indicator" based on matching values of
  ##"a_total", "b_total", and "sample_size"
  cppFunction('
  DataFrame add_posterior_values(DataFrame indicator, DataFrame posterior_matrix_filter, int n) {
  int nrows = indicator.nrows();
  NumericVector group = indicator["group"];
  NumericVector a_total = indicator["a_total"];
  NumericVector b_total = indicator["b_total"];
  NumericVector posterior = indicator["posterior"];

  int num_groups = posterior_matrix_filter.nrows();
  NumericVector filter_a_total = posterior_matrix_filter["a_total"];
  NumericVector filter_b_total = posterior_matrix_filter["b_total"];
  NumericVector filter_sample_size = posterior_matrix_filter["sample_size"];
  NumericVector filter_posterior = posterior_matrix_filter["posterior"];

  for (int i = 0; i < num_groups; ++i) {
    double cur_posterior = filter_posterior[i];
    for (int j = 0; j < nrows; ++j) {
      if (group[j] == (i + 1) && a_total[j] == filter_a_total[i] && b_total[j] == filter_b_total[i] && n == filter_sample_size[i]) {
        posterior[j] = cur_posterior;
      }
    }
  }

  return DataFrame::create(_["a1"] = indicator["a1"],
                           _["a2"] = indicator["a2"],
                           _["b1"] = indicator["b1"],
                           _["b2"] = indicator["b2"],
                           _["a_total"] = indicator["a_total"],
                           _["b_total"] = indicator["b_total"],
                           _["posterior"] = posterior,
                           _["group"] = indicator["group"]);
}

')
  indicator <- indicator %>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total) %>%
    mutate(sample_size = n)  %>%
    mutate(group = cur_group_id())

  posterior_matrix_filter <- posterior_matrix %>%
    filter(sample_size==n &a_total>r & b_total>r &a_total<b_total)%>%
    arrange(a_total, b_total) %>%
    group_by(a_total, b_total)

  indicator <- add_posterior_values(indicator, posterior_matrix_filter, n)



  indicator_p <- indicator %>%
    mutate(p.ra1rb1=dbinom(a1, n1, p=pa1)*dbinom(a2, n2, p=pa1)*
             dbinom(b1, n1, p=pb1)*dbinom(b2, n2, p=pb1),
           p.ra0rb0=dbinom(a1, n1, p=pa0)*dbinom(a2, n2, p=pa0)*
             dbinom(b1, n1, p=pb0)*dbinom(b2, n2, p=pb0))


  bW.aPbP.smry <- apply(indicator_p[indicator_p$posterior> delta,
                                    c("p.ra1rb1","p.ra0rb0")],2,sum)


  ## expected total sample size under p0
  ess <- n1*2*table0[1,1]+
    (n1*2+n2*2)*(table0[2,2] +table0[2,3]+table0[3,2] +table0[3,3])+
    (n1*2+n2)*(table0[1,2] +table0[1,3]+table0[2,1] +table0[3,1])



  ## power and type I error of pick the true winner
  bpower_1 <- bW.aPbP.smry["p.ra1rb1"]
  bpower_0 <- bW.aPbP.smry["p.ra0rb0"]

  power <-  table1[1,3]+table1[2,3]+bpower_1 ## pa_1 pb_1
  alpha <-  table0[1,3]+table0[2,3]+bpower_0 ## pa_0 pb_0

  table0<- as.data.frame(table0)
  table0 <- round(table0,3)

  table0[3,3] <- paste0(table0[3,3], " (B claims as winner: ", round(bpower_0,3), ")")
  return(table0)

}




## server
library(shiny)
library(shinycssloaders)
library(rhandsontable)
library(dplyr)
library(Rcpp)
library(clinfun)
server <- function(input, output, session){



  ## explanations of Parameters
  output$Parameters1 <- renderText({
    paste0("Pa1 must be less than Pb1, and Pa0 must be less or equal than Pb0\nMaximum sample size for a signle arm is set to be 100")
  })

  ## Hypothesis
  output$Hypothesis1 <- renderText({
    paste("Hypothesis test:\nH0: PB0 = rateB0 v.s PA0 = rateA0 \nH1: PB1 = rateB1 v.s PA1 = rateA1 ")
  })

  output$explain1 <- renderText({
    paste("n is the sample size for a signle arm. n1 is the sample size for stage 1. n2 is the sample size for\nstage 2. If number of response > r1, then the arm pass stage 1. If number of responses > r, then the\narm pass stage 2. En_p0 is the expected total sample size of two arms under H0.")
  })



  ## run the simulation

  results1 <- eventReactive(input$runButton1, {
    req(input$Pa_1, input$Pb_1, input$alpha1, input$beta1, input$Pa_0, input$Pb_0, input$delta1)
    optimal_minimax(as.numeric(input$Pa_1), as.numeric(input$Pb_1), as.numeric(input$alpha1), as.numeric(input$beta1), as.numeric(input$Pa_0), as.numeric(input$Pb_0), as.numeric(input$delta1), scenarios = scenarios, posterior_matrix = posterior_matrix)
  })

  output$results1 <- DT::renderDT({
    results1()
  })


  # Results for Page 2
  ## instruction
  output$Instruction <- renderText({
    print("This allows user to evaluate operating characteristic of at most 5 scenarios at a same time. If you have less than 5 scenarios, please leave parameters in the extra scenarios at zero.\nHypothesis test: H1: Pa1 = rateA1, Pb1 = rateB1 v.s H0: Pa0 = rateA0, Pb1 = rateB0")

  })

  output$Instruction2 <- renderText({
    print(paste0("A.ES.H0: Probability of arm A early stopping under H0\n",
          "B.ES.H0: Probability of arm B early stopping under H0\n",
          "AB.ES.H0: Probability of arm A or B  early stopping under H0\n",
          "A.P1.H1: Probability of arm A passing stage I under H1\n",
          "B.ES.H1: Probability of arm B passing stage I under H1\n",
          "AB.ES.H1: Probability of arm A or B passing stage I under H1"))

  })
  
  ## instruction for page 3
  output$Instruction_page3 <- renderText({
    print("This is the decision rule according to the Optimal/Minimax design got from page1,\nplease run the simulation in page1 ahead.\nYa, Yb represent the total number of response in arm A and B respectively.")
    
  })


  # initial table
  initial_data <- data.frame(
    Scenario = c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5"),
    n = rep(0, 5),
    n1 = rep(0, 5),
    r = rep(0, 5),
    r1 = rep(0, 5),
    Pa0 = rep(0, 5),
    Pb0 = rep(0, 5),
    Pa1 = rep(0, 5),
    Pb1 = rep(0, 5),
    delta = rep(0, 5)
  )

  # Create a reactive data frame for the editable table
  edited_data <- reactiveVal(initial_data)

  # Render the editable table using rhandsontable
  output$tableInput <- renderRHandsontable({
    rhandsontable(edited_data(), params = list(allowInsertRow = TRUE, allowRemoveRow = TRUE))
  })


  # Results calculations
  results2 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)

    # Initialize an empty data frame to store results
    result_df <- data.frame()

    for (i in 1:5) { # Adjust the loop range to consider only scenarios 1 and 2
      if (edited_data_data$n[i]!=0) {
        result <- oc(
          as.numeric(edited_data_data$n[i]),
          as.numeric(edited_data_data$n1[i]),
          as.numeric(edited_data_data$r[i]),
          as.numeric(edited_data_data$r1[i]),
          as.numeric(edited_data_data$Pa1[i]),
          as.numeric(edited_data_data$Pb1[i]),
          as.numeric(edited_data_data$Pa0[i]),
          as.numeric(edited_data_data$Pb0[i]),
          as.numeric(edited_data_data$delta[i]),
          posterior_matrix
        )

        # Check if the result is not empty
        if (!is.null(result) && nrow(result) > 0) {
          # Add the result to the result_df
          result_df <- bind_rows(result_df, result)
        }
      }
    }

    return(result_df)
  })

  output$results2 <- DT::renderDT({
    results2()
  })


  ## table0 output result 3
  results31 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)

    table <- octable0(edited_data_data$n[1], edited_data_data$n1[1], edited_data_data$r[1], edited_data_data$r1[1], edited_data_data$Pa1[1], edited_data_data$Pb1[1], edited_data_data$Pa0[1], edited_data_data$Pb0[1], edited_data_data$delta[1], posterior_matrix)

    colnames(table) <- c("Scenario1: B.fail.stage1", "Scenario1: B.fail.stage2", "Scenario1: B.pass")
    rownames(table) <- c("Scenario1: A.fail.stage1", "Scenario1: A.fail.stage2", "Scenario1: A.pass")
    table

  })

  output$results31 <- DT::renderDT({
    results31()
  })

  results32 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)
    if(edited_data_data$n[2]!=0){
      table <- octable0(edited_data_data$n[2], edited_data_data$n1[2], edited_data_data$r[2], edited_data_data$r1[2], edited_data_data$Pa1[2], edited_data_data$Pb1[2], edited_data_data$Pa0[2], edited_data_data$Pb0[2], edited_data_data$delta[2], posterior_matrix)

      colnames(table) <- c("Scenario2: B.fail.stage1", "Scenario2: B.fail.stage2", "Scenario2: B.pass")
      rownames(table) <- c("Scenario2: A.fail.stage1", "Scenario2: A.fail.stage2", "Scenario2: A.pass")
      table
      }
  })


  output$results32 <- DT::renderDT({
    results32()
  })

  results33 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)
    if(edited_data_data$n[3]!=0){
      table <- octable0(edited_data_data$n[3], edited_data_data$n1[3], edited_data_data$r[3], edited_data_data$r1[3], edited_data_data$Pa1[3], edited_data_data$Pb1[3], edited_data_data$Pa0[3], edited_data_data$Pb0[3], edited_data_data$delta[3], posterior_matrix)
      colnames(table) <- c("Scenario3: B.fail.stage1", "Scenario3: B.fail.stage2", "Scenario3: B.pass")
      rownames(table) <- c("Scenario3: A.fail.stage1", "Scenario3: A.fail.stage2", "Scenario3: A.pass")
      table

       }
  })

  output$results33 <- DT::renderDT({
    results33()
  })

  results34 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)
    if(edited_data_data$n[4]!=0){
      table <- octable0(edited_data_data$n[4], edited_data_data$n1[4], edited_data_data$r[4], edited_data_data$r1[4], edited_data_data$Pa1[4], edited_data_data$Pb1[4], edited_data_data$Pa0[4], edited_data_data$Pb0[4], edited_data_data$delta[4], posterior_matrix)
      colnames(table) <- c("Scenario4: B.fail.stage1", "Scenario4: B.fail.stage2", "Scenario4: B.pass")
      rownames(table) <- c("Scenario4: A.fail.stage1", "Scenario4: A.fail.stage2", "Scenario4: A.pass")
      table

      }
  })

  output$results34 <- DT::renderDT({
    results34()
  })

  results35 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)
    if(edited_data_data$n[5]!=0){
      table <- octable0(edited_data_data$n[5], edited_data_data$n1[5], edited_data_data$r[5], edited_data_data$r1[5], edited_data_data$Pa1[5], edited_data_data$Pb1[5], edited_data_data$Pa0[5], edited_data_data$Pb0[5], edited_data_data$delta[5], posterior_matrix)
      colnames(table) <- c("Scenario5: B.fail.stage1", "Scenario5: B.fail.stage2", "Scenario5: B.pass")
      rownames(table) <- c("Scenario5: A.fail.stage1", "Scenario5: A.fail.stage2", "Scenario5: A.pass")
      table
      }
  })

  output$results35 <- DT::renderDT({
    results35()
  })

  ## table0 output result 4
  results41 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)

    table <- octable1(edited_data_data$n[1], edited_data_data$n1[1], edited_data_data$r[1], edited_data_data$r1[1], edited_data_data$Pa1[1], edited_data_data$Pb1[1], edited_data_data$Pa0[1], edited_data_data$Pb0[1], edited_data_data$delta[1], posterior_matrix)

    colnames(table) <- c("Scenario1: B.fail.stage1", "Scenario1: B.fail.stage2", "Scenario1: B.pass")
    rownames(table) <- c("Scenario1: A.fail.stage1", "Scenario1: A.fail.stage2", "Scenario1: A.pass")
    table

  })

  output$results41 <- DT::renderDT({
    results41()
  })

  results42 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)
    if(edited_data_data$n[2]!=0){
      table <- octable1(edited_data_data$n[2], edited_data_data$n1[2], edited_data_data$r[2], edited_data_data$r1[2], edited_data_data$Pa1[2], edited_data_data$Pb1[2], edited_data_data$Pa0[2], edited_data_data$Pb0[2], edited_data_data$delta[2], posterior_matrix)

      colnames(table) <- c("Scenario2: B.fail.stage1", "Scenario2: B.fail.stage2", "Scenario2: B.pass")
      rownames(table) <- c("Scenario2: A.fail.stage1", "Scenario2: A.fail.stage2", "Scenario2: A.pass")
      table
    }
  })

  output$results42 <- DT::renderDT({
    results42()
  })


  results43 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)
    if(edited_data_data$n[3]!=0){
      table <- octable1(edited_data_data$n[3], edited_data_data$n1[3], edited_data_data$r[3], edited_data_data$r1[3], edited_data_data$Pa1[3], edited_data_data$Pb1[3], edited_data_data$Pa0[3], edited_data_data$Pb0[3], edited_data_data$delta[3], posterior_matrix)
      colnames(table) <- c("Scenario3: B.fail.stage1", "Scenario3: B.fail.stage2", "Scenario3: B.pass")
      rownames(table) <- c("Scenario3: A.fail.stage1", "Scenario3: A.fail.stage2", "Scenario3: A.pass")
      table

    }
  })

  output$results43 <- DT::renderDT({
    results43()
  })

  results44 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)
    if(edited_data_data$n[4]!=0){
      table <- octable1(edited_data_data$n[4], edited_data_data$n1[4], edited_data_data$r[4], edited_data_data$r1[4], edited_data_data$Pa1[4], edited_data_data$Pb1[4], edited_data_data$Pa0[4], edited_data_data$Pb0[4], edited_data_data$delta[4], posterior_matrix)
      colnames(table) <- c("Scenario4: B.fail.stage1", "Scenario4: B.fail.stage2", "Scenario4: B.pass")
      rownames(table) <- c("Scenario4: A.fail.stage1", "Scenario4: A.fail.stage2", "Scenario4: A.pass")
      table

    }
  })

  output$results44 <- DT::renderDT({
    results44()
  })

  results45 <- eventReactive(input$runButton2, {
    edited_data_data <- hot_to_r(input$tableInput)
    if(edited_data_data$n[5]!=0){
      table <- octable1(edited_data_data$n[5], edited_data_data$n1[5], edited_data_data$r[5], edited_data_data$r1[5], edited_data_data$Pa1[5], edited_data_data$Pb1[5], edited_data_data$Pa0[5], edited_data_data$Pb0[5], edited_data_data$delta[5], posterior_matrix)
      colnames(table) <- c("Scenario5: B.fail.stage1", "Scenario5: B.fail.stage2", "Scenario5: B.pass")
      rownames(table) <- c("Scenario5: A.fail.stage1", "Scenario5: A.fail.stage2", "Scenario5: A.pass")
      table
    }
  })


  output$results45 <- DT::renderDT({
    results45()
  })
  
  ## results for page3
  results_dr <- eventReactive(input$runButton3, {
    req(results1(), input$Pa_1, input$Pb_1, input$Pa_0, input$Pb_0, input$delta1, input$design_select)
    decision_table(input$design_select, results1(),
                   as.numeric(input$Pa_1), as.numeric(input$Pb_1), as.numeric(input$Pa_0), as.numeric(input$Pb_0), as.numeric(input$delta1), posterior_matrix = posterior_matrix)
  })
  
  # Render results_dr as a datatable
  output$results_dr <- DT::renderDT({
    results_dr()
  })
  

}


