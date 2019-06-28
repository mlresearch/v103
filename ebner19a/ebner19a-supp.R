# Ebner, Schwaferts, Augustin 2019:
# Robust Bayes Factor for independent two-sample comparisons
# Supplementary Material

library(ggplot2)
library(metR)

##### Functions #####

getBF_NP <- function(x1,x2,mu_d=0, var_d=1){
  #get Bayes Factor with Normal test-relevant effect size Prior
  
  # Takes the two samples and the hyperparameters for a normal test-relevant
  # effect size prior in a two-sample case with identical variance and the
  # hypotheses H0: delta = 0 vs. H1: delta ~ N(mu_d, var_d) and returns the
  # Bayes Factor BF_10.
  # Formula according to Gonan et al. 2005.
  
  t <- t.test(x1,x2)$statistic
  n1 <- length(x1)
  n2 <- length(x2)
  eff <- 1 / ((1/n1) + (1/n2)) #effective sample size
  
  postv <- 1 + (eff * var_d)
  nc <- sqrt(eff) * (mu_d / sqrt(1 + (eff * var_d)))
  
  BF = dt(t,n1+n2-2)/(dt(t/sqrt(postv), n1+n2-2,nc)/sqrt(postv))
  names(BF) <- NULL
  
  return(1/BF)
}

# The function to be optimized needs the parameters as first argument
getBF_to_optimize = function(par, x1, x2){
  #par: vector of parameters
  mu_d = par[1]
  var_d = par[2]
  BF = getBF_NP(x1=x1,x2=x2,mu_d = mu_d, var_d = var_d)
  return(BF)
}


##### Example #####

# True Parameters (so that the data is already standardized)
mu = 0 
sigma_sq = sigma = 1
alpha = delta = 0.2

n = m = 10
n2 = m2 = 30

# Generate Samples
set.seed(2102)
x = rnorm(n, mean = mu, sd = sigma)
y = rnorm(m, mean = mu + alpha, sd = sigma)
x2 = c(x, rnorm(n2 - n, mean = mu, sd = sigma))
y2 = c(y, rnorm(m2 - m, mean = mu, sd = sigma))

# Hypotheses
mu_d_L = 0
mu_d_U = 0.5
var_d_L = 0.5
var_d_U = 3

# Analysis
t.test(x,y)$statistic # t-value
1 / ((1/n) + (1/m)) # effective sample size

BF_L = optim(c(0,1), getBF_to_optimize, x1=x, x2=y, method= "L-BFGS-B",
             lower = c(mu_d_L,var_d_L), upper=c(mu_d_U,var_d_U))$value

BF_U = optim(c(0,1), getBF_to_optimize, x1=x, x2=y, method= "L-BFGS-B",
             lower = c(mu_d_L,var_d_L), upper=c(mu_d_U,var_d_U),
             control=list(fnscale=-1))$value

rBF = c(BF_L, BF_U)
rBF

# Analysis 2
t.test(x2,y2)$statistic # t-value
1 / ((1/n2) + (1/m2)) # effective sample size

BF_L2 = optim(c(0,1), getBF_to_optimize, x1=x2, x2=y2, method= "L-BFGS-B",
             lower = c(mu_d_L,var_d_L), upper=c(mu_d_U,var_d_U))$value

BF_U2 = optim(c(0,1), getBF_to_optimize, x1=x2, x2=y2, method= "L-BFGS-B",
             lower = c(mu_d_L,var_d_L), upper=c(mu_d_U,var_d_U),
             control=list(fnscale=-1))$value

rBF2 = c(BF_L2, BF_U2)
rBF2

1/c(0.18, 0.42)

##### Plots #####

N = 100   # number of grid points
mu_ds = seq(mu_d_L,mu_d_U, length.out=N+1)
var_ds = seq(var_d_L,var_d_U, length.out=N+1)

# initialize data frame
plot_df <- data.frame(BF=double(),
                      mu_d=double(),
                      var_d=double(),
                      data=factor())

# calculate BF values
for (mu_d in mu_ds) {
  for (var_d in var_ds) {
    BF1 = getBF_NP(x1=x,x2=y,mu_d = mu_d, var_d = var_d)
    BF2 = getBF_NP(x1=x2,x2=y2,mu_d = mu_d, var_d = var_d)
    new_rows = data.frame(BF=c(BF1,BF2),
                          mu_d=c(mu_d,mu_d),
                          var_d=c(var_d,var_d),
                          data=c("Analysis 1","Analysis 2"))
    plot_df = rbind(plot_df, new_rows)
  }
}

# plot 1
ggplot(data = plot_df[plot_df$data=="Analysis 1",], aes(x = mu_d, y = var_d, z=BF)) +
  geom_tile(aes(fill=BF)) +
  geom_contour(color="white") +
  geom_text_contour(stroke = 0.3, size=4.5) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14)) +
  xlab("Mean of Effect Size Prior") +
  ylab("Variance of Effect Size Prior") +
  ggtitle("Example: Analysis 1")

ggsave("plot1.pdf", width = 6, height = 6)


# plot 2
ggplot(data = plot_df[plot_df$data=="Analysis 2",], aes(x = mu_d, y = var_d, z=BF)) +
  geom_tile(aes(fill=BF)) +
  geom_contour(color="white") +
  geom_text_contour(stroke = 0.3, size=4.5) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14)) +
  xlab("Mean of Effect Size Prior") +
  ylab("Variance of Effect Size Prior") +
  ggtitle("Example: Analysis 2")

ggsave("plot2.pdf", width = 6, height = 6)
