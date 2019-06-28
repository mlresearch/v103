### Schwaferts & Augustin 2019              ###

### Imprecise Hypothesis-based Bayesian     ###
### Decision Making with Simple Hypotheses. ###

### Electronic Appendix: Example            ###


# parameters
n = 10
m = 10
N = n + m
x = 9
y = 5
z = x + y
# hypotheses: 0.5 vs [0.75 , 0.9]

# lower bounds of LR
LR_x_U = dbinom(x,n, 0.5) / dbinom(x,n, 0.75)
LR_y_U = dbinom(y,m, 0.5) / dbinom(y,m, 0.9)
LR_z_U = dbinom(z,N, 0.5) / dbinom(z,N, 0.9)

# upper bounds of LR
LR_x_L = dbinom(x,n, 0.5) / dbinom(x,n, 0.9)
LR_y_L = dbinom(y,m, 0.5) / dbinom(y,m, 0.75)
LR_z_L = dbinom(z,N, 0.5) / dbinom(z,N, 0.75)

# imprecise likelihood ratios
c(LR_x_L, LR_x_U)
c(LR_y_L, LR_y_U)
c(LR_z_L, LR_z_U)

# updating inconsistency
c(LR_y_L * LR_x_L, LR_y_U* LR_x_U)
c(LR_z_L, LR_z_U)
c(4.214 * 0.025, 165.4 * 0.052) # for paper: using rounded values

# imprecise loss function
k_L = 8
k_U = 20

# imprecise prior odds
pi_L = 1
pi_U = 4

# imprecise ratio of expected posterior losses for x
c(pi_L * k_L * LR_x_L, pi_U * k_U * LR_x_U)


# imprecise ratio of expected posterior losses for x and y
# with inconsistencies
c(pi_L * k_L * LR_y_L * LR_x_L, pi_U * k_U * LR_y_U * LR_x_U) 
c(pi_L * k_L * 4.214 * 0.025, pi_U * k_U * 165.4 * 0.052) # for paper


# imprecise ratio of expected posterior losses for z
c(pi_L * k_L * LR_z_L, pi_U * k_U * LR_z_U)


# arbitrary precise parameters
pi = 1
k = k_L
# hypotheses: 0.5 vs. 0.8

# precise likelihood ratio for x
LR_x =  dbinom(x,n, 0.5) / dbinom(x,n, 0.8)
LR_x

# precise ratio of expected posterior losses for x
pi * k * LR_x

