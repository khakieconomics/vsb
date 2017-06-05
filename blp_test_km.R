library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


blpdir <- "KM_blp_data/"

cdid <- read_delim(paste(blpdir,"/cdid.txt",sep=""), "\t",col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
iv <- read_delim(paste(blpdir,"/IV1.txt",sep=""), "\t",col_names = FALSE, escape_double = FALSE, trim_ws = TRUE, col_types = "nnnnnnnnnnnnnnn")
x1 <- read_delim(paste(blpdir,"/x1.txt",sep=""), "\t",col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
x2 <- read_delim(paste(blpdir,"/x2.txt",sep=""), "\t",col_names = FALSE, escape_double = FALSE, trim_ws = TRUE)
shares <- read_delim(paste(blpdir,"/share.txt",sep=""), col_names = FALSE,"\t", escape_double = FALSE, trim_ws = TRUE, col_types = "c") %>% 
  mutate(X1 = as.numeric(X1))
outshare <- read_delim(paste(blpdir,"/outshr.txt",sep=""), col_names = FALSE,"\t", escape_double = FALSE, trim_ws = TRUE)
mkt_prod_counts <- as.numeric(table(cdid))
mkt_ids <- as.factor(cdid$X1)

product_numbers <- cdid %>% 
  group_by(X1) %>% 
  mutate(N = n())
# lnprice <- x1[,1]
# X_plain <- cbind(scale(x1[,2:ncol(x1)]))
X_plain <- cbind(x1[,2:ncol(x1)], product_n = product_numbers$N/100)
X_rc <- x2[,2:ncol(x2)]
Z_iv <- iv[,-1:-4]

# Indices
T <- length(levels(mkt_ids))
P_plain <- ncol(X_plain)
P_rc <- ncol(X_rc)
P_iv <- ncol(Z_iv)

J_m <- nrow(x1)
NS <- 100

market_size <- round(rpois(T, 30000))

df <- data.frame(cdid = cdid[[1]], price = x1[[1]], shares=shares[[1]], outshare = outshare[[1]])


df$market_size <- rep(market_size, mkt_prod_counts)

xx <- function(x) {
  outside_share = first(x$outshare)
  bind_rows(x, data_frame(cdid =x$cdid[1], market_size = market_size[1], price = 0, shares = outside_share, outshare = outside_share)) 
}

aug_price_and_sales <- df %>% 
  group_by(cdid) %>% 
  do(xx(.))

aug_price_and_sales$sales <- round(aug_price_and_sales$shares * aug_price_and_sales$market_size)




# random shocks
beta_shocks <- matrix(rnorm(NS*P_rc), NS, P_rc)
alpha_shocks <- rnorm(NS)

# random shocks. (alpha_shocks, beta_shocks) = z_t
z <- cbind(alpha_shocks, beta_shocks)

price_mod_data <- as.data.frame(cbind(price = as.vector(df$price), X_plain, Z_iv[,-1]))
names(price_mod_data) <- c("Price", paste0("X", 1:15))

price_mod <- lm(log(Price) ~., data = price_mod_data)

residual_sd <- summary(price_mod)$sigma

data_list <- list(NS = NS, 
                  J_m = J_m, 
                  T = T, 
                  P_rc = P_rc, 
                  P_plain = P_plain, 
                  P_iv = P_iv-1,
                  price = as.vector(df$price), 
                  sales = as.vector(aug_price_and_sales$sales),
                  X_rc = X_rc,
                  X_plain = X_plain,
                  Z_iv = Z_iv[,-1],
                  mkt_size = as.vector(market_size),
                  mkt_numProds = as.vector(mkt_prod_counts),
                  mkt_id = as.vector(cdid$X1),
                  z = z,
                  nu = 3,
                  run_estimation=1,
                  residual_sd = residual_sd,
                  P0 = 4)




library(rstan)
options(mc.cores = parallel::detectCores())
compiled_model <- stan_model("vsb_flex_iv_horseshow.stan")

# Uncomment the below based on what estimates you want

#system.time(test_optim <- optimizing(compiled_model, data = data_list, iter=10000, verbose = T)) 

#hmc_test <- sampling(compiled_model, data = data_list, iter = 400, chains = 4)
#save(hmc_test, file = "blp_hmc_20170601.RData")
#shinystan::launch_shinystan(hmc_test)
# test_optim$value
# 
# vsb_xi <- test_optim$par[grepl("xi", names(test_optim$par))]
# hist(vsb_xi)
# test_optim$par[grepl("alpha", names(test_optim$par))]
# test_optim$par[grepl("eta", names(test_optim$par))]
# test_optim$par[grepl("lambda", names(test_optim$par))]
# test_optim$par[grepl("price", names(test_optim$par))]
# 
# hist(test_optim$par[grepl("mu", names(test_optim$par))])
# plot(vsb_xi,log(data_list$price))
# 
# plot(get_posterior_mean(hmc_test, pars = "xi")[,5], log(data_list$price))
# hist(get_posterior_mean(hmc_test, pars = "xi")[,5])
# alpha_estimates <- data_frame(alpha = (extract(hmc_test, pars = "alpha")[[1]]))
# km_estimates <- data_frame(KM_estimates = c(-.573, -.644, -.666, -.725, -.664, -.587, -.593, -.898, -.896, -.646),
#                            GMM_obj = c(195.6, 156, 147.5, 158.3, 148.2, 292.3, 230.2, 99.9, 154.7, 141.6))
# library(ggplot2) 
# 
# alpha_estimates %>% 
#   ggplot(aes(x = alpha)) +
#   geom_density() +
#   ggthemes::theme_hc() +
#   scale_colour_grey() +
#   geom_point(data = km_estimates, aes(x = KM_estimates, y = 0.5, size = 1/GMM_obj), alpha = 0.6) +
#   labs(x = "Value of alpha", y = "Posterior density")
# 
# Omegas <- test_optim$par[grepl("Omega", names(test_optim$par))]
# Scales <- test_optim$par[grepl("^scale", names(test_optim$par))]
# pars <- c(test_optim$par[1:(P_plain+1)], Omegas, Scales)
# alpha <- test_optim$par[grepl("alpha", names(test_optim$par))]
# 
# data_frame(Estimates = pars, 
#            Parameter = c("alpha", rep("beta", P_plain), rep("Omega", (P_rc+1)^2), rep("scale", P_rc+1) )) %>% 
#   ggplot(aes(x = seq_along(Estimates), y = Estimates)) + 
#   geom_point(aes(colour = Parameter)) + 
#   ggthemes::theme_economist() +
#   labs(title = "Estimates and true values\nof structural parameters")
