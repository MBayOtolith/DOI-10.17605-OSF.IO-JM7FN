library(MASS)

par(mfrow = c(3,3))
# evaluate power vs log/exp
# MARINE
# log/exp model
sim_joy_marine$pred_vars <- as.numeric(sim_joy_marine$pred_vars)
plot(sim_joy_marine$pred_vars,sim_joy_marine$value, col = "blue4")
exp_mod <- lm(log(value) ~ pred_vars, data = sim_joy_marine)
summary(exp_mod)
exp_pred <- predict(exp_mod, newdata = data.frame(pred_vars = seq(1,30,1))) %>% exp()

lines(exp_pred, col = "blue4", lwd = 2)
#points(sim_joy_marine$pred_vars,predict(exp_mod) %>% exp(), col = "red")

# power function
power_mod <- nls(value~b*pred_vars^z,start = list(b = .1, z = 1),data=sim_joy_marine, control = list(maxiter = 1000))
summary(power_mod)
pow_pred <- predict(power_mod, newdata = data.frame(pred_vars = seq(1,30,1)))
lines(seq(1,30,1),pow_pred, col = "blue4", lty = "dashed", lwd = 2)
#points(sim_joy_marine$pred_vars,predict(power_mod), col = 'blue')

plot(predict(power_mod),sim_joy_marine$value)
abline(a = 0,b=1, col = "red")

plot(predict(exp_mod) %>% exp(),sim_joy_marine$value)
abline(a = 0,b=1, col = "red")


# evaluate fits
mean(residuals(exp_mod)^2) %>% sqrt() # this isnt back transformed
mean(residuals(power_mod)^2) %>% sqrt()

# r2
cor(sim_joy_marine$value,predict(power_mod)) ^ 2
cor(sim_joy_marine$value,predict(exp_mod) %>% exp()) ^ 2

#mse
mean((sim_joy_marine$value - predict(power_mod))^2)
mean( (sim_joy_marine$value - (predict(exp_mod) %>% exp()) )^2)

#mae
mean((sim_joy_marine$value - predict(power_mod)))
mean( (sim_joy_marine$value - (predict(exp_mod) %>% exp()) ))


# evaluate power vs log/exp
# LAND
# log/exp model
sim_joy_land$pred_vars <- as.numeric(sim_joy_land$pred_vars)
plot(sim_joy_land$pred_vars,sim_joy_land$value, col = "green4", cex = .5)
exp_mod <- lm(log(value) ~ pred_vars, data = sim_joy_land)
summary(exp_mod)
lines(predict(exp_mod, newdata = data.frame(pred_vars = seq(1,30,1))) %>% exp(), col = "green4", lwd = 2)
#points(sim_joy_land$pred_vars,predict(exp_mod) %>% exp(), col = "red")

# power function
power_mod <- nls(value~b*pred_vars^z,start = list(b = .1, z = 1),data=sim_joy_land, control = list(maxiter = 1000))
summary(power_mod)
pow_pred <- predict(power_mod, newdata = data.frame(pred_vars = seq(1,30,1)))
lines(seq(1,30,1),pow_pred, col = "green4", lty = "dashed", lwd = 2)
#points(sim_joy_land$pred_vars,predict(power_mod), col = 'blue')

plot(predict(power_mod),sim_joy_land$value)
abline(a = 0,b=1, col = "red")

plot(predict(exp_mod) %>% exp(),sim_joy_land$value)
abline(a = 0,b=1, col = "red")

# eval fits
# r2
cor(sim_joy_land$value,predict(power_mod)) ^ 2
cor(sim_joy_land$value,predict(exp_mod) %>% exp()) ^ 2

#mse
mean((sim_joy_land$value - predict(power_mod))^2)
mean( (sim_joy_land$value - (predict(exp_mod) %>% exp()) )^2)

#mae
mean((sim_joy_land$value - predict(power_mod)))
mean( (sim_joy_land$value - (predict(exp_mod) %>% exp()) ))



# 
# plotting both log models on top
sim_joy_land$pred_vars <- as.numeric(sim_joy_land$pred_vars)
sim_joy_marine$pred_vars <- as.numeric(sim_joy_marine$pred_vars)

plot(sim_joy_land$pred_vars ,sim_joy_land$value%>% log(), col = "green4", cex = 1)

points(sim_joy_marine$pred_vars,sim_joy_marine$value%>% log(), col = "blue4", cex = .5)

exp_mod_land <- lm(log(value) ~ pred_vars, data = sim_joy_land)
exp_mod_mar <- lm(log(value) ~ pred_vars, data = sim_joy_marine)

lines(predict(exp_mod_land, newdata = data.frame(pred_vars = seq(1,30,1))) , col = "green4", lwd = 2)
lines(predict(exp_mod_mar, newdata = data.frame(pred_vars = seq(1,30,1))) , col = "blue4", lwd = 2)

plot(predict(exp_mod_mar, newdata = data.frame(pred_vars = seq(1,30,1))) , col = "blue4", lwd = 2, type = "l")
lines(predict(exp_mod_land, newdata = data.frame(pred_vars = seq(1,30,1))) , col = "green4", lwd = 2, type = "l")


ggplot()+
  geom_histogram(aes(x=sim_joy_land$value%>% log()), fill = "green4", bins = 100)+
  geom_histogram(aes(x=sim_joy_marine$value%>% log()), fill = "blue4", bins = 100, alpha = .5, color = "black")+
  coord_flip()+
  themeo

# same as above just on a log
plot(1:20,exp(1:20))
points(1:20,10^(1:20))


# evaluate power vs log/exp
sim_joy_land$pred_vars <- as.numeric(sim_joy_land$pred_vars)

plot(sim_joy_land$pred_vars ,sim_joy_land$value%>% log())
exp_mod <- lm(log(value) ~ pred_vars, data = sim_joy_land)
summary(exp_mod)
lines(predict(exp_mod, newdata = data.frame(pred_vars = seq(1,30,1))) , col = "red")

points(sim_joy_land$pred_vars,predict(exp_mod) , col = "red")

mean(exp_mod$residuals^2)

library(MASS)
power_mod <- nls(value~b*pred_vars^z,start = list(b = .1, z = 1),data=sim_joy_land, control = list(maxiter = 1000))
summary(power_mod)
pow_pred <- predict(power_mod, newdata = data.frame(pred_vars = seq(1,30,1)))
lines(seq(1,30,1),pow_pred %>% log(), col = "blue")
points(sim_joy_land$pred_vars,predict(power_mod)%>% log(), col = 'blue')










# evaluate power vs log/exp
sim_joy_marine$pred_vars <- as.numeric(sim_joy_marine$pred_vars)

plot(sim_joy_marine$pred_vars ,sim_joy_marine$value%>% log())
exp_mod <- lm(log(value) ~ pred_vars, data = sim_joy_marine)
summary(exp_mod)
lines(predict(exp_mod, newdata = data.frame(pred_vars = seq(1,30,1))) , col = "red")

points(sim_joy_marine$pred_vars,predict(exp_mod) , col = "red")

mean(exp_mod$residuals^2)

library(MASS)
power_mod <- nls(value~b*pred_vars^z,start = list(b = .1, z = 1),data=sim_joy_marine, control = list(maxiter = 1000))
summary(power_mod)
pow_pred <- predict(power_mod, newdata = data.frame(pred_vars = seq(1,30,1)))
lines(seq(1,30,1),pow_pred %>% log(), col = "blue")
points(sim_joy_marine$pred_vars,predict(power_mod)%>% log(), col = 'blue')
















