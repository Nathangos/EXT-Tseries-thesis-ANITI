### TEST EGPD implementation
rm(list=ls())
set.seed(133)
require(extRemes)
require(mev)
require(goftest)
M<-500

### Target: true shape parameter
shapeF<-1
Pareto_sim<-extRemes::revd(n = M,scale = 1,threshold = 1,shape = shapeF,
                           type = "GP")
print(goftest::ad.test(x = Pareto_sim,null = extRemes::"pevd",
                       threshold=1,shape=shapeF,scale=1,
                       type="GP"))
### EGPD=GPD if kappa=1
Init<-mev::fit.extgp(data = Pareto_sim,init = c(1,1,shapeF),
                     method = "pwm")
Theta_PWM<-Init$fit$pwm
mle_opt<-mev::fit.extgp(data = Pareto_sim,init =Theta_PWM,
               method = "mle")
print(mle_opt$fit$mle)
print("very different from the truth")

### with the log
LOG_value<-log(Pareto_sim)
Init<-mev::fit.extgp(data =LOG_value,init = c(1,1,0),
                     method = "pwm")
Theta_PWM<-Init$fit$pwm
opt_MLE<-mev::fit.extgp(data =LOG_value,init = Theta_PWM,
               method = "mle")
Theta_mle<-opt_MLE$fit$mle
print(Theta_mle)
print("Better results")

### Convert to Unif
Unif_<-mev::pextgp(q =LOG_value,kappa = Theta_mle[["kappa"]],
            sigma = Theta_mle[["sigma"]],
            xi = Theta_mle[["xi"]])
print(summary(Unif_))
print("AD test for the uniform")
goftest::ad.test(Unif_,null = "punif")
