one_iteration <- function(start, p_theta){
   step <- sample(c(-1, 1), 1, replace = TRUE, prob = c(0.5, 0.5))
   proposal <- start + step
   if (proposal < 1){
      proposal <- length(p_theta)
   }
   if (proposal > length(p_theta)){
      proposal <- 1
   }
   accept_rate <- min(1, p_theta[proposal] / p_theta[start])
   new <- start
   if(runif(1) < accept_rate){
     new <- proposal
   }
   return(new)

}

MH_discrete <- function(p_theta, Niter = 1000, burnin = 500){

  samples_temp <- rep(NA, Niter)
  samples_temp[1] <- sample(1:length(p_theta), 1, replace = TRUE, prob = rep(1, length(p_theta)))
  for(i in 2:Niter){
      samples_temp[i] <- one_iteration(samples_temp[i - 1], p_theta)
  }
  return(samples_temp[burnin:Niter])
}


p_theta <- 1:7
Niter <- 500
burnin <- floor(0.2 * Niter)
samples_MH_discrete <- MH_discrete(p_theta, Niter, burnin)


hist(samples_MH_discrete, 20, freq = FALSE)
plot(p_theta/sum(p_theta), 
     as.vector(table(samples_MH_discrete)/length(samples_MH_discrete)), 
     xlab = "target", ylab = "Metropolis Samples")
lines(seq(0, 0.4, length.out = 20), seq(0, 0.4, length.out = 20))

print(cbind(p_theta/sum(p_theta), 
            as.vector(table(samples_MH_discrete)/length(samples_MH_discrete))))
install.packages("tidyverse")
library(tidyverse)
a = c(1,2,3,1,2,4)
min_rank()

test= tibble(
   a = c(1,2,3,4,5,5,5,2,2),
   b = c(2,3,3,3,4,4,4,2,2)
)

c = test %>% top_n(5,a)
c

