library(ggplot2)
library(reshape2)


#Per time step, generates the new generation
#i is population of A, N is total population
#Ra and Rb are fitness of A and B
generator <- function (i, N, Ra, Rb) {
  #Prob of A giving birth
  p.bir <- (Ra*i)/((Ra*i) + (Rb*(N-i)))
  
  #Prob of A dying
  p.death <- (i/N)
  
  
  #Selecting who gives birth
  if (runif(1) < p.bir) {
    i <- i + 1 #If we select a, then a grows by 1
  } else {
    i <- i #If we select b, then a is steady
  }
  
  
  #Selecting who dies
  if (runif(1) < p.death) {
    i <- i - 1 #If we select a, then a reduces by 1
  } else {
    i <- i #If we select b, then a is steady
  }
  
  #Spit out number of a
  return(i)
}

#Takes the generator function and performs it across
#Many time steps for a specified number of reps
#Spits out a Mean with upper and lower confidence bounds
moran.process <- function (i, N, Ra, Rb, rep, step) {

#Sets up a matrix where each row is a time step
#Each column is a replicate
mat <- matrix(nrow = step, ncol = rep)

#Recording the origial i for each rep reset
orig.i <- i

#Loop, for each p (a rep), calculate the needed generations
#Then reset back the original i and repeat for a new column
for (p in 1:rep) {
  i <- orig.i
  for (q in 1:step) {
    mat[q, p] <- i
    i <- generator(i, N = N, Ra = Ra, Rb = Rb)
  }
}

#Generating a summary matrix with mean and confidence intervals
sum.mat <- matrix(nrow = step, ncol = 3)

#First column is mean, and second column are standard error
sum.mat[,1] <- apply(mat, 1, mean)
sum.mat[,2] <- apply(mat, 1, function (x) (2*(sd(x))/sqrt(rep)))
sum.mat[,3] <- apply(mat, 1, function (x) (2*(sd(x))/sqrt(rep)))

#Lower Bound (95%)
sum.mat[,2] <- sum.mat[,1] - sum.mat[,2]

#Upper Bound (95%)
sum.mat[,3] <- sum.mat[,1] + sum.mat[,3]

#Naming columns
colnames(sum.mat) <- c("mean", "lower_bound", "upper_bound")

return(as.data.frame(sum.mat))
}

#testing this out
test <- moran.process(1, 10, 1, 2, 10, 10000)

#test plot
ggplot(test, aes(x = as.numeric(row.names(test)), y = mean)) +
    geom_line()+
    geom_ribbon(aes(ymin=lower_bound,ymax=upper_bound), alpha=0.3)


#Set of populations
pop <- c(10, 100, 1000, 10000)

#Set of initial proportions of A
init.prop <- c(.3, .5, .7)

#Number of time steps
step <- 5000

vars <- as.data.frame(expand.grid(
                              N = pop,
                              i = init.prop, 
                              Ra = c(0.5, 0.75, 0.9, 0.99, 1, 1.01, 1.1, 1.5, 2),
                              Rb = 1, 
                              rep = 20, 
                              step = step))

#True initial i based off proportion of N
vars$i <- vars$i*vars$N

#applying all the variables to function
master <- mapply(moran.process, i = vars$i, 
                                N = vars$N, 
                                Ra = vars$Ra, 
                                Rb = vars$Rb, 
                                rep = vars$rep, 
                                step = vars$step, SIMPLIFY = FALSE)

#creating a dataframe 
master.sum <- do.call(rbind, master)

#Appending original variable  to the results
for (q in 1:nrow(vars)) {
  master.sum$i[(1+step*(q-1)):(step*q)] <- vars$i[q]
  master.sum$N[(1+step*(q-1)):(step*q)] <- vars$N[q]
  master.sum$R[(1+step*(q-1)):(step*q)] <- vars$Ra[q]
  master.sum$time[(1+step*(q-1)):(step*q)] <- 1:step
}

#Normalizing initial population and mean
master.sum$mean <- master.sum$mean/master.sum$N
master.sum$upper_bound <- master.sum$upper_bound/master.sum$N
master.sum$lower_bound <- master.sum$lower_bound/master.sum$N
master.sum$i <- master.sum$i/master.sum$N

#Converting initial pop, total pop size, and fitness ratio to factors
master.sum$i <- as.factor(master.sum$i)
master.sum$N <- as.factor(master.sum$N)
master.sum$R <- as.factor(master.sum$R)

#Naming Things Well
levels(master.sum$i) <- c("Init. Proportion of A = .3",
                          "Init. Proportion of A = .5",
                          "Init. Proportion of A = .7")

levels(master.sum$N) <- c("Total Population = 10",
                          "Total Population = 100",
                          "Total Population = 1000",
                          "Total Population = 10000")





#Plotting everything
fig.1 <- ggplot(master.sum, aes(x = time, y = mean, color=R,
                       ymin=lower_bound, ymax=upper_bound, fill = R)) +
  geom_line() +
  geom_ribbon(alpha=0.3, color = NA) +
  facet_grid(i ~ N) +
  labs(color = "Fitness Ratio \nof A/B",
       fill = "Fitness Ratio \nof A/B",
       x = "Time Steps",
       y = "Proportion of A in Population",
       subtitle = "")  + 
  guides(fill = guide_legend(reverse = TRUE),
         color = guide_legend(reverse = TRUE))

ggsave("Moran_Eash.png", plot = fig.1,
       scale = 1, width = 12, height = 8, units = "in",
       dpi = 400)
 



