library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
####Functions####
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

#Repeats generator function until i = 0 or N
#Returns number of time steps to extinction for x reps
extinguisher <- function (i, N, Ra, Rb, rep) {
  #Matrix to capture time steps to extinction
  mat <- matrix(nrow = rep, ncol = 1)
  
  #Retaining original i
  orig.i <- i
  

  
  #Run Moran Process Until i = 0 or i = N
  for (x in 1:rep) {
    #Start counter at 1, resets with each loop
    #logic is that if at i = N or 0, takes 0 steps to fix
    count <- 0
    
    #Reset i to original value
    i <- orig.i
    
    while (i != 0 & i != N)
    {
      i <- generator(i, N = N, Ra = Ra, Rb = Rb)
      count <- count + 1
    }
    
    #Record value in capture matrix
    mat[x, 1] <- count
    
  }
  
return(mat)  
  
}

#Repeats generator function until i = N
#If i = 0, restarts the function until i = N
#Returns average number of time steps to fixation of i
extinguisher.mut <- function (i, N, Ra, Rb, rep) {
  #Matrix to capture time steps to extinction
  mat <- matrix(nrow = rep, ncol = 1)
  
  #Retaining original i
  orig.i <- i
  
  #Run Moran Process Until i = 0 or i = N
  for (x in 1:rep) {
    #Start counter at 0, resets with each loop
    #logic is that if at i = N or 0, takes 0 steps to fix
    count <- 0
    
    #Reset i to original value
    i <- orig.i
    
    while (i != N) { #Keep loop going until i = N
      i <- generator(i, N = N, Ra = Ra, Rb = Rb) #Generate new i
      if (i == 0) { #if i goes extinct, gotta reset everything
        count <- 0 #count resets to 0
        i <- orig.i #i resets to original value
      } else {
        count <- count + 1 #If i didn't go to 0, add 1 to count
      }
    }
    #Record value in capture matrix
    mat[x, 1] <- count
    
  }
  
  #Returns average and standard deviation
  return(c(mean(mat), sd(mat)))
  
}


#Repeats generator function for a specified number of steps
#Asks if extinction was reached
#Reports proportion of reps that reached extinction 
fixer <- function (i, N, Ra, Rb, rep, step) {
#Sets counter to 0
count <- 0

#Record original i
orig.i <- i

  for (p in 1:rep) { #For how many ever reps
    i <- orig.i #Ensure i is set to original value
    for (q in 1:step) {
      i <- generator(i, N = N, Ra = Ra, Rb = Rb) #Run Moran Process for specified steps
    }
    if (i == N | i == 0) { #if i reaches extinaction, add 1 to counter
      count <- count + 1} else { #if it didn't, add zero to counter
        count <- count + 0
      }
  }

prop <- (count/rep) #Proportion of reps that reached extinction
return(prop)
}


#Repeats generator function for a specified number of steps
#Asks specifically if i (the mutant) has fixed
#Reports proportion of reps that reached fixation of i
fixer.mut <- function (i, N, Ra, Rb, rep) {
  #Sets counter to 0
  count <- 0
  
  #Record original i
  orig.i <- i
  
  for (p in 1:rep) { #For how many ever reps
    i <- orig.i #Ensure i is set to original value
    while (i != 0 & i != N)
    {
      i <- generator(i, N = N, Ra = Ra, Rb = Rb)
    }
    if (i == N) { #if i reaches fixation, add 1 to counter
      count <- count + 1} else { #if it didn't, add zero to counter
        count <- count + 0
      }
  }
  
  prop <- (count/rep) #Proportion of reps that reached extinction
  return(prop)
}


####Question 1: Basic Process over different N, N_0, and r####
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

#Changing order of levels for plot
master.sum$i <- factor(master.sum$i, levels = c("Init. Proportion of A = .7",
                                                "Init. Proportion of A = .5",
                                                "Init. Proportion of A = .3")) 





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

#Saving as PNG
ggsave("Moran_Eash.png", plot = fig.1,
       scale = 1, width = 12, height = 8, units = "in",
       dpi = 400)

####Question 2: Extinction Times over different N, N_0 and r####
#Set of populations
pop <- c(10, 50, 100, 500)

#Set of initial proportions of A
init.prop <- c(.3, .5, .7)

rep <- 100

vars <- as.data.frame(expand.grid(
  N = pop,
  i = init.prop, 
  Ra = c(0.5, 0.75, 0.9, 0.99, 1, 1.01, 1.1, 1.5, 2),
  Rb = 1, 
  rep = 100))

#True initial i based off proportion of N
vars$i <- vars$i*vars$N

#applying all the variables to function
master <- mapply(extinguisher, i = vars$i, 
                 N = vars$N, 
                 Ra = vars$Ra, 
                 Rb = vars$Rb, 
                 rep = vars$rep, SIMPLIFY = FALSE)

master.sum <- do.call(rbind, master)
master.sum <- as.data.frame(master.sum)

#Appending original variable  to the results
for (q in 1:nrow(vars)) {
  master.sum$i[(1+rep*(q-1)):(rep*q)] <- vars$i[q]
  master.sum$N[(1+rep*(q-1)):(rep*q)] <- vars$N[q]
  master.sum$R[(1+rep*(q-1)):(rep*q)] <- vars$Ra[q]
  master.sum$rep[(1+rep*(q-1)):(rep*q)] <- 1:rep
}


#Converting i back to porpotion 
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
                          "Total Population = 50",
                          "Total Population = 100",
                          "Total Population = 500")

#Changing order of levels for plot
master.sum$i <- factor(master.sum$i, levels = c("Init. Proportion of A = .7",
                                                "Init. Proportion of A = .5",
                                                "Init. Proportion of A = .3")) 


fig.2 <- ggplot(master.sum, aes(x = R, y = V1, color = i)) +
  scale_y_continuous(trans='log10') +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(~N) +
  labs(color = "Initial Proportion of A",
       x = "Relative Fitness of A over B",
       y = "Average Extinction Time (LOG SCALE)",
       subtitle = "")  

#Saving as PNG
ggsave("extinct_time_Eash.png", plot = fig.2,
       scale = 1, width = 16, height = 8, units = "in",
       dpi = 400)

#####Question 3: Recreating Graphs####
#Graph 1 Fix probability vs Selective Advantage of Mutant
#N = 50
#Sims = 1500
#A = 1
#Selective Advtange ranges from 1 to 5 by .2
#Standard Error of Proportion
vars <- as.data.frame(expand.grid(
  N = 50,
  i = 1, 
  Ra = seq(1, 5 ,.2),
  Rb = 1, 
  rep = 1000))

fig3a.dat <- mapply(fixer.mut, i = vars$i, N = vars$N, 
                    Ra = vars$Ra, Rb = vars$Rb, 
                    rep = vars$rep, SIMPLIFY = TRUE)


fig3a.dat.form <- as.data.frame(cbind(vars, fig3a.dat))
fig3a.dat.form$err <- sqrt((fig3a.dat.form$fig3a.dat*(1-fig3a.dat.form$fig3a.dat))/fig3a.dat.form$rep)


fig.3a <- ggplot(fig3a.dat.form, aes(x = Ra, y = fig3a.dat, 
                                     ymin=fig3a.dat-err, ymax=fig3a.dat+err)) +
  geom_point() + 
  geom_ribbon(alpha=0.3, color = NA) +
  labs(x = "Relative Fitness of Mutant A over B",
       y = "Fixation Probability \n(1000 Simulations)",
       subtitle = "") + 
  annotate(geom="text", x=4, y=.1, label="N = 50 \nInitial Mutant A Population = 1",
                                 color="black")

fig.3a



#Graph 2 Fix Time vs Selective Advantage of Mutant
#N = 50
#Rep = 50
#A = 1
#Selective Advtange ranges from 1 to 5 by .2
#Standard Error of Mean
vars <- as.data.frame(expand.grid(
    N = 50,
    i = 1, 
    Ra = seq(1, 5 ,.2),
    Rb = 1, 
    rep = 1000))

fig3b.dat <- mapply(extinguisher.mut, i = vars$i, N = vars$N, 
                    Ra = vars$Ra, Rb = vars$Rb, 
                    rep = vars$rep, SIMPLIFY = TRUE)

fig3b.dat.form <- as.data.frame(cbind(vars, t(fig3b.dat)))
colnames(fig3b.dat.form) <- c("N", "i", "Ra", "Rb", 
                              "rep", "mean", "sd") 

fig3b.dat.form$sderr <- fig3b.dat.form$sd


fig.3b <- ggplot(fig3b.dat.form, aes(x = Ra, y = mean, 
                                     ymin=mean-sderr, ymax=mean+sderr)) +
  scale_y_continuous(trans='log10') +
  geom_point() + 
  geom_ribbon(alpha=0.3, color = NA) +
  labs(x = "Relative Fitness of Mutant A over B",
       y = "Fixation Time if \nMutant A Fixes (1000 Simulations)",
       subtitle = "") + 
  annotate(geom="text", x=4.2, y=2000, label="N = 50 \nInitial Mutant A Population = 1",
           color="black")

fig.3b

#Graph 3 Fix probability vs Population Size
#Fitness = 1.1
#A = 1
#Sims = 50
#Population Size from 0 to 100 by 5
#Standard Error of Proportion
vars <- as.data.frame(expand.grid(
  N = seq(5, 100 , 5),
  i = 1, 
  Ra = 1.1,
  Rb = 1, 
  rep = 2000))

fig3c.dat <- mapply(fixer.mut, i = vars$i, N = vars$N, 
                    Ra = vars$Ra, Rb = vars$Rb, 
                    rep = vars$rep, SIMPLIFY = TRUE)


fig3c.dat.form <- as.data.frame(cbind(vars, fig3c.dat))
fig3c.dat.form$err <- sqrt((fig3c.dat.form$fig3c.dat*(1-fig3c.dat.form$fig3c.dat))/fig3c.dat.form$rep)


fig.3c <- ggplot(fig3c.dat.form, aes(x = N, y = fig3c.dat, 
                                     ymin=fig3c.dat-err, ymax=fig3c.dat+err)) +
  geom_point() + 
  geom_ribbon(alpha=0.3, color = NA) +
  labs(x = "Total Population",
       y = "Fixation Probability \n(2000 Simulations)",
       subtitle = "") + 
  annotate(geom="text", x=75, y=.45, label="R = 1.1 \nInitial Mutant A Population = 1",
           color="black")

fig.3c

#Graph 3 Fix Time vs Population Size
#Fitness = 1.1
#A = 1
#Rep = 50
#Population Size from 0 to 100 by 5
#Standard Error of Mean
vars <- as.data.frame(expand.grid(
  N = seq(5, 100 , 5),
  i = 1, 
  Ra = 1.1,
  Rb = 1, 
  rep = 100))

fig3d.dat <- mapply(extinguisher.mut, i = vars$i, N = vars$N, 
                    Ra = vars$Ra, Rb = vars$Rb, 
                    rep = vars$rep, SIMPLIFY = TRUE)

fig3d.dat.form <- as.data.frame(cbind(vars, t(fig3d.dat)))
colnames(fig3d.dat.form) <- c("N", "i", "Ra", "Rb", 
                              "rep", "mean", "sd") 

fig3d.dat.form$sderr <- fig3d.dat.form$sd


fig.3d <- ggplot(fig3d.dat.form, aes(x = N, y = mean, 
                                     ymin=mean-sderr, ymax=mean+sderr)) +
  scale_y_continuous(trans='log10') +
  geom_point() + 
  geom_ribbon(alpha=0.3, color = NA) +
  labs(x = "Total Population",
       y = "Fixation Time if \nMutant A Fixes (100 Simulations)",
       subtitle = "") + 
  annotate(geom="text", x=75, y=300, label="R = 1.1 \nInitial Mutant A Population = 1",
           color="black")

fig.3d

#Correcting scales
fig.3a <- fig.3a + ylim(0, 1)
fig.3c <- fig.3c + ylim(0, 1)
fig.3b <- fig.3b + scale_y_continuous(trans='log10', limits = c(3, 9000))
fig.3d <- fig.3d + scale_y_continuous(trans='log10', limits = c(3, 9000))


png(filename = "Eash_frean_traulsen_fig_recreate.png",
    width = 14, height = 8, units = "in", res = 600)

grid.arrange(fig.3a, fig.3c, fig.3b, fig.3d, ncol = 2)

dev.off()

####Timing Code####
library(microbenchmark)
#Benchmarking A moran process that lasts for 5000 generations repeated 10 times
#Repeat Benchmark 10 times and take average
moran.process(i = 5, N = 10, Ra = 1.1, Rb = 1, rep = 10, step = 5000)


bench.N10 <- microbenchmark(moran.process(i = 1, N = 10, 
                              Ra = 1.1, Rb = 1, 
                              rep = 10, step = 5000), times = 10)
print(bench.N10)

bench.N100 <- microbenchmark(moran.process(i = 1, N = 100, 
                                          Ra = 1.1, Rb = 1, 
                                          rep = 10, step = 5000), times = 10)
print(bench.N100)

bench.N1000 <- microbenchmark(moran.process(i = 1, N = 1000, 
                                           Ra = 1.1, Rb = 1, 
                                           rep = 10, step = 5000), times = 10)
print(bench.N1000)

