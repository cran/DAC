#input
# data = vector of the data you have
# bench.min = lower bound for uniform benchmark prior
# bench.max = upper bound for uniform benchmark prior
# experts = matrix 1 row per expert, columns: mean, sd and skewness parameter
# lb.area = lower bound of the parameter space that is to be evaluated
# ub.area = upper bound of the parameter space that is to be evaluated
# step = the step size or precision by which the area of the parameterspace is evaluated
DAC_Uniform <- function(data, bench.min,bench.max, experts,lb.area,ub.area,step=0.01){


  ##################################################
  # Metropolis - Hastings normal data uniform prior#
  ##################################################

  MH_normal_unif <- function(iterations=10000, a1, b1, prop.sd, y,n.chains=4){

    b0.matrix <- matrix(NA,ncol=n.chains,nrow=iterations)
    sd.matrix <- matrix(NA,ncol=n.chains,nrow=iterations)
    for(i in 1:n.chains){

      ##EACH CHAIN

      #Create empty matrices
      b0 <- matrix(NA, iterations)
      sd <- matrix(NA, iterations)

      #Initial values
      b0[1,] <- runif(1,a1,b1) #get random starting value for chain within space of prior
      sd[1,] <- runif(1,0,(max(data)-min(data))) #get random starting value for chain sd
      # sd will not be smaller than 0 and not larger than distance between biggest and smallest value of data.

      #Creating the likelihood function
      likelihood <- function(b0, sd){
        singlelikelihoods <- dnorm(y, mean=b0, sd=sd, log=TRUE)
        sumll <- sum(singlelikelihoods)
        return(sumll)}

      #ALTERNATIEF: -N*(log(sqrt(2*pi)*sd))-((1/(2*sd^2))*(sum((y-b0)^2)))

      #Creating the prior function
      prior1 <- function(b0){
        b0.prior <- dunif(b0, min=a1, max=b1, log=T)
        return(b0.prior)
      }

      prior2 <- function(sd){
        sd.prior <- dunif(sd, min=0, max=(max(data)-min(data)), log=T)
        return(sd.prior)
      }

      #ALTERNATIEF: -log(b-a)

      #Creating the posterior function
      posterior1 <- function(b0, sd){
        return(likelihood(b0, sd) + prior1(b0))
      }

      posterior2 <- function(b0, sd){
        return(likelihood(b0, sd) + prior2(sd))
      }


      #Create for loop

      for(h in 2:iterations){

        #Sample beta0 with a random walk MH-step

        #Sample a value from the normal proposal density, with as mean the current value for b0 and the proposal sd
        proposal <- rnorm(1, mean=b0[h-1], sd=prop.sd)


        #Fill the proposal value in, in the posterior distribution
        post.prop <- posterior1(proposal, sd[h-1])

        #Fill the current value in, in the posterior distribution
        post.b0 <- (posterior1(b0[h-1], sd[h-1]))

        #Divide the first value by the second value, which is equal to subtracting the logs of the posteriors
        probab <- exp(post.prop-post.b0) #logposterior(proposal) - logposterior(b1[h-1]))


        #Sample a random value from a uniform distribution between 0 and 1
        u <- runif(1, min=0, max=1)

        #If the random value is lower than the difference between the logs, the proposal value will be taken, if not the current value of b1 is
        #taken again.
        if(u < probab){
          b0[h] <- proposal
        }else{
          b0[h] <- b0[h-1]
        }

        #Sample sd
        #Sample a proposal value
        proposal.sd <- sd[h-1]+runif(1, min=-0.1, max=0.1) #deze proposal heb ik uit Lynch

        #Fill the proposal value in, in the posterior distribution
        post.prop.sd <- (posterior2(b0[h], proposal.sd))

        #Fill the current value in, in the posterior distribution
        post.sd <- (posterior2(b0[h], sd[h-1]))

        #Divide the first value by the second value, which is equal to subtracting the logs of the posteriors
        probab2 <- exp(post.prop.sd-post.sd) #logposterior(proposal) - logposterior(b1[h-1]))

        #Sample a random value from a uniform distribution between 0 and 1
        u2 <- runif(1, min=0, max=1)

        #If the random value is lower than the difference between the logs, the proposal value will be taken, if not the current value of b1 is
        #taken again.
        if(u2 < probab2){
          sd[h] <- proposal.sd
        }else{
          sd[h] <- sd[h-1]
        }

      }

      b0.matrix[,i] <- b0
      sd.matrix[,i] <- sd
    } # end chain loop


    #burnin iterations
    burn <- (iterations/2)+1

    #Take half of the first chain --> obtain list with iterations after burnin in first chain
    b0.i <- b0.matrix[burn:iterations,]
    sd.i <- sd.matrix[burn:iterations,]

    #Obtain estimates from the combined chains
    #Means
    meanb0 <- mean(b0.i)
    meansd <- mean(sd.i)

    #Standard deviations
    sdb0 <- sd(b0.i)
    sd.sd <- sd(sd.i)

    #MC error
    MC.b0 <- sdb0/sqrt(iterations)
    MC.sd <- sd.sd/sqrt(iterations)


    #Acceptance rate
    acceptance1 <- 1-mean(duplicated(b0.i)) #Acceptance rate for the first chain

    acceptance1.sd <- 1-mean(duplicated(sd.i)) #Acceptance rate for the first chain


    acceptance.list <- data.frame(acceptance1, acceptance1.sd)
    colnames(acceptance.list)<- c("Acceptance b0", "Acceptance sd")

    frame <- data.frame(c(meanb0, meansd), c(sdb0, sd.sd), c(MC.b0, MC.sd), c(burn), c(iterations))

    row.names(frame) <- c("b0", "sd")
    colnames(frame)<- c("mean", "sd", "MC error", "start", "sample")

    output <- list(frame, acceptance.list,b0.i,sd.i)
    return(output)


  }

  out <- MH_normal_unif(a1=bench.min,b1=bench.max,y=data,prop.sd=(sd(data)/1.5))
  # get posterior from data and benchmark prior. note that the prop.sd is arbitrary and if the function fails this is a likely cause
  posterior <- unlist(out[1])[c(1,3)]


 ## KL divergences
  # get the hyper parameters for the distribution of the mean
  x.axis <- seq(lb.area,ub.area,length.out = ((ub.area-lb.area)/step))
  post <- dnorm(x.axis,posterior[1],posterior[2])


  #KL benchmark and posterior
  benchmark <- dunif(x.axis,bench.min,bench.max)

  #Creating a matrix of the posterior and benchprior for the use of the KL D function
  matrix1 <- cbind(post, benchmark)

  #Kullback-Leibler Divergence for posterior and benchmark prior
  KL.bench <- KLdiv(matrix1, method= c("continuous"), eps=10^-20)[1,2]



  experts.densities <- matrix(ncol=nrow(experts),nrow=length(x.axis))
  #create a storage matrix to store densities for each expert in a separate column
  #length rows equals the length of the x axis that is to be evaluated

  for(i in 1:nrow(experts)){
  experts.densities[,i] <- dsnorm(x.axis,experts[i,1],experts[i,2],experts[i,3])
  #get from input for each expert (row) the values for the mean (col1) sd(col2) and skewness(col3)
  }
  # for each epxert store their density in a column, density is based on skewed normal distribution
  # input for the parametrization of the expert densities is provided by user in the function.

    #KL benchmark and experts

  #Creating a matrix of the posterior and experts for the use of the KL D function
  matrix2 <- cbind(post, experts.densities)



  KL.experts <- matrix(NA,ncol=ncol(experts.densities))
  for(ii in 1:ncol(experts.densities)){
  KL.experts[,ii] <- KLdiv(cbind(matrix2[,1],matrix2[,ii+1]),method= c("continuous"), eps=10^-20)[1,2]
  # get KL diverence between posterior (column 1 of matrix 2) for each expert (expert 1 = column 2 in matrx 2 etc.)
  }

  DAC <- KL.experts / KL.bench
  #calculate DAC values
  names.DAC <- c()
  for(q in 1:length(DAC)){
    names.DAC[q] <-  paste("DAC Expert",q,sep = " ")
  }
  colnames(DAC) <- names.DAC
  print(DAC)
  #print DAC values

  #prepare output for output list
  #form a matrix with KL divergence values for each prior including the benchmark
  KL_divergences <- matrix(c(KL.experts,KL.bench),nrow=length(KL.experts)+1)
  colnames(KL_divergences) <- c("KL Divergence")
  names <- c()
  for(qq in 1:length(KL.experts)){
    names[qq] <-  paste("Expert",qq,sep = " ")
  }
  rownames(KL_divergences) <- c(names,"Benchmark")
  #name the output columns and rows

  #clean up the posterior output
  posterior.out <- matrix(as.numeric(posterior),nrow=1)
  colnames(posterior.out) <- c("Mean","SD")

  #create output list to store all relevant information so this can be obtained.
  out <- list(DAC.scores=DAC,posterior.distribution = posterior.out,KL.divergences = KL_divergences)
  return(out)
  #return list with relevant information
}


