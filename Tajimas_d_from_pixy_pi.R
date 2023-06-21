###############################################################################################
###############################################################################################
##############    Calculate Tajima's D using the output of pixy    ############################
###############################################################################################
###############################################################################################


### Tajima's D: D = (pi - S/a1)/sqrt(V)
## info from calculate_tajimas_D.pdf
# pi (pi value): pi = avg_pi_pixy *  no_sites_pixy
# S (number of segregating sites): "need do be calculated from the VCF"
# n (number of samples): 

## here some variables defined as functions:
a1 <- function(n){return(sum(sapply(1:(n-1), function(i){return(1/i)})))}
a2 <- function(n){return(sum(sapply(1:(n-1), function(i){return(1/(i^2))})))}
b1 <- function(n){return((n+1)/(3*(n-1)))}
b2 <- function(n){return((2*(n^2+n+3))/(9*n*(n-1)))}
c1 <- function(n){return(b1(n) - (1/a1(n)))}
c2 <- function(n){return(b2(n) - ((n+2)/(a1(n)*n)) + (a2(n)/(a1(n)^2)  ))}
e1 <- function(n){return(c1(n)/a1(n))}
e2 <- function(n){return(c2(n)/(a1(n)^2 + a2(n)))}

## here the function for Tajima's D 
Tajimas_D <- function(pi,S,n){return((pi - (S/a1(n)))/sqrt(e1(n)*S + e2(n)*S*(S-1)) )}

## defining some fixed variables:
n_A_Italia <- 66*2
n_A_East <- 94*2
n_B_West <- 72*2
n_B_East <- 73*2
n_C <- 27*2



## read in data:
## bed files:
Bdistachyon_314_v3.2.gene.Bd1 <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdistachyon_314_v3.2.gene.Bd1.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdistachyon_314_v3.2.gene.Bd2 <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdistachyon_314_v3.2.gene.Bd2.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdistachyon_314_v3.2.gene.Bd3 <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdistachyon_314_v3.2.gene.Bd3.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdistachyon_314_v3.2.gene.Bd4 <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdistachyon_314_v3.2.gene.Bd4.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdistachyon_314_v3.2.gene.Bd5 <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdistachyon_314_v3.2.gene.Bd5.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

## Pixy results
Bdis_332_Bd1_genes_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_genes_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Bdis_332_Bd2_genes_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_genes_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Bdis_332_Bd3_genes_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_genes_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Bdis_332_Bd4_genes_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_genes_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Bdis_332_Bd5_genes_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_genes_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## number of segregating
Bdis_332_Bd1_pixy_A_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_A_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_A_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_A_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_A_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_A_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_A_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_A_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_A_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_A_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

Bdis_332_Bd1_pixy_A_Italia_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_A_Italia_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_A_Italia_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_A_Italia_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_A_Italia_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_A_Italia_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_A_Italia_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_A_Italia_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_A_Italia_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_A_Italia_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

Bdis_332_Bd1_pixy_B_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_B_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_B_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_B_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_B_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_B_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_B_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_B_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_B_East_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_B_East_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

Bdis_332_Bd1_pixy_B_West_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_B_West_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_B_West_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_B_West_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_B_West_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_B_West_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_B_West_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_B_West_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_B_West_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_B_West_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

Bdis_332_Bd1_pixy_C_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_C_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_C_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_C_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_C_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_C_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_C_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_C_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_C_number_of_segregating_sites <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_C_number_of_segregating_sites.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

## calculate Tajima's D:
Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D <- cbind(Bdis_332_Bd1_genes_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D)[10] <- "tajimas_D"
Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[,10] <- unlist(sapply(1:length(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,4]
  my_pi <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,5] * Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D)[11] <- "number_of_segregating_site"
Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[,11] <- unlist(sapply(1:length(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))

Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D <- cbind(Bdis_332_Bd2_genes_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D)[10] <- "tajimas_D"
Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[,10] <- unlist(sapply(1:length(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,4]
  my_pi <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,5] * Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D)[11] <- "number_of_segregating_site"
Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[,11] <- unlist(sapply(1:length(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))

Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D <- cbind(Bdis_332_Bd3_genes_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D)[10] <- "tajimas_D"
Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[,10] <- unlist(sapply(1:length(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,4]
  my_pi <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,5] * Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D)[11] <- "number_of_segregating_site"
Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[,11] <- unlist(sapply(1:length(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))

Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D <- cbind(Bdis_332_Bd4_genes_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D)[10] <- "tajimas_D"
Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[,10] <- unlist(sapply(1:length(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,4]
  my_pi <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,5] * Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D)[11] <- "number_of_segregating_site"
Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[,11] <- unlist(sapply(1:length(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))

Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D <- cbind(Bdis_332_Bd5_genes_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D)[10] <- "tajimas_D"
Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[,10] <- unlist(sapply(1:length(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,4]
  my_pi <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,5] * Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D)[11] <- "number_of_segregating_site"
Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[,11] <- unlist(sapply(1:length(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[,10]), function(x){
  my_clade <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,1]
  my_scaffold <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,2]
  my_start <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,3] -1
  my_end <- Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))


quantile(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D[,10], na.rm = TRUE)
quantile(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D[,10], na.rm = TRUE)
quantile(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D[,10], na.rm = TRUE)
quantile(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D[,10], na.rm = TRUE)
quantile(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D[,10], na.rm = TRUE)


## write out data:
write.table(rbind(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D, Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D, Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D, Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D, Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D), file = "/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_genes_pi_pixy_and_tajimas_D.txt", quote = FALSE)






#######################################
#######################################
###### genome wide
#######################################
#######################################


## read in data:
## bed files:
Bd1_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bd1_5000_windows.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bd2_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bd2_5000_windows.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bd3_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bd3_5000_windows.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bd4_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bd4_5000_windows.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bd5_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bd5_5000_windows.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

## Pixy results
Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## number of segregating
Bdis_332_Bd1_pixy_A_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_A_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_A_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_A_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_A_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_A_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_A_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_A_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_A_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_A_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

Bdis_332_Bd1_pixy_A_Italia_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_A_Italia_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_A_Italia_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_A_Italia_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_A_Italia_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_A_Italia_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_A_Italia_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_A_Italia_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_A_Italia_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_A_Italia_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

Bdis_332_Bd1_pixy_B_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_B_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_B_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_B_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_B_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_B_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_B_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_B_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_B_East_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_B_East_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

Bdis_332_Bd1_pixy_B_West_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_B_West_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_B_West_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_B_West_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_B_West_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_B_West_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_B_West_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_B_West_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_B_West_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_B_West_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

Bdis_332_Bd1_pixy_C_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd1_pixy_C_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd2_pixy_C_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd2_pixy_C_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd3_pixy_C_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd3_pixy_C_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd4_pixy_C_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd4_pixy_C_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
Bdis_332_Bd5_pixy_C_number_of_segregating_sites_5000_windows <- read.table("/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_Bd5_pixy_C_number_of_segregating_sites_5000_windows.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

## calculate Tajima's D:
Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D_5000_windows <- cbind(Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D_5000_windows)[10] <- "tajimas_D"
Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D_5000_windows[,10] <- unlist(sapply(1:length(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  my_pi <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,5] * Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D_5000_windows)[11] <- "number_of_segregating_site"
Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D_5000_windows[,11] <- unlist(sapply(1:length(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd1_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))

Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D_5000_windows <- cbind(Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D_5000_windows)[10] <- "tajimas_D"
Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D_5000_windows[,10] <- unlist(sapply(1:length(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  my_pi <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,5] * Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D_5000_windows)[11] <- "number_of_segregating_site"
Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D_5000_windows[,11] <- unlist(sapply(1:length(Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd2_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))

Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D_5000_windows <- cbind(Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D_5000_windows)[10] <- "tajimas_D"
Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D_5000_windows[,10] <- unlist(sapply(1:length(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  my_pi <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,5] * Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D_5000_windows)[11] <- "number_of_segregating_site"
Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D_5000_windows[,11] <- unlist(sapply(1:length(Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd3_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))

Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D_5000_windows <- cbind(Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D_5000_windows)[10] <- "tajimas_D"
Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D_5000_windows[,10] <- unlist(sapply(1:length(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  my_pi <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,5] * Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D_5000_windows)[11] <- "number_of_segregating_site"
Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D_5000_windows[,11] <- unlist(sapply(1:length(Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd4_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))

Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows <- cbind(Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi, NA, NA)
names(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows)[10] <- "tajimas_D"
Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows[,10] <- unlist(sapply(1:length(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  my_pi <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,5] * Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,6]
  my_n <- get(paste0("n_",my_clade))
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  if(length(my_number_of_segregating_site) != 1){
    return("not_found")
  } else {
    return(Tajimas_D(pi = my_pi, S = my_number_of_segregating_site, n = my_n))
  }
}))
names(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows)[11] <- "number_of_segregating_site"
Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows[,11] <- unlist(sapply(1:length(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows[,10]), function(x){
  my_clade <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,1]
  my_scaffold <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,2]
  my_start <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,3]
  my_end <- Bdis_332_Bd5_genome_wide_5000_windows_pi_pixy_output_pi[x,4]
  
  my_number_of_segregating_site <- get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[
    get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,2] == my_start &
      get(paste0("Bdis_332_", my_scaffold, "_pixy_", my_clade, "_number_of_segregating_sites_5000_windows"))[,3] == my_end,
    4]
  
  return(my_number_of_segregating_site)
  
}))


quantile(Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows[,11], na.rm = TRUE)


## write out data:
write.table(rbind(Bdis_332_Bd1_genes_pi_pixy_and_tajimas_D_5000_windows, Bdis_332_Bd2_genes_pi_pixy_and_tajimas_D_5000_windows, Bdis_332_Bd3_genes_pi_pixy_and_tajimas_D_5000_windows, Bdis_332_Bd4_genes_pi_pixy_and_tajimas_D_5000_windows, Bdis_332_Bd5_genes_pi_pixy_and_tajimas_D_5000_windows), file = "/Users/roberthorvath/Desktop/Projects/X) Slim_for_Nikos/4) Results/Pixy/Bdis_332_genome_wide_5000_windows_pi_pixy_and_tajimas_D.txt", quote = FALSE)






