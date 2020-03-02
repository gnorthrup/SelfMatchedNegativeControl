setwd("~/Desktop/berkeley/Lewnard Lab/Non-Parametric Estimation NAI/Rotavirus") 

load("india.diars.Rdata");load("mexico.diars.Rdata");load("india.infects.Rdata");load("mexico.infects.Rdata")

#########################
###get indian A B C D ###
#########################

# duprows <- rownames(india.diars) %in% rownames(india.infects)
# total <- rbind(subset(india.diars,select=-c(modsev)), subset(india.infects,select=-c(num))[!duprows,])

total <- unique(rbind(subset(india.diars,select=-c(modsev)),subset(india.infects,select=-c(num))))
india.total <- total[order(total$id, total$age),]

india.people <- unique(india.total["id"])

india.two.rota <- c()

for (person in 1:nrow(india.people)){
  #print(india.people[person,"id"])
  if(sum(subset(india.total,id==india.people[person,"id"],select=rv))<2){
    next
  } else {
    india.two.rota <- c(india.two.rota,india.people[person,"id"])
  }
}
india.times.list <- c()
person.num <- 1
for (person in india.two.rota){
  cases <- subset(india.total,id==person)
  if (length(cases[,2])<5){
    next
  }
  if (length(rle(cases[,2])$values)<4){
    if (length(rle(cases[,2])$values) < 3 | rle(cases[,2])$values[2]==0){
      next
    }
  }
  if (cases[1,2]==1){
    next
  }
  time.zero <- cases[1,3]
  num.rota <- 0
  A.empty <- TRUE
  B.empty <- TRUE
  C.empty <- TRUE
  D.empty <- TRUE
  for (i in 2:nrow(cases)){
    if (A.empty){
      if (cases[i,2]==1){
        A <- cases[i,3]-time.zero
        new.time <- cases[i,3]
        new.point <- i
        A.empty <- FALSE
      }
    }
    if (B.empty){
      if (cases[i,2]==0){
        B <- cases[i,3]-time.zero
        new.time <- cases[i,3]
        new.point <- i
        B.empty <- FALSE
      }
    }
    if (!A.empty & !B.empty){
      break
    }
  }
  next.start <- new.point+1
  if(next.start>=nrow(cases)){
    next
  }
  for (i in next.start:nrow(cases)){
    if (C.empty){
      if (cases[i,2]==1){
        C <- cases[i,3]-new.time
        C.empty <- FALSE
      }
    }
    if (D.empty){
      if (cases[i,2]==0){
        D <- cases[i,3]-new.time
        D.empty <- FALSE
      }
    }
    if (!C.empty & !D.empty){
      break
    }
  }
  person.vec <- c("id"=person,"A"=A,"B"=B,"C"=C,"D"=D)
  india.times.list[[person.num]] <- person.vec
  person.num <- person.num + 1
}

##################
### Get Mexico ###
##################

total <- unique(rbind(subset(mexico.diars),subset(mexico.infects,select=-c(num))))
mexico.total <- total[order(total$id, total$age),]

mexico.people <- unique(mexico.total["id"])

mexico.two.rota <- c()

for (person in 1:nrow(mexico.people)){
  #print(india.people[person,"id"])
  if(sum(subset(mexico.total,id==mexico.people[person,"id"],select=rv))<2){
    next
  } else {
    mexico.two.rota <- c(mexico.two.rota,mexico.people[person,"id"])
  }
}
mexico.times.list <- c()
person.num <- 1
for (person in mexico.two.rota){
  cases <- subset(mexico.total,id==person)
  if (length(cases[,2])<5){
    next
  }
  if (length(rle(cases[,2])$values)<4){
    if (length(rle(cases[,2])$values) < 3 | rle(cases[,2])$values[2]==0){
      next
    }
  }
  if (cases[1,2]==1){
    next
  }
  time.zero <- cases[1,3]
  num.rota <- 0
  A.empty <- TRUE
  B.empty <- TRUE
  C.empty <- TRUE
  D.empty <- TRUE
  for (i in 2:nrow(cases)){
    if (A.empty){
      if (cases[i,2]==1){
        A <- cases[i,3]-time.zero
        new.time <- cases[i,3]
        new.point <- i
        A.empty <- FALSE
      }
    }
    if (B.empty){
      if (cases[i,2]==0){
        B <- cases[i,3]-time.zero
        new.time <- cases[i,3]
        new.point <- i
        B.empty <- FALSE
      }
    }
    if (!A.empty & !B.empty){
      break
    }
  }
  next.start <- new.point+1
  if(next.start>=nrow(cases)){
    next
  }
  for (i in next.start:nrow(cases)){
    if (C.empty){
      if (cases[i,2]==1){
        C <- cases[i,3]-new.time
        C.empty <- FALSE
      }
    }
    if (D.empty){
      if (cases[i,2]==0){
        D <- cases[i,3]-new.time
        D.empty <- FALSE
      }
    }
    if (!C.empty & !D.empty){
      break
    }
  }
  person.vec <- c("id"=person,"A"=A,"B"=B,"C"=C,"D"=D)
  mexico.times.list[[person.num]] <- person.vec
  person.num <- person.num + 1
}

#################
### Get Theta ###
#################

times.list <- c(india.times.list,mexico.times.list)

num_sum <- 0
den_sum <- 0

for (i in 1:length(times.list)){
  if (times.list[[i]]['A']<times.list[[i]]['B']){
    Gi <- 1
    Hi <- 0
  } else {
    Gi <- 0
    Hi <- 1
  }
  
  if (times.list[[i]]['C']<times.list[[i]]['D']){
    Ei <- 1
    Fi <- 0
  } else {
    Ei <- 0
    Fi <- 1
  }
  
  num_sum <- num_sum + Ei*Hi
  den_sum <- den_sum + Fi*Gi
}

theta <- num_sum/den_sum

real.est <- theta

print(theta)

#################
### bootstrap ###
#################

numstrap <- 10000

sample.size <- length(times.list)

est.list <- c()

for (i in 1:numstrap){
  samples <- sample(1:sample.size,sample.size,replace=TRUE)
  
  num_sum <- 0
  den_sum <- 0
  
  for (j in samples){
    if (times.list[[j]]['A']<times.list[[j]]['B']){
      Gi <- 1
      Hi <- 0
    } else {
      Gi <- 0
      Hi <- 1
    }
    
    if (times.list[[j]]['C']<times.list[[j]]['D']){
      Ei <- 1
      Fi <- 0
    } else {
      Ei <- 0
      Fi <- 1
    }
    
    num_sum <- num_sum + Ei*Hi
    den_sum <- den_sum + Fi*Gi
  }
  
  theta <- num_sum/den_sum
  
  est.list <- c(est.list,theta)
  
}

upper <- quantile(est.list,.975)
lower <- quantile(est.list,.025)

print(quantile(est.list,.975))
print(quantile(est.list,.025))

baseProtect <- c('lower'=lower,'est'=real.est,'upper'=upper)

save(baseProtect,file='baseProtect.Rdata')
