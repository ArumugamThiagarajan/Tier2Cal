library(foreach)

get.values <- function(RunMap){
out.df<- foreach(i=1:nrow(RunMap),.combine='rbind') %do% {
  mu=RunMap[i,]$init.soc.mu
  sd=RunMap[i,]$init.soc.upper-RunMap[i,]$init.soc.mu
  values = rnorm(100,mean=mu,sd=sd)
  val.df <- t(as.data.frame(values))
  for(j in 1:ncol(val.df)){
    names(val.df)[j]<- paste0('init.iter.',j)
  }
  val.df$RunIndex=RunMap[i,]$RunIndex
  return(val.df)
} 
rownames(out.df) <- NULL
out.df <- data.frame(out.df)
column_index <- which(names(out.df) == "RunIndex")
out.df <- out.df[, c(column_index, 1:(column_index-1))]

return(out.df)
}
RunMap<-out.table
init.soc.df<- get.values(RunMap)

init.soc.df <- apply(init.soc.df,2,as.character)

write.csv(init.soc.df,
          paste0(main.dir,
                 "//data//prepared data//RunMap_iterative_initSOC_sept25.csv"),
          row.names = F)

