
rm(list = ls())
set.seed(251093)
# Load working directory either at home or at CSS
if(substring(getwd(),0,1) == 'C'){
  setwd('C:/Users/fhk798/OneDrive/R Codes/') # PC
} else{
  setwd('/Users/ditlevkf/OneDrive/R Codes/') # Mac
}


library('dplyr')
library('ggplot2')
library('gridExtra')

load("Load your image here")

  
#Plot "real" profit at different prices. Add Profitmaximising prices for WLB and Bayesian Lasso.
Graph <- as.data.frame(cbind(price_vec,profit_real, prof, prof2))
Plot1 <- ggplot(Graph, aes(x= price_vec, y=profit_real), label=price_vec) + 
          geom_point() + 
          geom_point(data = Graph[which.max(prof),], aes(x=price_vec,y=profit_real), colour = "red", size = 3.5)+ 
          geom_text(data = Graph[which.max(prof),],aes(x=price_vec,y=profit_real), label = paste("WLB - Price =", Graph[which.max(prof),1],sep=" "), hjust=0, vjust=0) +        
          geom_point(data = Graph[which.max(prof2),], aes(x=price_vec,y=profit_real), colour = "blue", size = 3.5)+ 
          geom_text(data = Graph[which.max(prof2),],aes(x=price_vec,y=profit_real), label = paste("BLasso - Price =", Graph[which.max(prof2),1],sep=" "), hjust=0, vjust=-0) +        
          theme_minimal() + 
          labs(title = "Profit at different prices",y="Profit",x="Prices") 
Plot1

#Plot Densities of optimal Price
Graph2 <- data.frame(x=Price_Opt)
Graph3 <- data.frame(y=price_opt_lasso)
Plot2 <- ggplot() + 
            geom_density(data = Graph2, aes(x = x), 
               fill = "Red", color = "black", alpha = 0.7) + 
            geom_density(data = Graph3, aes(x = y),
               fill = "#56B4E9", color = "black", alpha = 0.7) +
            labs(title = "Density of Profit-maximising prices", x = "Profit-maximising Price", y="Density")
Plot2

#Plot Elasticities
ElasDATA <- as.data.frame(cbind(price_vec,elas_WLB))  
PlotElas <- ggplot(ElasDATA, aes(x=price_vec, y= elas_WLB)) + 
  geom_point(color = "blue") +
  theme_minimal() +
  labs(title = "Elasticity over 100 Bootstraps at different prices",y="MSE",x="Prices") 
PlotElas  

ElasDATA2 <- data.frame(cbind(price_vec,elas_BLasso))
PlotElas1 <- ggplot(ElasDATA2, aes(x=price_vec, y = elas_BLasso)) + 
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "Elasticity over 2000 draws at different prices",y="MSE",x="Prices") 
PlotElas1  

### Plot MSE  
MSEDATA <- as.data.frame(cbind(price_vec,MSE_WLB))  
PlotMSE <- ggplot(MSEDATA, aes(x=price_vec, y= MSE_WLB)) + 
  geom_point(color = "blue") +
  theme_minimal() +
  labs(title = "MSE over 100 Bootstraps at different prices",y="MSE",x="Prices") 


MSEDATA2 <- data.frame(cbind(price_vec,MSE_BLasso))
PlotMSE1 <- ggplot(MSEDATA2, aes(x=price_vec, y = MSE_BLasso)) + 
  geom_point(color = "red") +
  theme_minimal() +
  labs(title = "MSE over 2000 draws at different prices",y="MSE",x="Prices") 

grid.arrange(PlotElas,PlotElas1,PlotMSE,PlotMSE1)





