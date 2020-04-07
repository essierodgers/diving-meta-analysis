
#Functions

## Plots
plot_func <- function(data, x, y){
  plot(log(data[,y]) ~log(data[,x]), pch = 16, ylab = "log(y)", xlab = "log(x)") 
}
