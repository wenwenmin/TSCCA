library("ggplot2")

load("Result_ku10_kv100_kw20_J50.RData")

objs = out$d_iter
lines.num = length(objs)
figs = list()
for(i in 1:lines.num){
  obj = objs[[i]]
  Iterations = 1:length(obj)
  dat_temp = as.data.frame(cbind(Iterations,obj)) 
  figs[[i]] = ggplot(dat_temp, aes(x=Iterations, y=obj)) + geom_line() + geom_point() + 
    xlab("Iterations") + ggtitle(NULL) + ylab("Obj") + 
    annotate("text",  x = -Inf, y = Inf, label = paste(" M", i, sep=""), 
             vjust=1, hjust=0) + theme_classic()
}
# please see what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot
fig = plot_grid(figs[[1]],figs[[2]],figs[[3]],figs[[4]],figs[[5]],
                figs[[6]],figs[[7]],figs[[8]],figs[[9]],figs[[10]], 
                figs[[11]],figs[[12]],figs[[13]],figs[[14]],figs[[15]],
                figs[[16]],figs[[17]],figs[[18]],figs[[19]],figs[[20]],
                nrow =4,ncol=5)

ggsave("Fun1_Obj.png",width = 10, height = 7.5) 


