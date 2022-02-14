##cited from https://jef.works/blog/2017/09/21/connected-barplot-for-series-data-visualization/

connectedBarplot <- function(dat,color=rainbow(nrow(dat)), space=1, alpha=0.4, ...) {  
  b <- barplot(dat, col=color, space = space, ...)                     
  
  for (i in seq_len(ncol(dat) - 1)) {     
    lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## bottom line       
    
    for (j in seq_len(nrow(dat))) {     
      if (j == 1) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))                       
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(0, dat[j,i], dat[j,i+1], 0),               
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j == 2) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),                       
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j > 2) {                    
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], colSums(dat[1:(j-1),])[i+1]),              
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
    }          
  }              
}     