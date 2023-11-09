
#
###
######
########
###########
#############
##############
###############
################
##################
####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
###################
################
##############
############
##########
#######
#####
###
##
##     cumulative plots of proportion of patients with gain or loss 
##
###
####
#####
######
#######
########
#########
##########
###########
#############
##############
###############
###############
###############
###############
###############





# x = binned Copy Number Variant (CNV) dataframe 
# y = dataset name -- Sample
# color1 = color for gain
# color2 = color for loss 

# function to plot 

cumulative_plot <- function(x,y,color1,color2) {
  length(unique(x$sample))
  cnbins_plot <- x[,c('sample','chr','startPos','ptotal')]
  head(cnbins_plot)
  #### rows with duplicate sample,chr,pos will have the ptotal averaged 
  cnbins_plot$key <- str_c(cnbins_plot$sample,'@',cnbins_plot$chr,'@',cnbins_plot$startPos)
  cnbins_plot_mean <- aggregate(ptotal ~ key, data=cnbins_plot, FUN=mean)
  #creating a category column
  cnbins_plot_mean$category <- 'neutral'
  cnbins_plot_mean[cnbins_plot_mean$ptotal > 2,]$category <- 'gain'
  cnbins_plot_mean[cnbins_plot_mean$ptotal < 2,]$category <- 'loss'
  head(cnbins_plot_mean)
  #splitting up the key by -
  cnbins_plot_mean2 <- cbind(cnbins_plot_mean,str_split_fixed(cnbins_plot_mean$key,'@',3))
  head(cnbins_plot_mean2)
  cnbins_plot_mean3 <- cnbins_plot_mean2[,c(4:6,2:3)]
  head(cnbins_plot_mean3)
  colnames(cnbins_plot_mean3) <- c('sample','chr','startPos','ptotal','category')
  str(cnbins_plot_mean3)
  head(cnbins_plot_mean3)
  length(unique(cnbins_plot_mean3$sample))
  ## converting to numeric
  cnbins_plot_mean3$chr <- as.numeric(cnbins_plot_mean3$chr)
  cnbins_plot_mean3$startPos <- as.numeric(cnbins_plot_mean3$startPos)
  head(cnbins_plot_mean3)
  cnbins_plot_mean4 <- cnbins_plot_mean3[,c('chr','startPos','category')]
  #getting rid of neutral rows 
  head(cnbins_plot_mean4)
  getwd()
  
  pdf(str_c(y,'_cumulative_plot_',Sys.Date(),'.pdf'),height=5,width=8)
  #creating layout of plots 
  #1 row x 22 columns (22 chromosomes)
  layout.matrix <- matrix(c(1:22), ncol = 22,nrow=1,byrow = TRUE)
  layout.matrix
  layout(mat = layout.matrix, # Heights of the two rows
         #using chr_para to specify widths relative to chromosome lengths
         widths = unique(chr_para[,c('chr','chrLength')])[1:22,'chrLength']) # Widths of the columns / we want them to be scaled / divide by chr length
  #bottom, left, top, right 
  par(oma=c(5,5,5,1),
      mar=c(1,0,0,0),
      xpd=NA)
  for (w in 1:22) {
    for (i in w:w) {
      print(i)
      a <- cnbins_plot_mean4[cnbins_plot_mean4$chr==i,]
      a
      b <- as.matrix(prop.table(table(a$category,a$startPos),2))
      b
      #str(b)
      c <- as.data.frame(b)
      head(c)
      colnames(c) <- c('category','pos','freq')
      head(c)
      #str(c)
      c$category <- as.character(c$category)
      c$pos <- as.numeric(as.character(c$pos))
      #c$freq <- as.numeric(c$freq)
      head(c)
      #str(c)
      #making the loss % negative
      c[c$category == 'loss',]$freq <- c[c$category == 'loss',]$freq*-1
      head(c)
      #get rid of neutral rows 
      d <- c[c$category != 'neutral',]
      head(d)
      
      if(nrow(d[d$category == 'gain',]) > 0) {
        plot(d[d$category == 'gain',]$pos,d[d$category == 'gain',]$freq,type='h',
             ylim = c(-1,1),
             col=color1,
             lwd=1,
             yaxt='n',
             xaxt='n',
             xlab='',
             ylab='',
             bty='n')
        
        
        
        par(new=TRUE)
      } else {
        plot(0,0,
             ylim = c(-1,1),
             col='white',
             lwd=1,
             yaxt='n',
             xaxt='n',
             xlab='',
             ylab='',
             bty='n')
        par(new=TRUE)
        }
      
      
      
      if(nrow(d[d$category == 'loss',]) > 0) {
        plot(d[d$category == 'loss',]$pos,d[d$category == 'loss',]$freq,type='h',
             ylim = c(-1,1),
             col=color2,
             lwd=1,
             yaxt='n',
             xaxt='n',
             xlab='',
             ylab='',
             bty='n')
      } else {
        plot(0,0,type='h',
             ylim = c(-1,1),
             col='white',
             lwd=1,
             yaxt='n',
             xaxt='n',
             xlab='',
             ylab='',
             bty='n')
        }
      
      #top and bottom solid lines
      #axis(3,labels = FALSE,lwd.ticks = 0,pos = 1,lwd = 2, xlim = c(0,unique(chr_para[,c('chr','chrLength')])[i,'chrLength'])) #top
      #axis(1,labels = FALSE,lwd.ticks = 0,pos = -1,lwd=2,xlim = c(0,unique(chr_para[,c('chr','chrLength')])[i,'chrLength']))  #bottom
      
      #x-axis chromsome label
      # odd numbers above the even
      if (i%%2 == 0) {
        title(xlab=paste(i),line=-0.4,cex.lab=1,xpd=TRUE)
        
      } else { title(xlab=paste(i),line=-1,cex.lab=1,xpd=TRUE) }
      
      
    }

    # adding the lefthand axis and y-axis label 
    if (w == 1) {
      axis(2,labels = c(1,0.5,0,0.5,1),at = c(-1,-0.5,0,0.5,1),line = 0, lwd = 2,cex.axis = 2, las=2)
      title(ylab = 'Proportion of Patients',line = 3.2,outer=TRUE,cex.lab = 2.2)
    }
    
    #right border of the plot - solid
    if (i == 22) {
      axis(4,labels = FALSE,lwd.ticks = 0,lwd = 2)
    } else axis(4,labels = FALSE,lwd.ticks = 0,lwd = 1,lty = 1) # middle dashed line between chromosomes  
    
    #drawing dashed line between p/q arms
    #xpos_arm <- chr_para[chr_para$chr==i & chr_para$arm == 'q',]$startPos
    #segments(x0=xpos_arm,y0=1,x1=xpos_arm,y1=-1,lty = 2)
    #xpos_end <- chr_para[chr_para$chr==i & chr_para$arm == 'q',]$endPos
    #drawing line segment on top & bottom of plot 
    #segments(x0=0,x1=xpos_end,y0=-1,y1=-1,lty = 1)
    #segments(x0=0,x1=xpos_end,y0=1,y1=1,lty = 1)
    #box(lty = 'dashed')
  }
  title(main = y,outer = TRUE,line = 1,cex.main = 3)
  title(xlab='Chromosome',line=1,cex.lab=2.2,outer = TRUE)
  #bottom, left, top, right 
  par(fig = c(0,1,0,1),oma = c(0, 3, 32,0), mar = c(0, 0, 0, 0), new = TRUE)
  legend('bottom',
         legend = c('Gain','Loss'),
         fill = c(color1,color2),
         xpd=TRUE,bty='n',
         border='white',
         horiz=TRUE,
         cex = 1.5)
  dev.off()
  
}



####################
# Import sample data set // 
####################

bb <- rio::import(bb,'1_mb_hg19_sample_binned_genome_cnv.txt')

# working directory to deposit the plots 

setwd("<INSERT DIRECTORY>")

####################
# Generate plot for sample data set 
####################
cumulative_plot(x=bb, y ='Sample Data Set',color1="forestgreen", color2="gold")
