## Please do not use this as an example for good graphing practices.
## While the end product is pretty, the code is very hastily put together.
## Because of this, at times it is redundant and inconsistent.
## I was strapped for time!

library('Rmisc')
library('phyloseq')
library('ggplot2')
library('reshape2')
library('grid')
library('gridExtra')
library('gtable')

source('permutationTest.r')

## Create "plus-or-minus" operator
'%+-%' = function(x,y) return(c(x+y,x-y))

se = function(x) {
## function se
## parameters: <vector> x: a numeric vector
## returns: <double> the standard deviation of x
    return(sd(x)/sqrt(length(x)))
}

getRelevantStats = function(x) {
## function getRelevantStats
## parameters: <vector> x: a numeric vector
## returns: <vector> the mean, standard deviation, standard error, and 95%
##          confidence interval of x
    ## Put together vector of mean, sd, se, and 95 percent confidence interval
    out = c(mean(x), sd(x), se(x), mean(x) %+-% (2*se(x)))
    
    ## Assign names to out vector
    names(out) = c('Mean', 'Standard Deviation', 'Standard Error',
                   '95% Confidence Upper', '95% Confidence Lower')
    
    return(out)
}

matrixify_list = function(x) {
    mat = matrix(nrow=0,ncol=length(ci[[1]]))
    for(i in 1:length(x)) {
        mat = rbind(mat, x[[i]])
    }
    rownames(mat) = names(x)
    return(mat)
}

addConfidenceToStripChart = function(x,groups, p) {
## function addConfidenceToStripChart
## parameters: <vector> x: vector of values plotted; <vector> y: groupings of x
## outputs: Prints mean dot and confidence interval bars to current plot.
##          Prints 95% confidence interval via getRelevantStats()
    ## Coerce groups to a factor variable
    group_factor = as.factor(groups)
    
    ## Create confidence interval
    ci = tapply(x,group_factor,getRelevantStats)
    
    ## Plot the means
    #points(sapply(ci, function(ci) return(ci[1])),cex=2,pch=20)
    
    ## Get cartesian coords for confidence lines
    x0 = c(0,30) #as.numeric(unique(group_factor))
    yu = sapply(ci, function(ci) return(ci[4])) #  Upper
    yl = sapply(ci, function(ci) return(ci[5])) #  Lower
    y0 = c(yu,yl)
    
    print(yl)
    
    ## Draw 95% confidence interval and print confidence interval to screen
    #segments(x0 - 0.05, y0, x0 + 0.05, y0, lwd=2)
    #segments(x0, yu, x0, yl, lwd=2)
    
    for(i in 1:length(unique(group_factor)))
    {
        
        p = p + geom_segment(x=x0[i]-1, xend=x0[i]+1, y=yu[i], yend=yu[i]) +
        geom_segment(x=x0[i]-1, xend=x0[i]+1, y=yl[i], yend=yl[i]) +
        geom_segment(x=x0[i], xend=x0[i], y=yu[i], yend=yl[i]) + 
        geom_segment(x=x0[i]-2, xend=x0[i]+2, y=sapply(ci, function(ci) return(ci[1]))[i], yend=sapply(ci, function(ci) return(ci[1]))[i],size=1.5)
    }
    return(p)
}

vectors = list.files('../exported-metrics','*vector.tsv')

titles = c("Observed SVs",
           "Pielou's E",
           "Shannon",
           "Faith's PD")
yticks = data.frame(obs = c(100,110,120,130,140,150,160,170,180,190,200),
              pielou = c(0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79, 0.80,0.81,0.82),
              shannon = c(5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2),
              faith = c(8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13))
columnnames = c('obs','pielou','shannon','faith')
              
f_physeq = prune_samples(rownames(sample_data(physeq))[sample_data(physeq)$trmnt_type_s == "CamSa"],physeq)
f_physeq = prune_taxa(taxa_sums(f_physeq) > 0, f_physeq)
plots = list()

d.f = data.frame(groups=sample_data(f_physeq)[,"trmnt_day_s"])

for (i in 1:length(vectors))
{
    fname = paste0('../exported-metrics/',vectors[i])
    tab = read.table(fname,row.names=1,header=TRUE)
    tab = data.frame(index=tab[match(sample_names(f_physeq),rownames(tab)),])
    colnames(tab) = columnnames[i]
    d.f = cbind(d.f,tab)
}

d.f = na.omit(d.f)
tests = list(kt=vector(mode='double'))

sink('alphaDiversityStats.txt')
for(i in 2:ncol(d.f))
{
    print(colnames(d.f)[i])
    ## Shapiro-Wilk test
    sh = tapply(d.f[,i], d.f$trmnt_day_s, shapiro.test)
    print(sh)
    ## Kruskal-Wallace test
    kt = kruskal.test(d.f[,i], d.f$trmnt_day_s)
    print(kt)
    ## Welch's unpaired t-test
    tt = t.test(d.f[,i]~d.f$trmnt_day_s)
    print(tt)
    tests$kt = c(tests$kt,kt$p.value)
    permt = permutation.test(d.f[,i], d.f$trmnt_day_s,nperm=5000)
    print(permt)
    print("==============================")
}
sink()
    
sig = data.frame(variable=colnames(d.f)[2:ncol(d.f)],tests$kt)

#print(sig)
#ylims = c(yticks[[i]][1],yticks[[i]][length(yticks[[i]])])

means = matrix(nrow=0, ncol=ncol(d.f)-1)
se.lower = matrix(nrow=0, ncol=ncol(d.f)-1)
se.upper = matrix(nrow=0, ncol=ncol(d.f)-1)
for(i in 1:length(unique(d.f$trmnt_day_s)))
{
    stats = apply(d.f[d.f$trmnt_day_s == unique(d.f$trmnt_day_s)[i],2:ncol(d.f)],2,getRelevantStats)
    means = rbind(means,stats[1,])
    se.upper = rbind(se.upper, stats[5,])
    se.lower = rbind(se.lower, stats[4,])
    rownames(means)[i] = unique(d.f$trmnt_day_s)[i]
    rownames(se.lower)[i] = unique(d.f$trmnt_day_s)[i]
    rownames(se.upper)[i] = unique(d.f$trmnt_day_s)[i]
}

means = data.frame(t(means),check.names=FALSE)
se.upper = data.frame(t(se.upper),check.names=FALSE)
se.lower = data.frame(t(se.lower),check.names=FALSE)

means$variable = rownames(means)
se.lower$variable = rownames(se.lower)
se.upper$variable = rownames(se.upper)

means = melt(means)
means[,2] = as.integer(as.integer(as.character(means[,2])))
colnames(means)[2] = 'trmnt_day_s'
se.lower = melt(se.lower)
se.lower[,2] = as.integer(as.character(se.lower[,2]))
colnames(se.lower)[2] = 'trmnt_day_s'

se.upper = melt(se.upper)
se.upper[,2] = as.integer(as.character(se.upper[,2]))
colnames(se.upper)[2] = 'trmnt_day_s'

errorbar = cbind(means, se.lower[,3], se.upper[,3])
colnames(errorbar) = c('variable','trmnt_day_s','mean','se.lower','se.upper')

errorbar$variable = factor(errorbar$variable,levels=columnnames)

mdf = melt(d.f)
mmdf = mdf[mdf$variable != "trmnt_day_s",]
mmdf$trmnt_day_s = d.f$trmnt_day_s
mmdf$variable = factor(mmdf$variable,levels=columnnames)

ylims = data.frame(ymin=t(yticks[1,]), ymax=t(yticks[nrow(yticks),]), xmin=-15, xmax=45, variable=factor(columnnames, levels=columnnames))

colnames(ylims)[1:2] = c('ymin','ymax')

breaks = melt(yticks)
breaks$variable = factor(breaks$variable,levels=columnnames)

sig$coords = melt(yticks[nrow(yticks),])[,"value"]
sig$label = paste0('p=',round(sig$tests.kt,3))

samlabs = as_labeller(c(`obs` = "Observed SVs", `pielou` = "Pielou's E", `shannon` = "Shannon", `faith` = "Faith's PD"))
plots = lapply(columnnames, function(f) {
p = ggplot() +
facet_wrap(~variable, scales='free',nrow=1, labeller=samlabs) +
geom_errorbar(data=errorbar[errorbar$variable==f,], aes(x=trmnt_day_s, ymin=se.lower, ymax=se.upper), width=6, size=1) +
geom_errorbar(data=errorbar[errorbar$variable==f,], aes(x=trmnt_day_s, ymin=mean, ymax=mean), width=8, size=2) +
geom_text(data=sig[sig$variable==f,], aes(x=-15, y=coords, label=label), size=(5/14)*18, family='sans',hjust=-0.1,vjust=1.5) +
geom_blank(data=ylims[ylims$variable==f,], aes(ymin=ymin, ymax=ymax,xmin=-15,xmax=45)) +
geom_point(data=mmdf[mmdf$variable==f,], aes(x=trmnt_day_s, y=value, fill=as.factor(trmnt_day_s)), size=4, pch=21) +
scale_fill_manual(values=c('#FEBF00','#006FBF')) +
theme_bw() +
scale_x_continuous(labels=c('Day 0','Day 30'),breaks=c(0,30), expand = c(0, 0)) +
scale_y_continuous(breaks=breaks[breaks$variable==f,"value"],expand = c(0,1e-8)) +
theme(legend.position = 'none',
      axis.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color='gray'),
      axis.line = element_line(color='gray'),
      axis.ticks = element_line(color='gray'),
      axis.ticks.length = unit(0.3, "lines"),
      strip.background = element_blank(), #element_rect(color='black',fill='transparent'),
      panel.border = element_blank(),
      plot.background = element_rect(color='black', fill=NA),
      text = element_text(color='black', size=24, family='sans'),
      axis.text = element_text(color='black', size=18),
      strip.text = element_text(color='black', size=18),
      panel.spacing = unit(2,'lines'),
      # top,
      plot.margin = margin(5.5,44,5.5,22)
     )
    #lw = 3
    #ci = tapply(d.f$index,d.f$trmnt_day_s,getRelevantStats)
    #ci = matrixify_list(ci)
    #ci.size = max(as.vector(ci[,4] - ci[,5]))
    #window.size = ci.size * 3/2
    
    #window.ymin = max(0,round(min(as.vector(ci[,4] - ci[,5])) - (window.size/2),digits=-1*(floor(log10(window.size))+1)))
    #window.ymax = round(max(as.vector(ci[,4] - ci[,5])) + (window.size/2))
    #print(window.size)
    #print(c(window.ymin,window.ymax))

    # p = ggplot(data=d.f)
    # p = addConfidenceToStripChart(x=d.f$index, group=d.f$trmnt_day_s,p=p) +
    # geom_point(pch=21, size=4, aes(fill=as.factor(trmnt_day_s),
                               # x=abs(as.integer(trmnt_day_s)),
                               # y=index)) +
    # #geom_segment(x=0-lw, xend=0+lw, y=means[1], yend=means[1],size=2) +
    # #geom_segment(x=30-lw, xend=30+lw, y=means[2], yend=means[2],size=2) +
    # geom_text(x=-14,y=yticks[[i]][length(yticks[[i]])], size=(5/14)*18, family='sans', label=paste0("P=", round(sig,digits = 3)), hjust=0,vjust=1.5) +
    # theme_bw() +
    # #ylab(titles[i]) +
    # theme(legend.position="none",
          # axis.line = element_blank(), #element_line(color='gray', size=1),
          # panel.grid.minor = element_blank(),
          # panel.grid.major.y = element_blank(), #element_line(color='gray',size=1),
          # panel.grid.major.x = element_blank(),
          # panel.border = element_blank(),
          # plot.background = element_rect(color='black',fill='transparent'),
          # #element_rect(colour = "gray", fill = NA,size=1.5),
          # axis.title.x=element_blank(),
          # axis.ticks.x=element_line(color='gray',size=1),
          # axis.ticks.y=element_line(color='gray',size=1),
          # axis.ticks.length = unit(0.1, "in"),
          # #element_line(size=1.5,color='gray'),
          # axis.text.y = element_text(color="black",family="sans", size=18),
          # axis.text.x = element_text(color="black",family="sans", size=18,hjust=0.5),
          # axis.title = element_blank(), #element_text(colour="black",family="sans", size=24),
          # plot.title = element_text(color="black",family="sans", size=24,hjust=0.5)) +
    # scale_fill_manual(values=c('#FEBF00','#006FBF')) +
    # coord_cartesian(xlim=c(-14,46),ylim=ylims) +
    # scale_y_continuous(breaks=yticks[[i]], labels=yticks[[i]], expand = c(1e-8*(ylims[1]-ylims[2]), 1e-8*(ylims[2]-ylims[1]))) +
    # scale_x_continuous(labels=c('','Day 0','Day 30',''),breaks=c(-15,0,30,45), expand = c(0, 0)) +
    # ggtitle(titles[i])
    
    # for (j in 1:length(yticks[[i]]))
    # {
        # p = p + geom_segment(x=-15,xend=45, y = yticks[[i]][j], yend = yticks[[i]][j], color='gray',size=1)
    # }
    
    # #if (sig < 0.05)
    # #{
        # #p = p + geom_signif(comparisons = list(c("0", "30")), annotations="*", y_position = max(d.f$index)*1.1, tip_length = 0, size=2)
    # #}
    # #p = addConfidenceToStripChart(x=d.f$index, group=d.f$trmnt_day_s,p=p)
    
    # plots[[i]] = p;
    # assign(paste0("ind",i),d.f)
# }
return(p)
})
#svglite(file=paste0('alphaDiversity.svg'),bg='transparent',width=6.29*2,height=3.0*2,pointsize=24)
svg(filename=paste0('alphaDiversity.svg'),bg='transparent',width=7.09*2,height=6.09,pointsize=24, family='sans')
    multiplot(plotlist = plots,cols=length(vectors))
    
dev.off()
