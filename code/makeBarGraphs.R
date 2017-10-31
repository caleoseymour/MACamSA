library('dplyr')
library('ggplot2')
library('RColorBrewer')
library('knitr')
library('stringr')
library('svglite')

tex_dir = 'textures_2x/'

texture_list = list.files(tex_dir)
texture_list = texture_list[which(!vapply(texture_list, function(x) return(grepl('ZZ',x)),logical(1)))]
#texture_list = sample(texture_list)
textures = vapply(texture_list,function(x) return(image_uri(paste0(tex_dir,x))),character(1))
names(textures) = gsub(' ','_',names(textures))

makeBarPlotTextured = function(f, pal)
{
    pal = unique(pal)
    text_keep = textures
    names(text_keep) = gsub('.png','',names(text_keep))
    text_keep[pal == '#000000'] = text_keep["ZZBlack"]
    text_keep[pal == '#FFFFFF'] = text_keep["ZZWhite"]
    text_keep = text_keep[which(names(text_keep) != "ZZBlack")]
    text_keep = text_keep[which(names(text_keep) != "ZZWhite")]
    text_keep = text_keep[1:length(pal)]
    
    l = scan(f, what='character', sep='\n', quiet=TRUE)
    #print(l)
    #l = gsub('%','',l)
    #string_deletethis_editing <<- print(l[325])
    for (i in 1:length(l))
    {
        for(j in 1:length(text_keep))
        {
            if ((pal[j] != '#FFFFFF') & (pal[j] != '#000000'))
            {
                palstring = paste0('fill: ',pal[j],';')
                replacestring = paste0('fill:url(#texture_',names(text_keep)[j],'); fill-opacity:1.0')
                #print(palstring)
                #print(replacestring)
                l[i] = sub(palstring,replacestring,l[i])
            }
        }
        #print(hexToRGBPercent(pal[i]))
    }
    
    styleLine = which(l == "  ]]></style>")[1]
    buildTex = vector()
    for(j in 1:length(text_keep))
        {
            buildTex[j] = buildTextureSVG(text_keep[j])
        }
    
    #print(l)
    sink(paste0(f,'.modified.svg'))
    cat(paste(l[1:styleLine],collapse='\n'),'\n')
    cat(buildTex,collapse='\n\n','\n')
    cat(paste(l[(styleLine+1):length(l)],collapse='\n'))
    sink()
}

hexToRGBPercent = function(x)
{
    l = paste0('0x',unlist(strsplit(gsub("(.{2})", "\\1 ", sub('#','',x)),' ')))
    r = round((strtoi(l)*100)/255,6)
    o = paste('fill.rgb(',r[1],'%,',r[2],'%,',r[3],'%);',sep='')
    return(o)
}

buildTextureSVG = function(x)
{
# base = strsplit(x,',')
# base[[1]][2] = paste(base[[1]][2:length(base[[1]])],collapse=',')
# base[[1]] = base[[1]][1:2]
# base[[2]] = strsplit(gsub("(.{76})", "\\1 ", base[[1]][2]),' ')
# base[[2]] = paste(unlist(base[[2]]),collapse='\n')
# base[[1]] = base[[1]][1]
# base = unlist(base)
# base = paste(base,collapse=',')
    return(    paste('<pattern',
       'patternUnits="userSpaceOnUse"',
       'width="64"',
       'height="64"',
       paste0('id="texture_',names(x),'">'),
      '<image',
        'width="64"',
         'height="64"',
         'preserveAspectRatio="align"',
         paste0('xlink:href="',x,'"'),
         paste0('id="img_',names(x),'"'),
         'x="0"',
         'y="0" />',
    '</pattern>',
    sep='\n'))
}

plot_bar_ext = function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
    title = NULL, facet_grid = NULL, cutoff=0.05) 
{   
    rank_picked = which(rank_names(physeq) == fill)
    taxtable = as.matrix(attributes(physeq)[[2]])
    for(i in 2:ncol(taxtable))
    {
        unclassified = c(paste(paste0('D_',i-1,'__'),c('uncultured','uncultured bacterium',paste0('Unknown ',fill),''),sep=''),'unclassified','uncultured','uncultured bacterium')
        taxtable[(taxtable[,i] %in% unclassified), i] = paste('Unclassified',taxtable[(taxtable[,i] %in% unclassified), i-1],sep=' ')
        taxtable[,i] = apply(taxtable[,i],1,function(x) return (gsub('Unclassified [Uu]nclassified','Unclassified',x)))
        taxtable[,i] = apply(taxtable[,i],1,function(x) return (gsub('group','',x)))
    }

    attributes(physeq)[[2]] = taxtable
    physeq = tax_glom(physeq, taxrank=fill)
    physeq = transform_sample_counts(physeq, function(x) return(x/sum(x)))
    physeq = merge_taxa(physeq, taxa_names(physeq)[taxa_sums(physeq) < cutoff])
    #print(tax_table(physeq))
    
	tax_table(physeq)[is.na(tax_table(physeq)[,rank_picked])] = 'ZZ_LessThan'
    labs = gsub('MA','',sample_names(physeq))
    labs = gsub('0','0-',labs)
    sample_names(physeq) = labs
    print(labs)
    print(sample_names(physeq))
    mdf = psmelt(physeq)
    taxa = levels(with(mdf,get(fill)))
    print(taxa)
    taxa = gsub("D_[0-9]__","",taxa)
    #taxa = gsub("\\[","",taxa)
    #taxa = gsub("\\]","",taxa)
    eval(parse(text=paste0('levels(mdf$',fill,') = taxa')))
    pallette =  rep(c(
          rev(brewer.pal(9,"Set1")[2:8]),
          brewer.pal(7,"Set2"),
          brewer.pal(9,"Pastel1"),
          brewer.pal(8,"Dark2")
    ),times=30)
    
    ## Pallette modifications
    pallette[which(sort(unique(as.character(taxa))) == 'ZY_Unclassified')] = '#000000'
    pallette[which(sort(unique(as.character(taxa))) == 'ZZ_LessThan')] = '#FFFFFF'
    print(taxa)
    taxlabels = gsub('ZY_Unclassified', 'Unclassified',unique(sort(taxa)))
    taxlabels = gsub('ZZ_LessThan', paste0('Others (<', cutoff*100,'%)'),taxlabels)
    
    if (fill == 'Genus')
    {
        pallette = c('#c0017e','#a65628','#9ecae1','#081d58','#984ea3','#4daf4a','#ff0000','#878787','#fc8d62','#014636','#e78ac3','#a6d854','#ffd92f','#a50026','#fbb4ae','#ff7f00','#ccebc5','#decbe4','#80cdc1','#ffffcc','#ffffff')
    }

    mdf[,colnames(mdf) == fill] = as.character(mdf[,colnames(mdf) == fill])
    
    
    samlabs = as_labeller(c(`0` = 'Day 0', `30` = 'Day 30'))
    
    p = ggplot(arrange(mdf,get(fill)), aes_string(x = 'Sample', y = 'Abundance', fill = fill)) +
    geom_bar(stat = "identity", position = "stack", color = 'black', width=0.8, size=1) +
    facet_wrap( ~ trmnt_day_s, scales = "free_x", strip.position="bottom", labeller = samlabs) +
    theme_bw() +
    theme(#text = element_text(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, colour="black",size=24),
          axis.text.y = element_text(hjust = 1, vjust=0.5, colour="black",size=24),
          axis.title.x=element_blank(),
          axis.title.y=element_text(vjust=1,colour="black",size=24),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.x = element_text(colour="black",size=24),
          strip.placement = "outside",
          legend.text=element_text(color='black', size=24),
          legend.title=element_text(color='black', size=24),
          strip.background = element_blank(),
          panel.border = element_blank(),
          #rect(color='black',size=1),
          plot.background = element_blank(),
          panel.background = element_blank(),
          legend.key.size =  unit(1.75, "lines"),
          legend.position=c(1.0,1.01),
          legend.key = element_blank(),
          legend.justification = c("left","top"),
          plot.margin = unit(c(0.1,12,0.02,0.62), "null"),
          legend.background = element_blank(),
          text=element_text(family="sans"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()
          #legend.position = 'bottom'
          ) +
          guides(fill = guide_legend(ncol = 1)) +
    #scale_x_discrete(labels=labs) +
    scale_fill_manual(values=pallette,labels=taxlabels) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0,0.25,0.5,0.75,1),labels=c('0%','25%','50%','75%','100%')) +
    labs(y = 'Relative Abundance (Percent)')

    gt = ggplotGrob(p)
    N <- mdf %>% group_by(trmnt_day_s) %>% 
         summarise(count = length(unique(Sample))) %>% 
         `[[`(2)

    # Get the column index in the gt layout corresponding to the panels.
    panelI <- gt$layout$l[grepl("panel", gt$layout$name)]

    # Replace the default panel widths with relative heights.
    gt$widths[panelI] <- unit(N, "null")

    # Add extra width between panels (assuming two panels)
    gt$widths[panelI[1] + 1] = unit(0.2, "cm")

    ## Draw gt
    
    svglite(file=paste0('barplots_',fill),bg='transparent',width=(6.89+0.89+0.5)*(9.49/7.12)+0.2,height=(6.89+0.89+0.5)*(9.49/7.12),pointsize=24)
    grid.newpage()
    grid.draw(gt)
    dev.off()
    dev.off()
    makeBarPlotTextured(paste0('barplots_',fill),pallette)
    file.rename(paste0('barplots_',fill), paste0('barplots_',fill,'.color.svg'))
    

    # Replace the default widths with relative widths:
    #treatments = as.character(sample_data(physeq)$trmnt_type_s)
    #print(treatments)
    #nulls = which(as.logical(vapply(z$widths, function(x) return(as.character(x) == "1null"),1)))
    #print(nulls)
    #print(gt)
    #for (i in 1:length(treatments))
    #{
    #    gt$widths[nulls[i]] = unit(sum(treatments == treatments[i]), "null")
    #}
    #print(gt)
    

    # Draw the plot
    #grid.newpage()
    #grid.draw(gt)
}
