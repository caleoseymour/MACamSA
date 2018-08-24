library('phyloseq')

## importPhyseqFromPlaintext: Import and join physeq objects from filenames.
importPhyseqFromQiime2 = function(otutable=NULL,samdata=NULL,taxtable=NULL)
{
    ## Parse files using specialized functions
    ## Prioritize for OTU table completeness.
    otumat = NULL
    samdf = NULL
    taxmat = NULL
    if (!is.null(otutable))
    {
        otumat = readTSVasOTUs(otutable)
        OTU = otu_table(otumat, taxa_are_rows = TRUE)
    }
    if (!is.null(taxtable))
    {
        taxmat = readTSVasTAX(taxtable)
        TAX = tax_table(taxmat)
    }
    if (!is.null(samdata))
    {
        samdf = readTSVasSAMPLEDATA(samdata)
        SAM = sample_data(samdf)
    }
    
    ## Create phyloseq object
    physeq = merge_phyloseq(OTU,TAX,SAM)
    #taxa_names(physeq) = paste0('feature',c(1:ntaxa(physeq)))
    physeq = prune_taxa(taxa_sums(physeq) > 0, physeq)
    #attributes(physeq) = list(otu_table = attributes(physeq)[[1]],
    #                          tax_table = attributes(physeq)[[2]],
    #                          sam_data = attributes(physeq)[[3]],
    #                          phy_tree = attributes(physeq)[[4]],
    #                          ref_seq = attributes(physeq)[[5]],
    #                          class = attributes(physeq)[[6]],
    #                          feature_names = taxa_names(OTU))
    return(physeq)
}

## Read in a Qiime2 TSV feature table as a matrix
readTSVasOTUs = function(otutable)
{
    OTU = as.matrix(read.table(otutable,sep='\t',row.names=1,header=FALSE))
    txt = scan(otutable,what='character',sep='\n',n=2,quiet=TRUE)[2]
    splitxt = unlist(strsplit(txt,'\t'))
    samnames = splitxt[2:length(splitxt)]
    colnames(OTU) = samnames
    
    return(OTU)
}

## Read in Qiime2 compatible sample metadata as a dataframe
readTSVasSAMPLEDATA = function(samdata)
{
    SAM = read.table(samdata,sep='\t',row.names=1,header=FALSE)
    txt = scan(samdata,what='character',sep='\n',n=2,quiet=TRUE)[1]
    splitxt = unlist(strsplit(txt,'\t'))
    samnames = splitxt[2:length(splitxt)]
    colnames(SAM) = samnames
    
    return(SAM)
}

## Read in a Qiime2 taxonomy table as a 7-column taxonomy matrix
readTSVasTAX = function(taxtable,nlevel=7)
{
    tt = read.table(taxtable,sep=c('\t'),row.names=1,header=FALSE)
    CONF = tt[ncol(tt)]
    tt = tt[-ncol(tt)]
    taxa_list = strsplit(as.character(tt[,1]),';')
    taxa_len = max(vapply(taxa_list, function(x) return(length(x)), FUN.VALUE=integer(1)))
    taxa_fixed = lapply(taxa_list, function(x) return(c(x,unlist(replicate(taxa_len-length(x),"unclassified")))))
    taxa_7level = lapply(taxa_fixed, function(x) return(x[1:nlevel]))
    
    TAX = matrix(unlist(taxa_7level), nrow = length(taxa_list), ncol = 7, byrow=TRUE)
    rownames(TAX) = rownames(tt)
    colnames(TAX) = c('Domain','Phylum','Class','Order','Family','Genus','Species')
    
    return(TAX)
}