require(oligo)
require(pd.mogene.2.1.st)
require(Heatplus)
require(reshape2)
require(ReportingTools)
require(hwriter)
require(ggplot2)
require(limma)
require(clusterProfiler)
require(XML)
require(R.utils)

#Given a vector of sample names in the form: X_Y produce colors based off of the value of X
make.colors.boxplot <- function(sample.names)
{
    col.strains <- sapply(strsplit(sample.names, "_"), "[", 1)
    use.cols <- darkColors(length(unique(col.strains)))
    names(use.cols) <- unique(col.strains)
    return(as.character(use.cols[col.strains]))
}

#A convienience function for making expression boxplots 
#Input:
#       raw.exprs: A GeneFeatureSet or similar object which has an rma method returning an ExpressionSet
#                   Note that if the GeneFeatureSet has pData with a 'RIN' or 'Weight' column, these will be added to the plot
#       base.name: The base name of the PDF to generate (if make.pdf == T)
#       type:  One of raw or norm indicating whether to rma normalize ("norm") the expression prior to generating the boxplot or not ("raw")
#       make.pdf: A logical value indicating whether to generate a PDF or simply plot it to the screen
#       highlight.names: A character vector containing the sample names to be colored differently
#Output:
#       A PDF named as base.name_type.pdf
make.boxplot <- function(raw.exprs, base.name="test", type=c("raw", "norm"), make.pdf=T, highlight.names=NULL)
{
    type <- match.arg(type)
    
    use.exprs <- switch(type, raw=rma(raw.exprs, normalize=FALSE, target="core"), norm=rma(raw.exprs, normalize=TRUE, target="core"))

    basic.ord <- do.call("order", pData(use.exprs)[,c("Lab", "Tissue", "Mating", "Inf_Status", "Timepoint")])
    
    use.exprs <- use.exprs[,basic.ord]
    
    if (make.pdf)
    {
        pdf(file=paste0(reports.dir,"/QA_QC/", base.name, "_", type, ".pdf"), width=16, height=10) # updated to save PDF in directory (GC)
    }
    
    layout(c(1,1,1,1,2,3,3,3))
    par(mar=c(0,4,4,2))
    
    # Update title (GC)
    boxplot(exprs(use.exprs), xaxt="n", ylab="Expression", main=paste(capitalize(type),"Expression", sep=" "), outline=FALSE, col=make.colors.boxplot(colnames(use.exprs)))
    
    #RIN subplot
    if ('RIN' %in% names(pData(use.exprs)))
    {
        par(mar=c(0,4,0,2))
        plot(x=1:ncol(use.exprs), y=pData(use.exprs)$RIN, xlim=c(1, ncol(use.exprs)) + c(-.5, .5), ylim=range(pData(use.exprs)$RIN, na.rm=T), xaxt="n", type="p", ylab="RIN", xlab="", col=make.colors.boxplot(colnames(use.exprs)))
        abline(h=7, lty="dashed")
    }
    
    #Weight subplot
    if ('Weight' %in% names(pData(use.exprs)))
    {
        par(mar=c(15,4,0,2))
        plot(x=1:ncol(use.exprs), y=pData(use.exprs)$Weight, xlim=c(1, ncol(use.exprs)) + c(-.5, .5), ylim=range(pData(use.exprs)$Weight, na.rm=TRUE), xaxt="n", type="b", ylab="Weight", xlab="", col=make.colors.boxplot(colnames(use.exprs)))
    }
    
    #lables and such
    if (missing(highlight.names) || is.null(highlight.names))
    {
        Axis(at=1:ncol(use.exprs), side=1, labels=colnames(use.exprs), las=2)
    }
    else if(is.character(highlight.names) && all(highlight.names %in% colnames(use.exprs)))
    {
        should.col <- ifelse(colnames(use.exprs) %in% highlight.names, 'red', 'black')
        rle.col <- Rle(should.col)
        
        lo.pos <- 1
        
        for (i in 1:nrun(rle.col))
        {
            positions <- lo.pos:(runLength(rle.col)[i] + lo.pos-1)
            Axis(at=positions, side=1, labels=colnames(use.exprs)[positions], las=2, col.axis=runValue(rle.col)[i])
            lo.pos <- max(positions) + 1
        }
        
        #Axis(at=1:ncol(use.exprs), side=1, labels=colnames(use.exprs), las=2, col.axis=ifelse(should.col, 'red', 'black'))
        
        Axis(at=26, side=1, labels=colnames(use.exprs)[26], las=2, col.axis='red')
        
    }
    
    if (make.pdf)
    {
        dev.off()
    }
    
    
}

#Given a vector of time values in the form: 'prefix'timepoint (e.g. D10) produce an ordering such that the
#earliest timepoint comes first. Useful for producing ordered factors.
order.times <- function(time.vals, prefix="D")
{
    unique.times <- unique(as.character(time.vals))
    time.ord <- order(as.numeric(sub(prefix, "", unique.times)), decreasing=F)
    time.levs <- unique.times[time.ord]
    return(time.levs)
}

#Given a FeatureSet or ExpressionSet containing pData with the following columns: Mating,Lab,Tissue,Inf_Status and Timepoint
#produce a data.frame summary of the corresponding counts.
make.summary.table<- function(raw.exprs)
{
    tab.dta <- pData(raw.exprs)
    
    unique.times <- unique(as.character(tab.dta$Timepoint))
    time.ord <- order(as.numeric(sub("D", "", unique.times)), decreasing=F)
    time.levs <- unique.times[time.ord]
    
    tab.dta$Timepoint <- factor(tab.dta$Timepoint, levels=time.levs, labels=time.levs)
    full.table <- dcast(formula=Mating+Lab+Tissue+Inf_Status~Timepoint, data=tab.dta, fun.aggregate=length, value.var="Timepoint")
    return(full.table)
}


#A function that provides a consistent way of processing expression data for heatmaps etc.
#Input:
#   raw.exprs needs to be FeatureSet
#   num.genes should either be an integer value indicating the 'num.genes' most variable genes to return
#           or should not be specified at all.       
#Output:
#   An ExpressionSet
heatmap.process <- function(raw.exprs, num.genes)
{
    norm.exprs.basic <- rma(raw.exprs, target="core")
    norm.exprs.basic.mat <- exprs(norm.exprs.basic)
    if(missing(num.genes) || is.null(num.genes) || is.na(num.genes))
    {
        return(norm.exprs.basic.mat)
    }else{
         exprs.var <- apply(norm.exprs.basic.mat, 1, var)
        norm.exprs.basic.var <- norm.exprs.basic.mat[order(exprs.var, decreasing=TRUE),][1:num.genes,]
        return(norm.exprs.basic.var)
    }
   
}

#A function that provides a simple way to retrieve the group assignments for samples base on standard hclust
#Input:
#   raw.exprs: needs to be FeatureSet
#   num.groups: is the number of groups to form by cutting the tree
#   num.genes: is the number of genes to pass to heatmap.process 
#Output:
#   A list containing a data.frame of group assignmnets and the dendrogram from hclust
get.tree.groups <- function(raw.exprs, num.groups=2, num.genes=1000)
{
    norm.exprs.basic.var <- heatmap.process(raw.exprs, num.genes)
    
    cur.clust <- hclust(dist(t(norm.exprs.basic.var)))

    samp.assign <- as.data.frame(cutree(cur.clust, k=num.groups))
    samp.assign$sample <- rownames(samp.assign)
    names(samp.assign) <- c("group", "sample")
    return(list(groups=samp.assign[,c("sample", "group")], dendro=cur.clust))
}

#A convienience function for making expression heatmaps
#Input:
#       raw.exprs: A GeneFeatureSet or similar object which has an rma method returning an ExpressionSet
#       base.name: The base name of the PDF to generate (if make.pdf == T)
#       cut.dist: The cut height to use for coloring groups in annHeatmap2 (get.tree.groups is useful here)
#       num.genes: is the number of genes to pass to heatmap.process 
#       make.pdf: A logical value indicating whether to generate a PDF or simply plot it to the screen
#       rin.bin.func: A function to use to form RIN score quality bins, probably should use 'cut' internally as the default does.
#       weight.bin.func: Unimplemented
#       include.tissue: Should tissue be included in the annotations.
#Output:
#       A PDF named as base.name_heatmap.pdf

make.heatmap <- function(raw.exprs, base.name="test", cut.dist=NULL, num.genes=1000, make.pdf=T, rin.bin.func=function(x) cut(x, breaks=c(0,5, 7, 10)), weight.bin.func=NULL, include.tissue=F)
{
    
   norm.exprs.basic.var <- heatmap.process(raw.exprs, num.genes)
    
    all.dta <- pData(raw.exprs)
    
    stopifnot(all(rownames(all.dta) == colnames(norm.exprs.basic.var)))
    
    if (missing(cut.dist) || is.null(cut.dist) || is.na(cut.dist))
    {
        cluster <- NULL
    }else{
        cluster <- list(cuth=cut.dist)
    }
    
    if (missing(weight.bin.func) == F && is.null(weight.bin.func) == F)
    {
        stop("ERROR: Haven't implemented adding weight binning function yet...")
    }
    
    base.factors <- c("Mating", "Inf_Status", "Timepoint")
    
    if (include.tissue)
    {
        base.factors <- append(base.factors, 'Tissue')
    }
    
    if (is.function(rin.bin.func) && 'RIN' %in% names(all.dta))
    {
        all.dta$RIN <- rin.bin.func(all.dta$RIN)
        use.dta <- all.dta[,append(base.factors, "RIN")]
    }else{
        use.dta <- all.dta[,base.factors]
    }
    
    new.heat <- annHeatmap(norm.exprs.basic.var, annotation=use.dta, dendrogram = list(clustfun = hclust, distfun = dist, Col = list(status = "yes"), Row = list(status = "hidden")),
                           labels=list(Row=list(labels=rep("", nrow(norm.exprs.basic.var))), Col=list(nrow=0, labels=rep("", ncol(norm.exprs.basic.var)))), legend = TRUE, cluster = cluster)
    
    #cuth was chosen by examination of plot(new.heat$dendrogram$Col$dendro)
    
    if (make.pdf)
    {
        pdf(paste0(reports.dir,"/QA_QC/",base.name, "_heatmap.pdf"), width=16, height=10)
        plot(new.heat)
        dev.off()
    }
    else
    {
        plot(new.heat)
    }
}

#A convienience function for making bacterial spike plots
#Input:
#       raw.exprs: A GeneFeatureSet or similar object which has an rma method returning an ExpressionSet
#       base.name: The base name of the PDF to generate (if make.pdf == T)
#       pgf.file: probe group file from Affy for the appropriate array type.
#Output:
#       A PDF named as base.name_bac_spikes.pdf

plot.bac.spikes <- function(raw.exprs, base.name="test", make.pdf=T, pgf.file="./MoGene-2_1-st.pgf")
{
    require(ggplot2)
    require(affxparser)
    use.exprs <- rma(raw.exprs)
    
    pgf.list <- readPgf(pgf.file)
    probeset.types <- data.frame(fsetid=pgf.list$probesetId, fsetName=pgf.list$probesetName, type=pgf.list$probesetType, stringsAsFactors=FALSE)
    
    spike.probes <- probeset.types[probeset.types$type == 'control->affx->bac_spike',]

    spike.exprs <- exprs(use.exprs)[as.character(unique(spike.probes$fsetid)),]
    
    spike.dta <- melt(spike.exprs)
    spike.dta.merge <- merge(spike.dta, spike.probes, by.x="Var1", by.y="fsetid", all=F, incomparables=NA, sort=FALSE)
    
    spike.dta.merge$`Bac. Spike` <- sapply(strsplit(spike.dta.merge$fsetName, "-"), "[[", 4)
    #should be BioB<BioC<BioD<Cre
    spike.dta.merge$`Bac. Spike` <- factor(spike.dta.merge$`Bac. Spike`, levels=c("cre", "bioD", "bioC", "bioB"), labels=c("Cre", "BioD", "BioC", "BioB"), ordered=TRUE)
    
    bac.plot <- qplot(x=Var2, y=value, group=`Bac. Spike`, color=`Bac. Spike`, data=spike.dta.merge, stat="summary",
            fun.y=mean, geom="line", ylab="log2(Expression)", xlab="", main="Bacterial Spikes")
        bac.plot <- bac.plot + theme(axis.text.x=element_text(size=8, angle=90, hjust=1))
    
    if (make.pdf)
    {
        ggsave(bac.plot, file=paste0(reports.dir,"/QA_QC/",base.name, "_bac_spikes.pdf"), width=16, height=10)
    }else{
        plot(make.plot)
    }
}


#A convienience function for making polya spike plots
#Input:
#       raw.exprs: A GeneFeatureSet or similar object which has an rma method returning an ExpressionSet
#       base.name: The base name of the PDF to generate (if make.pdf == T)
#       plot.type: Either "AFFX-r2-Bs" or "AFFX" which are the two types of probesets, defaults to "AFFX-r2-Bs"
#       pgf.file: probe group file from Affy for the appropriate array type.
#Output:
#       A PDF named as base.name_polya_spikes.pdf
plot.polya.spikes <- function(raw.exprs, base.name="test", make.pdf=T, plot.type=c("AFFX-r2-Bs", "AFFX"), pgf.file="./MoGene-2_1-st.pgf")
{
    require(ggplot2)
    require(affxparser)
    
    plot.type <- match.arg(plot.type)
    
    use.exprs <- rma(raw.exprs)
    
    pgf.list <- readPgf(pgf.file)
    probeset.types <- data.frame(fsetid=pgf.list$probesetId, fsetName=pgf.list$probesetName, type=pgf.list$probesetType, stringsAsFactors=FALSE)
    
    poly.a.spikes <- probeset.types[probeset.types$type == "control->affx->polya_spike",]
    
    spike.exprs <- exprs(use.exprs)[as.character(unique(poly.a.spikes$fsetid)),]
    
    spike.dta <- melt(spike.exprs)
    spike.dta.merge <- merge(spike.dta, poly.a.spikes, by.x="Var1", by.y="fsetid", all=F, incomparables=NA, sort=FALSE)
    
    spike.dta.merge$type <- ifelse(grepl("AFFX-r2-Bs", spike.dta.merge$fsetName), "AFFX-r2-Bs", "AFFX")
    spike.dta.merge$pos <- sub("_*s*_[as]t", "", sapply(strsplit(spike.dta.merge$fsetName, "-"), function(x) x[length(x)]))
    spike.dta.merge$pos <- factor(spike.dta.merge$pos, levels=c("5", "M", "3"), ordered=T)
    
    spike.dta.merge$direction <- sapply(strsplit(spike.dta.merge$fsetName, "[-_]"), function(x) x[length(x)])
    
    
    split.fset <- strsplit(spike.dta.merge$fsetName, "-")
    spike.dta.merge$Spike <- ifelse(spike.dta.merge$type == "AFFX", sapply(split.fset, "[", 2), sapply(split.fset, "[", 4))
    spike.dta.merge$Spike <- sub("X", "", spike.dta.merge$Spike)
    substr(spike.dta.merge$Spike, 1, 1) <- toupper(substr(spike.dta.merge$Spike, 1, 1))
    
    #from the affy exon/gene array whitepaper: lys<phe<thr<dap, not sure about Trpn
    spike.dta.merge$Spike <- factor(spike.dta.merge$Spike, levels=c("Dap", "Thr", "Phe", "Lys", "Trpn"), ordered=T)
    
    pa.basic.plot <- qplot(x=Var2, y=value, group=Spike, color=Spike, data=spike.dta.merge[spike.dta.merge$type == plot.type,], stat="summary",
            fun.y=mean, geom="line", ylab="log2(Expression)", xlab="", main="PolyA Spikes", facets=pos~direction) + theme(axis.text.x=element_text(size=8, angle=90, hjust=1))
    
    if(make.pdf){
         ggsave(pa.basic.plot, file=paste0(reports.dir,"/QA_QC/",base.name, "_polya_spikes.pdf"), width=16, height=10)
    }else{
        plot(pa.basic.plot)
    }
    
    
}

#Generates a more attractively formatted HTML table for simple numeric summaries than ReportingTools or hwriter
#Input:
#       summary.df:  A data.frame which contains several columns of character or factor followed by numeric summaries such as counts such as that
#                   produced by 'make.summary.table'
#       cast.form: An option formula as needed by reshape2::dcast which if applied prior to conversion to HTML.
#       total.col: If cast.form is supplied and non-NULL, then total.col is used as the 'value.var' for reshape2::dcast
#       sep.header: Should headers of the form 'table.header' be interpreted as seperate portions of the table.
#       rec.limit: The number of columns to process in terms of <th> tags.  Any left over columns will be converted to a normal (<tr><td>) HTML table
#Output:
#       A string containing the HTML representation of summary.df after application of cast.form.

#assumes bootstrap CSS, which should be the case for reportingtools
make.pretty.df <- function(summary.df, cast.form=NULL, total.col=NULL, sep.header=F, rec.limit=10000)
{
    ret_str <- '<table class="table table-hover table-bordered table-striped", style="width:auto">'
    
    if (sep.header)
    {
        #seperate the header by names of the form: table.header
        
        split.heads <- strsplit(colnames(summary.df), "_")
        
        which.to.sep <- elementLengths(split.heads) > 1
        
        table.val <- character(ncol(summary.df))
        table.val[which.to.sep] <- sapply(split.heads[which.to.sep], "[", 1)
        
        col.val <- colnames(summary.df)
        col.val[which.to.sep] <- sapply(split.heads[which.to.sep], "[", 2)
        
        table.rle <- Rle(table.val)
        
        ret_str <- append(ret_str, paste('<thead><tr>', paste(paste0('<th colspan="',runLength(table.rle),'" style="text-align:center">', runValue(table.rle), '</th>'), collapse=""), '</tr>'))
        
        ret_str <- append(ret_str, paste('<tr>', paste(paste0('<th>', col.val, '</th>'), collapse=""), '</tr></thead>'))
    }else{
        ret_str <- append(ret_str, paste('<thead><tr>', paste(paste0('<th>', colnames(summary.df), '</th>'), collapse=""), '</tr></thead>'))
    }
    
    ret_str <- append(ret_str, '<tbody>')
    
    #need to also check total.col
    if (missing(cast.form) || is.null(cast.form) || is.na(cast.form))
    {
        use.df <- summary.df
    }else{
        use.df <- dcast(formula=as.formula(cast.form), data=summary.df, value.var=total.col)
    }
    
    col.ords <- sapply(1:ncol(use.df),function(x) is.character(use.df[,x]) || is.factor(use.df[,x]))
    
    #if all the character cols don't come before the numeric, reorder so it is so...
    if (identical(runValue(Rle(col.ords)), c(T, F)) == F)
    {
        use.df <- use.df[,order(col.ords, decreasing=T)]
    }
    
    html.list <- lapply(1:ncol(use.df), function(i)
    {
        if ((is.character(use.df[,i]) || is.factor(use.df[,i])) && (rec.limit > i))
        {
            df.runs <- Rle(use.df[,i])
            
            temp <- mapply(function(len, val){
                return(append(paste0('<th rowspan="',len,'" valign="top">', val,'</th>'), rep("", len-1)))
            }, runLength(df.runs), runValue(df.runs), SIMPLIFY=F, USE.NAMES=F)
            
            return(unlist(temp))
            
        }else{
            return(paste('<td>',use.df[,i],'</td>'))
        }
    })
    
    #also should order by size of the rowspan
    
    for (i in 1:nrow(use.df))
    {
        ret_str <- append(ret_str, '<tr>')
        
        for(j in 1:ncol(use.df))
        {
            ret_str <- append(ret_str, html.list[[j]][i])
        }
        
        ret_str <- append(ret_str, '</tr>')
    }
    
    ret_str <- append(ret_str, '</tbody>')
    ret_str <- append(ret_str, '</table>')
    
    return(paste(ret_str, collapse=""))
}


#Returns a clusterProfResult object resulting from significance results in sig.list
#Input:
#       sig.list:  A named list which contains an object (matrix or data.frame) with mogene21st probeIDs as rownames.  This should be subsetted for only the significant results.
#           Such as list is generated by contrast.summary.tc
#       univ: A data.frame as produced by select which contains the universe of Entrez ID as the univ$ENTREZID column.
#       make.readable: A logical value indicating whether compareCluster should output gene symbols (T) or entrez IDs (F)

#Output:
#       A clusterProfResult object 

basic.go.enrich.tc <- function(sig.list, univ, make.readable=F)
{
    use.list <- lapply(sig.list, function(x)
                       {
                            if (nrow(x) > 0)
                            {
                                annot.temp <- select(mogene21sttranscriptcluster.db, keys=rownames(x), columns=c("ENTREZID"), keytype="PROBEID")
                            
                                stopifnot(nrow(annot.temp) == nrow(x))
                                
                                return(annot.temp[,"ENTREZID"])
                            }else{
                                return(NULL)   
                            }
                       })
    
    use.list <- use.list[sapply(use.list, is.null)==F]
    
    ck <- compareCluster(use.list, fun="enrichGO", organism="mouse", ont="BP", pvalueCutoff = 0.05, pAdjustMethod = "BY", universe=univ$ENTREZID, qvalueCutoff = 1, minGSSize = 5, readable=make.readable)

    return(ck)
}


#Computes differential expression relative to mock samples seperately for each cross using limma.  Warnings will be emitted recording which comparisions cannot be carried out.
#Input:
#       use.exprs: ExpressionSet with a pData slot containing a data.frame with Mating, Timepoint, Inf_Status and optionally Tissue
#       contrast.correction: How to correct for multiple testing should be a valid value to provide to the limma::decideTests function
#       use.mock.time: If specified in the form: "DX", the mock timepoint to perform the comparisions relative to.  If not specified, it will choose the closest available mock.
#       subset.contrasts.by.tissue: A logical value indicating whether the ExpressionSet should be subsetted by tissue first before carrying out the comparisons.  Tissue needs to be a column in the pData slot of use.exprs
#Output:
#       A list containing the following elements:
#           -results: A list containing the following elements named by Mating:
#               1) coefs: The contrast coefficients
#               2) infecteds: The matrix of 0,1 and -1's formed from running limma::decideTests with method=contrast.correction
#           -fit: The MArrayLM object from limma.
#           -cont.counts: The number of samples for the two parts of each formed contrast or NA if the samples were not available

de.bycross.to.mock <- function(use.exprs, contrast.correction="separate", use.mock.time=NULL, subset.contrasts.by.tissue=NULL)
{
    if (missing(subset.contrasts.by.tissue) == F && is.null(subset.contrasts.by.tissue) == F)
    {
        use.exprs <- use.exprs[,pData(use.exprs)$Tissue == subset.contrasts.by.tissue]
    }
    
    colnames(use.exprs) <- paste0("M", colnames(use.exprs))
    
    exprs.dta <- pData(use.exprs)
    exprs.dta$Mating <- paste0("M", exprs.dta$Mating)
    
    time.exp <- with(exprs.dta, paste(Mating, Timepoint, Inf_Status, sep="."))
    
    non.mock <- setdiff(exprs.dta$Inf_Status, "M")
    stopifnot(length(non.mock)==1)
    
    if (missing(use.mock.time) || is.null(use.mock.time) || is.na(use.mock.time))
    {
        mock.days <- unique(use.exprs$Timepoint[use.exprs$Inf_Status == "M"])
    }else{
        mock.days <- use.mock.time
    }
    
    mock.days.int <- as.integer(sub("D", "", mock.days))
    
    mouse.time <- expand.grid(list(unique(exprs.dta$Mating), unique(exprs.dta$Timepoint)))
    
    inf.days <- as.integer(sub("D", "", as.character(mouse.time$Var2)))
    
    closest.mock <- Biobase::matchpt(inf.days, mock.days.int)
    
    mouse.time <- cbind(mouse.time, non.mock, mouse.time$Var1, mock.days[closest.mock$index], "M")
    
    mouse.conts <- paste(paste(mouse.time[,1], mouse.time[,2], mouse.time[,3], sep="."), paste(mouse.time[,4], mouse.time[,5], mouse.time[,6], sep="."), sep="-")
    
    paste.mouse.inf <- paste(mouse.time[,1], mouse.time[,2], mouse.time[,3], sep=".")
    paste.mouse.mock <- paste(mouse.time[,4], mouse.time[,5], mouse.time[,6], sep=".")
    
    cont.counts <- mapply(function(x,y) return(c(x,y)), as.numeric(table(time.exp)[paste.mouse.inf]), as.numeric(table(time.exp)[paste.mouse.mock]), SIMPLIFY=F)
    names(cont.counts) <- mouse.conts
    
    #also should count the number of samples involved for each contrast and potentially flag those with < 2 samples...
        #the flag would be part of res.list or the returned struture
    
    diff.mouse.inf <- setdiff(paste.mouse.inf, time.exp)
    diff.mouse.mock <- setdiff(paste.mouse.mock, time.exp)
    
    if(length(diff.mouse.inf) > 0 || length(diff.mouse.mock) > 0)
    {
        rm.contrasts <- which((paste.mouse.inf %in% diff.mouse.inf) | (paste.mouse.mock %in% diff.mouse.mock))
        warning(paste("Removing incomplete contrasts:", paste(mouse.conts[rm.contrasts], collapse=" , ")))
        
        mouse.conts <- mouse.conts[-rm.contrasts]
    }
    
    mod <- model.matrix(~0+time.exp)
    colnames(mod) <- sub("time\\.exp", "", colnames(mod)) 
    
    fit.1 <- lmFit(use.exprs, mod)
    
    conts <- makeContrasts(contrasts=mouse.conts, levels=mod)
    
    fit.2 <- contrasts.fit(fit.1, conts)
 
    fit.2 <- eBayes(fit.2)
    
    unique.matings <- unique(exprs.dta$Mating)
    
    res.list <- lapply(unique.matings, function(x)
           {
                print(x)
                where.line <- grep(x, colnames(fit.2$coefficients))
                if (length(where.line) > 0)
                {
                    return(list(infecteds=decideTests(fit.2[,where.line], method=contrast.correction), coefs=fit.2$coefficients[,where.line,drop=F]))
                }
                else
                {
                    return(NULL)
                }
                
           })
    
    names(res.list) <- unique.matings
    
    null.list <- sapply(res.list, is.null)
    
    return(list(results=res.list[null.list == F], fit=fit.2, cont.counts=cont.counts))
}

#Computes differential expression relative to other timepoints for each cross using limma.  Warnings will be emitted recording which comparisions cannot be carried out.
#Input:
#       use.exprs: ExpressionSet with a pData slot containing a data.frame with Mating, Timepoint, Inf_Status.  Tissue specific comparisons are not currently supported.
#       contrast.correction: How to correct for multiple testing should be a valid value to provide to the limma::decideTests function
#       run.mocks: A logical value indicating whether mock-only timepoint comparisons should be carried out as well.  If it is, then a 'mocks' element is added to the 'results' list.
#Output:
#       A list containing the following elements:
#           -results: A list containing the following elements named by Mating:
#               1) coefs: The contrast coefficients
#               2) infecteds: The matrix of 0,1 and -1's formed from running limma::decideTests with method=contrast.correction
#               3) mocks: Same as infecteds but for the mock comparisons if run.mocks = TRUE
#           -fit: The MArrayLM object from limma.
#           -cont.counts: The number of samples for the two parts of each formed contrast or NA if the samples were not available

de.bycross.to.timecourse <- function(use.exprs, contrast.correction="separate", run.mocks=F)
{
    colnames(use.exprs) <- paste0("M", colnames(use.exprs))
    
    exprs.dta <- pData(use.exprs)
    exprs.dta$Mating <- paste0("M", exprs.dta$Mating)
    
    time.exp <- with(exprs.dta, paste(Mating, Timepoint, Inf_Status, sep="."))
    
    time.combs <- t(combn(unique(exprs.dta$Timepoint), 2))
    time.combs <- t(apply(time.combs, 1, sort, decreasing=T))
    
    mouse.time <- expand.grid(list(unique(exprs.dta$Mating), 1:nrow(time.combs)))
    
    non.mock <- setdiff(exprs.dta$Inf_Status, "M")
    stopifnot(length(non.mock)==1)
    
    paste.mouse.1 <- paste(mouse.time[,1], time.combs[mouse.time[,2],1], non.mock, sep=".")
    paste.mouse.2 <- paste(mouse.time[,1], time.combs[mouse.time[,2],2], non.mock, sep=".")
    
    if (run.mocks)
    {
        paste.mouse.1 <- append(paste.mouse.1, paste(mouse.time[,1], time.combs[mouse.time[,2],1], "M", sep="."))
        paste.mouse.2 <- append(paste.mouse.2, paste(mouse.time[,1], time.combs[mouse.time[,2],2], "M", sep="."))
    }
    
    mouse.conts <- paste(paste.mouse.1, paste.mouse.2, sep="-")
    
    cont.counts <- mapply(function(x,y) return(c(x,y)), as.numeric(table(time.exp)[paste.mouse.1]), as.numeric(table(time.exp)[paste.mouse.2]), SIMPLIFY=F)
    names(cont.counts) <- mouse.conts
    
    diff.mouse.1 <- setdiff(paste.mouse.1, time.exp)
    diff.mouse.2 <- setdiff(paste.mouse.2, time.exp)
    
    if(length(diff.mouse.1) > 0 || length(diff.mouse.2) > 0)
    {
        rm.contrasts <- which((paste.mouse.1 %in% diff.mouse.1) | (paste.mouse.2 %in% diff.mouse.2))
        warning(paste("Removing incomplete contrasts:", paste(mouse.conts[rm.contrasts], collapse=" , ")))
        
        mouse.conts <- mouse.conts[-rm.contrasts]
    }
    
    mod <- model.matrix(~0+time.exp)
    colnames(mod) <- sub("time\\.exp", "", colnames(mod)) 
    
    fit.1 <- lmFit(use.exprs, mod)
    
    conts <- makeContrasts(contrasts=mouse.conts, levels=mod)
    
    fit.2 <- contrasts.fit(fit.1, conts)
 
    fit.2 <- eBayes(fit.2)
    
    unique.matings <- unique(exprs.dta$Mating)
    
    res.list <- lapply(unique.matings, function(x)
           {
                where.line <- grepl(x, colnames(fit.2$coefficients))
                if (sum(where.line) > 0)
                {
                    if (run.mocks)
                    {
                        where.mocks <- grepl("D\\d+\\.M", colnames(fit.2$coefficients))
                        
                        line.mocks <- where.line & where.mocks
                        line.inf <- where.line & (where.mocks == F)
                        
                        return(list(infecteds=decideTests(fit.2[,line.inf], method=contrast.correction),  coefs=fit.2$coefficients[,line.inf], mocks=fit.2[,line.mocks]))
                        
                    }else{
                        return(list(infecteds=decideTests(fit.2[,where.line], method=contrast.correction),  coefs=fit.2$coefficients[,where.line]))
                    }
                }
                else
                {
                    return(NULL)
                }
                
           })
    
    names(res.list) <- unique.matings
    
    null.list <- sapply(res.list, is.null)
    
    return(list(results=res.list[null.list == F], fit=fit.2, cont.counts=cont.counts))
}

#Creates a tabular summary of the differential expression results from de.bycross.to.mock or de.bycross.to.timecourse
#Input:
#       de.res.list: A list in the form of the 'results' element of de.bycross.to.mock or de.bycross.to.timecourse
#       direction.as.col: A logical indicating whether or not the Up or Down indicators should be their own columns or part of the Direction column
#       ...: Filters to apply to the data in the form: filter1=filter.vec, filter2=filter2.vec where filter.vec and filter2.vec are logical vectors named by probeset IDs.
#           TRUE should indicate that the probeset should be removed and FALSE should indicate that the value should be kept.
#Output:
#       A data.frame indicating the Cross (Mating), Comparison, Total or Filtered_total if direction.as.col == F or Up/Down if direction.as.col==T. 

tabular.summary.from.list <- function(de.res.list, direction.as.col=F, ...)
{
    
    flag.list <- list(...)
    
    ret.list <- mapply(function(val, name)
           {
                stopifnot("infecteds" %in% names(val))
                
                sum.dta <- as.data.frame(summary(val$infected))
                
                if (length(flag.list) > 0)
                {
                    
                    flag.mat <- sapply(flag.list, function(x)
                                       {
                                            x[rownames(val$infected)]
                                       })
                    
                    use.flags <- apply(flag.mat, 1, prod)
                    
                    new.sum.dta <- as.data.frame(summary(val$infected * (use.flags[rownames(val$infected)]==F)))
                    names(new.sum.dta)[3] <- "Filtered_Total"
                    sum.dta <- merge(sum.dta, new.sum.dta)
                    
                }
                
                sum.dta <- sum.dta[sum.dta$Var1 != 0,]
                sum.dta$Direction <- ifelse(sum.dta$Var1 == 1, "Up", "Down")
                sum.dta$Cross <- name
                sum.dta$Comparison <- gsub(paste0(name, "\\."), "", sum.dta$Var2)
                sum.dta$Total <- sum.dta$Freq
                
                sum.dta <- sum.dta[,names(sum.dta) %in% c("Var1", "Var2", "Freq") == F]
                
                if(length(flag.list) > 0)
                {
                    sum.dta <- sum.dta[,c("Cross", "Comparison", "Direction", "Filtered_Total")]
                }else{
                    sum.dta <- sum.dta[,c("Cross", "Comparison", "Direction", "Total")]
                }
                
                return(sum.dta)
                
           },de.res.list, names(de.res.list), SIMPLIFY=F)
    
    ret.dta <- do.call("rbind", ret.list)
    
    if (direction.as.col)
    {
        if (length(flag.list) > 0)
        {
            val.col <- "Filtered_Total"
        }else{
            val.col <- "Total"
        }
        
        ret.dta <- dcast(data=ret.dta, formula=Cross+Comparison~Direction, value.var=val.col)
    }
    
    return(ret.dta)
}


#Based on featureFilter from the genefilter package which seems to not work with these oligo expression sets, the point is to represent each gene by the probeset with the highest IQR.
#First the probesets without EntrezIDs are removed, then each probeset with the largest IQR for each Entrez ID is determined and used to subset the ExpressionSet.
#Input:
#   eset: ExpressionSet
#   annot.db: platform design database for the oligo arrays
#Output:
#   A new ExpressionSet with probesets representing a single gene
featureFilter <- function(eset, annot.db="mogene21sttranscriptcluster.db"){
    stopifnot(require(genefilter))
    stopifnot(require(annot.db,character.only=T))
    
    suppressWarnings(probe.ents <- select(get(annot.db), keys=featureNames(eset), columns=c("ENTREZID"), keytype="PROBEID"))
   
    use.probe.ents <- probe.ents[!is.na(probe.ents$ENTREZID),]
    
    sub.eset <- eset[unique(use.probe.ents$PROBEID),]
    
    unique.probes <- findLargest(featureNames(sub.eset), rowIQRs(sub.eset), annot.db)
    
    return(sub.eset[unique.probes,])
}

#Summarizes the result from de.bycross.to.mock
#Input:
#       res.list: A list as provided by the 'results' element of de.bycross.to.mock
#       by: Whether or not to annotated the 'summary.dta' data.frame with gene symbols or probeset IDs.
#       ...: Any filters to apply in the form: filter1=filter.vec, filter2=filter2.vec where filter.vec and filter2.vec are logical vectors named by probeset IDs.
#           TRUE should indicate that the probeset should be removed and FALSE should indicate that the value should be kept.
#           If not all the filters are FALSE for a given probeset, it will be removed from the 'summary.dta' element of the result, however it will only be flagged in the 'de.res' element.
#Output:
#       A list containing the following elements:
#           -de.res: A list named containing a data.frame containg the probeset, gene symbol, coeffecients (log2 scale) in the form: Inf_Timepoint.Inf_name.Mock_Timpoint.Mock_name.logFC
#                    and whether or not the comparison was significant named as before but with .Signif instead of .logFC.
#           -summary.dta: A data.frame containing either symbols or probeset IDs followed by a column or 1s and 0s for each Mating, indicating whether at least one contrast was significant over the timecourse
#                           after application of any supplied filters.

contrast.summary.tc <- function(res.list, by=c("symbol", "probe"), ...)
{
    #all rownames should be the same...
    
    by <- match.arg(by)
    
    filt.list <- list(...)
    
    annot.temp <- select(mogene21sttranscriptcluster.db, keys=rownames(res.list[[1]][[1]]), columns=c("SYMBOL"), keytype="PROBEID")
    stopifnot(anyDuplicated(annot.temp$SYMBOL) == 0)
    rownames(annot.temp) <- annot.temp$PROBEID
    
    cont.res.list <- lapply(res.list, function(x)
                            {
                                stopifnot(all(rownames(x$coefs) == rownames(x$infecteds)))
                                
                                coef.mat <- x$coefs
                                sig.mat <- x$infecteds
                                
                                any.sig <- apply(sig.mat, 1, function(x) any(x != 0))
                                
                                colnames(coef.mat) <- paste(gsub("M\\d+[xX]\\d+\\.", "", colnames(coef.mat)), "logFC", sep=" ")
                                colnames(sig.mat) <- paste(gsub("M\\d+[xX]\\d+\\.", "", colnames(sig.mat)), "Signif", sep=" ")
                                
                                temp.mat <- data.frame(ProbeId=rownames(coef.mat), Symbol=annot.temp[rownames(coef.mat),"SYMBOL"], coef.mat, sig.mat, stringsAsFactors=F)
                                
                                for(i in names(filt.list))
                                {
                                    temp.mat <- cbind(temp.mat, as.integer(filt.list[[i]][rownames(temp.mat)]))
                                    names(temp.mat)[ncol(temp.mat)] <- i
                                }
                                
                                return(temp.mat[any.sig,])
                            })
    
    #see if all columns are the same, if not add in the necessary columns with NAs
    
    cols <- lapply(cont.res.list, colnames)
    col.lens <- sapply(cols, length)
    col.max <- max(col.lens)
    
    exp.cols <- unique(unlist(cols))
    
    if (any(col.lens < col.max))
    {
        max.cols <- cols[[which.max(col.lens)]]
        
       for(i in which(col.lens < col.max))
        {
            diff.cols <- setdiff(max.cols, cols[[i]])
            fill.mat <- matrix(NA, ncol=length(diff.cols), nrow=nrow(cont.res.list[[i]]), dimnames=list(NULL, diff.cols))
            
            cont.res.list[[i]] <- cbind(cont.res.list[[i]], fill.mat)[,max.cols]
        }
    }
    
    cont.res.list <- cont.res.list[sapply(cont.res.list, nrow) > 0]
    
    cont.summary.dta <- data.frame(do.call("rbind", lapply(1:length(cont.res.list), function(x) cbind(cont.res.list[[x]], Mating=names(cont.res.list)[x]))))
    
    if (length(filt.list) > 0)
    {
        filt.mat <- apply(do.call("cbind", filt.list), 1, function(x) all(x == F))
        cont.summary.dta <- cont.summary.dta[cont.summary.dta$ProbeId %in% names(filt.mat)[filt.mat],]
    }
    
    use.form <- formula(paste0(switch(by, probe="ProbeId", symbol="Symbol"), "~Mating"))
    
    cont.summary <- dcast(use.form, data=cont.summary.dta, fill=0, value.var="Mating", fun.aggregate=length)
    
    return(list(de.res=cont.res.list, summary.dta=cont.summary))
}

#A function which produces timecourse images to accompany a ReportingTools report
#This should be run similar to: publish(result, cur.rep, cur.mating=cur.mating, eSet=use.exprs,limit=500, Tissue=cur.tissue, .modifyDF=list(timecourse.images))
#Note that cur.mating, eSet, limit and Tissue are used in this function in addition to the 'result' and 'cur.rep' arguments as defined below:
#cur.mating is the current Mating value
#eSet is the ExpressionSet used for the generation of the DE report
#limit is the limit on the number of images to produce assuming that df has been sorted prior to this function being called.
#Tissue is the tissue to limit the plot to.

timecourse.images <- function(df,...){
        
        rest.vars <- list(...)
        
        imagename <- character(nrow(df))
        
        if ('cur.mating' %in% names(rest.vars))
        {
            cur.mating <- sub("M", "", rest.vars$cur.mating)
        }else{
            cur.mating <- unique(unlist(regmatches(names(df), regexec("M(\\d+x\\d+)", names(df)))))[2]
        }
        
        if ('limit' %in% names(rest.vars))
        {
            max.rows <- rest.vars$limit
        }else{
            max.rows <- nrow(df)
        }
        
        stopifnot('eSet' %in% names(rest.vars))
        
        if ('Tissue' %in% names(rest.vars))
        {
            
            rest.vars$eSet <- rest.vars$eSet[,pData(rest.vars$eSet)$Tissue == rest.vars$Tissue]
        }
        
        use.dir <- file.path(dirname(path(rest.vars[[1]])), paste0("figures", cur.mating))
        if (file.exists(use.dir) == F)
        {
            dir.create(use.dir)
        }
        
        for (i in 1:min(nrow(df), max.rows)){
            imagename[i] <- paste0("plot_", df$ProbeId[i], ".png")
            use.sub <- rest.vars$eSet[df$ProbeId[i],grep(cur.mating, colnames(rest.vars$eSet))]
            cur.exprs <- cbind(pData(use.sub), exp=as.numeric(exprs(use.sub)))
            cur.exprs$Timepoint <- factor(cur.exprs$Timepoint, levels=order.times(cur.exprs$Timepoint), ordered=T)
            cur.exprs$Inf_Status <- sapply(as.character(cur.exprs$Inf_Status), function(x) switch(x, M="Mock", W="WNV", SV="SARS", IV="Flu"))
            cur.plot <- qplot(y=exp, x=Timepoint, data=cur.exprs, group=Inf_Status, color=Inf_Status, stat="summary", fun.y=mean, geom="line", ylab="Expression") + stat_summary(fun.data="mean_se", geom="errorbar", width=.2)
           
            ggsave(cur.plot, file=file.path(use.dir, imagename[i]), width=4, height=3, units="in")
        }
        
        df$Expression <- ifelse(imagename != "", hwriteImage(file.path(basename(use.dir), imagename), link=file.path(basename(use.dir), imagename), table=FALSE, width=100), "N/A")
        
        return(df)
    }

#A helper function which facilitates automatic figure/table numbering for ReportingTools reports
#Input:
    #report: An HTMLReportRef as created by HTMLReport in the ReportingTools package
    #node.type: Whether the report should be searched for images or tables present in the report.
#Output:
    #The current number of 'node.type' elements found in the report
report.entities <- function(report, node.type=c("img", "table"))
{
    search.str <- switch(node.type, img="a.img.src", table="table.tbody.tr.td")
    
    ent.count <- 0
    
    cur.pos <- 1
    while(inherits(tryCatch(report[[cur.pos]], error=function(e) NA), "XMLAbstractNode"))
    {
        if(any(names(unlist(xmlToList(report[[cur.pos]]))) == search.str))
        {
            ent.count <- ent.count + 1
        }
        
        cur.pos <- cur.pos+1
    }

    return(ent.count)
}

#A function to simplify adding a figure to an HTMLReportRef object.  Needs the 'convert' utility to be present.
#Input:
    #report: An HTMLReportRef as created by HTMLReport in the ReportingTools package
    #base.name: The prefix of the file as in make.boxplot
    #reports.dir: The dirname of the report to add the image to (in PNG form).  If not specified defaults to the main path of the report
    #plot.title: If specified, the title of the plot.  If not specified, it is set to a default if one currently exists
    #plot.text: If specified, the legend to add to the plot.  If not specified, it is set to a default if one currently exists
    #plot.type: A choice of one of the types of plot, assumed to be named as paste0(base.name, plot.type, '.png')
#Output:
    #An HTMLReportRef object with the specified figure inserted into it along with the title and legend.

add.image.to.report <- function(report, base.name, reports.dir=NULL, plot.title=NULL, plot.text=NULL, plot.type=c('raw', 'norm', 'heatmap', 'colored_heatmap', 'bac_spikes', 'polya_spikes'))
{
   plot.type <- match.arg(plot.type)
   figure.num <- report.entities(report, node.type="img") + 1
   
    if (missing(reports.dir) || is.null(reports.dir) || is.na(reports.dir))
    {
         out.path <- file.path(dirname(path(report)), paste0(base.name, plot.type, '.png'))
    }else{
         out.path <- file.path(reports.dir,paste0("QA_QC/",paste(base.name, plot.type,sep="_"), '.png'))
    }
   
   if (missing(plot.title) || is.null(plot.title))
    {
        plot.title <- switch(plot.type, raw="Overall Raw intensities",
                                    bac_spikes="Performance of the pre-labeled bacterial spike controls",
                                    colored_heatmap="Clustering of samples and genes based on experimental variables",
                                    heatmap = "Clustering of samples and genes based on experimental variables",
                                    norm="Boxplots of the overall intensity distribution are shown after RMA pre-processing",
                                    polya_spikes="Performance of the poly-adenylated spikes.")
    }
   
   if (missing(plot.text) || is.null(plot.text))
    {
        plot.text <- switch(plot.type, raw="(Top) Boxplot summaries of the expression intensity distribution of each of the samples colored by RIX cross and ordered by infection status and time within each cross.
                                        (Bottom) RIN scores for each sample with a cutoff indicated by the dashed line indicating whether the sample would be acceptable or not based on previous analyses. ",
                                    bac_spikes="Expression values for the bacterial spikes are shown.  All spikes should be observed and order should be as labeled in the legend. From the bottom: BioB, BioC, BioD, and then Cre.",
                                    colored_heatmap="Shown is a heatmap of the 1000 most variable genes.  On the bottom is a legend that indicates which timepoint (Timepoint=), RIX cross (Mating=) and infection status (Inf_status=)
                                    a given sample belongs to (black bar) as well as plots of the RIN score.  Coloring is applied to highlight outlier samples.",
                                    polya_spikes="Expression values for the poly-a spikes are shown, plotted by sense (st) or antisense (at) and location (3', 5' or middle (M)).  Rank ordering and expression for each category
                                    of spike is important with the ordering of Dap, Thr, Phe, Lys as indicated by position in the legend.",
                                    norm="")
    }
   
   
   
   system(paste("convert", paste0(reports.dir,"/QA_QC/",base.name, "_", plot.type, ".pdf"), out.path))
   himg <- hwriteImage(basename(out.path), link=basename(out.path), width=400)
   #himg <- hwriteImage(paste0(reports.dir,"/QA_QC/",basename(out.path)), link=basename(out.path), width=400)
   publish(hwrite(himg, br=TRUE), report)
   publish(hwrite(paste0('Figure ', figure.num, '. ', plot.title), heading=5), report)
   publish(hwrite(plot.text), report)
   return(report)
}


#A utility function, used in conjunction with .modifyDF in the publish method to add gene symbols to a data.frame already containing a ProbeId column
add.symbols <- function(df, ...){
    
        rest.vars <- list(...)
        annot.temp <- select(mogene21sttranscriptcluster.db, keys=df$ProbeId, columns=c("SYMBOL"), keytype="PROBEID")
        rownames(annot.temp) <- annot.temp$PROBEID
        use.df <- cbind(Symbol=annot.temp[df$ProbeId,"SYMBOL"], df)
        stopifnot(all(rownames(use.df) == df$ProbeId))
        #remove probeID
        #use.df <- use.df[,names(use.df) != "ProbeId"]
        return(use.df)
    }

#A utility function, used in conjunction with .modifyDF in the publish method to reorder and rename columns to be easier to view.
fix.contrasts.mock <- function(df,...){
        
        cur.mating <- unique(unlist(regmatches(names(df), regexec("M\\d+x\\d+", names(df)))))
        colnames(df) <- gsub(paste0(cur.mating, "\\."), "", colnames(df))
        #order them
        mating.cols <- colnames(df)[grep("D\\d+", colnames(df))]
        mating.ord <- order(as.numeric(sub("D", "", sapply(strsplit(mating.cols, "\\."), "[", 1))), decreasing=F)
        df <- df[,c("ProbeId", "Symbol", "Expression", mating.cols[mating.ord], "Adjusted p-Value")]
        names(df)[names(df) == "Adjusted p-Value"] <- "FDR Adjusted P-Value"
        
        return(df)
    }

#A utility function, used in conjunction with .modifyDF in the publish method to order the significance and log fold change columns by the day they were carried out
order.signif.lfc <- function(df,...)
{
    sig.df <- order.tc.cols(df,'\\.Signif')
    
    sig.lfc.df <- order.tc.cols(sig.df,'\\.logFC')
    
    return(sig.lfc.df)
}

#Orders a subset of the columns in a data.frame determined by the pattern in 'regex' by the experimental day assuming it is X in columns named by: DX.Y.Z
#Input:
    #df: The data.frame to be modified
    #regex: A regular expression to be used to locate the columns via grep.
#Output:
    #The modified data.frame
order.tc.cols <- function(df,regex, ...){
    
    signif.cols <- grep(regex, names(df))
    signif.names <- names(df)[signif.cols]
    
    init.cols <- 1:ncol(df)
    init.cols[init.cols %in% signif.cols] <- signif.cols[order(as.numeric(sub("D", "", sapply(strsplit(signif.names, "\\."), "[", 1))), decreasing=F)]
    
    return(df[,init.cols])
}

####Experimental, adapted or non-analysis functions

#adapted from ReportingTools:::.GOhyperG.to.html for use with the output from ClusterProfiler as opoosed to GOstats
go.df.compareCluster <- function(ck.sum, report.path="reports/GOPages")
{
    require(org.Mm.eg.db)
    
    ck.sum$GOLink <- paste("<a href=\"http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=",ck.sum$ID,"\">", ck.sum$ID,"</a>", sep = "")
    ck.sum$geneID <- as.character(ck.sum$geneID)
    ck.sum$goName <- ""
    #make the individual pages for the GO terms
    
    for (i in 1:nrow(ck.sum))
    {
        Oname <- ck.sum[i, "ID"]
        found.genes <- strsplit(ck.sum[i,"geneID"], "/")[[1]]
        countTable <- getNamesAndSymbols(found.genes, "org.Mm.eg.db")
        countTable$entrezLink <- paste("<a href=\"http://www.ncbi.nlm.nih.gov/gene/", 
            countTable$entrez, "\">", countTable$entrez, "</a>", 
            sep = "")
        table_title <- paste("Included genes in ", Oname, sep = "")
        GeneTable <- data.frame(countTable$entrezLink, countTable$symbol, 
            countTable$name)
        colnames(GeneTable) <- c("GeneEntrezId", "GeneSymbol", 
            "GeneName")
        pageFile <- as.character(strsplit(Oname, ":")[[1]][2])
        ck.sum[i,"goName"] <- pageFile
        Report <- HTMLReport(shortName = pageFile, title = table_title, 
            reportDirectory = report.path)
        publish(GeneTable, Report)
        finish(Report)
    }
    
    ck.sum$Count<- paste("<a href=\"", basename(report.path), "/", ck.sum$goName, 
        ".html", "\">", ck.sum$Count, "</a>", sep = "")
   
    ret <- data.frame(ck.sum$GOLink, ck.sum$Description, ck.sum$GeneRatio, ck.sum$BgRatio, ck.sum$Count, signif(ck.sum$pvalue, 3), signif(ck.sum$p.adjust, 3), stringsAsFactors = FALSE)
    colnames(ret) <- c("Accession", "GO Term", "GeneRatio", "BgRatio", "Overlap", "P-value", "BY-Adjusted P-value")
    
    return(ret)
}

#from ReportingTools:::getNamesAndSymbols
getNamesAndSymbols <- function (entrez, annotation.db) 
{
    require(annotate)
    geneids <- entrez
    entrezmap <- NULL
    tryCatch(entrezmap <- getAnnMap("ENTREZID", annotation.db), 
        error = function(e) {
        })
    if (is(entrezmap, "AnnDbBimap")) {
        entrez <- unlist(mget(geneids, entrezmap, ifnotfound = NA))
    }
    symbol <- unlist(mget(geneids, getAnnMap("SYMBOL", annotation.db), 
        ifnotfound = NA))
    name <- unlist(mget(geneids, getAnnMap("GENENAME", annotation.db), 
        ifnotfound = NA))
    countTable <- list(entrez = entrez, symbol = symbol, name = name)
    if (is(entrezmap, "AnnDbBimap")) {
        countTable$geneids <- geneids
    }
    return(countTable)
}


#to overwrite the version in ReportingTools and utilize .writeHTMLTable2
setMethod("objectToHTML", signature("data.frame"), function(object, report, .modifyDF, tableTitle = "",  filter.columns = sapply(object, is.numeric), colClasses = NULL, ...)
          {
            if (!missing(.modifyDF) && !is.null(.modifyDF)) {
            if (!is.list(.modifyDF)) 
                .modifyDF = list(.modifyDF)
            for (f in .modifyDF) object = f(object, report, object = object, ...)
            }
            if (nrow(object) == 0) 
                stop("No rows available in data.")
            if (ncol(object) == 0) 
                stop("No columns available in data.")
            if (is.null(colClasses)) {
                colClasses <- ReportingTools:::getDefaultColumnClasses(object, filter.columns = filter.columns)
            }
            col.specs <- data.frame(column = seq_along(object), label = colnames(object), 
                class = colClasses, stringsAsFactors = FALSE)
            numeric.columns <- which(unlist(lapply(object, class) == 
                "numeric"))
            for (col in numeric.columns) {
                object[, col] <- signif(object[, col], 3)
            }
            html = .writeHTMLTable2(object, tableTitle = tableTitle, 
                col.specs)
            list(html = html, object = object)
    })

#More or less copied from ReportingTools:::.writeHTMLTable, but modified for more interesting CSS...
.writeHTMLTable2 <- function(df, tableTitle, column.specs = NULL, p = NULL)
{
    df <- df[, column.specs$column, drop = FALSE]
    colnames(df) <- column.specs$label
    if (any(is.na(column.specs$class))) {
        column.specs[is.na(column.specs$class), ]$class <- ""
    }
    column.specs$class <- gsub("^\\s+", "", column.specs$class, 
        perl = TRUE)
    column.specs$class <- paste(column.specs$class, "top-header-row", 
        sep = " ")
    col.class <- data.frame(do.call(cbind, lapply(column.specs$class, 
        function(z) c(z, rep("", nrow(df))))), stringsAsFactors = FALSE)
    names(col.class) <- column.specs$label
   
    titleHtml <- hwrite(tableTitle, heading = 2)
    titleHtml <- sub("<.*?>", "<p class=\"page-header\">", titleHtml)
    titleHtml <- sub("</.*?>", "</p>", titleHtml)
    tableHtml <- hwrite(df, col.class = as.list(col.class), row.names = FALSE, 
        table.class = "dataTable table table-hover table-striped table-bordered")
    tableHtml <- sub("border=\"1\"", "", tableHtml)
    tableHtml <- sub("<tr>", "<thead><tr role=\"row\">", tableHtml)
    tableHtml <- sub("</tr>", "</tr></thead><tbody>", tableHtml)
    topTableHtml <- sub("(.*?)</thead>.*", "\\1", tableHtml)
    topHtml <- sub("(.*?)<thead>.*", "\\1", topTableHtml)
    topHeaderRow <- sub(".*<thead>(.*?)", "\\1", topTableHtml)
    topHtml <- sub("table ", "table cellpadding=\"0\" cellspacing=\"0\" border=\"0\"", 
        topHtml)
    topHeaderRow <- gsub("<td", "<th ", topHeaderRow)
    topHeaderRow <- gsub("</td>", "</th>", topHeaderRow)
    bottomHeaderRow <- gsub("top-header-row", "bottom-header-row sorting", 
        topHeaderRow)
    topHeaderRow <- gsub("top-header-row", "top-header-row no-print", 
        topHeaderRow)
    topHeaderRow <- gsub("(filter-num.*?top-header-row) no-print", 
        "\\1", topHeaderRow)
    headHtml <- paste(topHtml, "<thead>", topHeaderRow, bottomHeaderRow, 
        "</thead>", sep = "")
    bottomHtml <- sub(".*?</thead>(.*)", "\\1", tableHtml)
    bottomHtml <- sub("</table>", "</tbody></table>", bottomHtml)
    html <- paste("<div class=\"container\" style=\"margin-top: 10px\"> ", 
        titleHtml, headHtml, bottomHtml, "<foot>", bottomHeaderRow, 
        "</foot>", "</div>", sep = "")
    if (!missing(p)) {
        cat(html, file = p, sep = "\n")
        return(p)
    }
    else {
        return(html)
    }

}

#from toupper docs
capwords <- function(s, strict = FALSE) {
         cap <- function(s) paste(toupper(substring(s, 1, 1)),
                       {s <- substring(s, 2); if(strict) tolower(s) else s},
                                  sep = "", collapse = " " )
         sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
     }
     
#Experimental function meant to be more of a quick enrichment comparison function rather than a workhorse...
make.go.plot <- function(res.list)
{
    require(clusterProfiler)
    require(mogene21sttranscriptcluster.db)
    
    use.list <- lapply(res.list, function(x)
                           {
                                temp.list <- lapply(colnames(x$infecteds), function(y)
                                       {
                                            sig.x <- x$infecteds[,y]
                                            sig.x <- sig.x[sig.x != 0]
                                            
                                            if(length(sig.x) > 0)
                                            {
                                                annot.temp <- select(mogene21sttranscriptcluster.db, keys=names(sig.x), columns=c("ENTREZID"), keytype="PROBEID")
                                                return(unique(na.omit(annot.temp[,"ENTREZID"])))
                                            }
                                            else
                                            {
                                                return(character())
                                            }
                                            
                                       })
                                
                                names(temp.list) <- colnames(x$infecteds)
                                
                                return(temp.list)
                           })
        
    unl.list <- unlist(use.list, recursive=F)
    
    split.names <- strsplit(names(unl.list), "\\.")
    names(unl.list) <- mapply(function(name, mating)
                              {
                                    paste0(name, ".", gsub(paste0(name, "\\."), "", mating))
                              }, sapply(split.names, "[", 1), names(unl.list))
    
    ck <- compareCluster(unl.list, fun="enrichGO", organism="mouse", ont="BP", qvalueCutoff=1, minGSSize = 10)
    
    return (plot(ck) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
}

#Experimental function to make a heatmap ordered by time
#need to integrate tissue and lab information into this
make.timecourse.heatmap <- function(eset, use.genes=NULL, use.samples=NULL, plot.path=NULL)
{
    if (missing(use.genes) || is.null(use.genes))
    {
        use.genes <- featureNames(eset)
    }
    
    if (missing(use.samples) || is.null(use.samples))
    {
        use.samples <- sampleNames(eset)
    }
    
    sub.exprs <- exprs(eset[use.genes, use.samples])
    
    #order by infection status and mock then by day
    
    ord.dta <- data.frame(t(sapply(strsplit(colnames(sub.exprs), "_"), "[", 6:7)), stringsAsFactors=F)
    ord.dta$X1 <- factor(ord.dta$X1, levels=c(setdiff(ord.dta$X1, "M"), "M"), ordered=T)
    ord.dta$X2 <- factor(ord.dta$X2, levels=order.times(ord.dta$X2), ordered=T)
    
    col.ord <- do.call("order", ord.dta)
    
    sub.exprs <- sub.exprs[,col.ord]
    colnames(sub.exprs) <- sapply(strsplit(colnames(sub.exprs), "_"), function(x) paste(x[c(6,2,7)], collapse="_"))
    
    if (nrow(sub.exprs) > 500)
    {
        use.row.labels <- rep("", nrow(sub.exprs))
    }else{
        
        exprs.genes <- select(mogene21sttranscriptcluster.db, keys=rownames(sub.exprs), columns=c("SYMBOL"), keytype="PROBEID")
        rownames(exprs.genes) <- exprs.genes$PROBEID
        
        use.row.labels <- exprs.genes[rownames(sub.exprs),"SYMBOL"]
    }
    
    use.hm <- regHeatmap(sub.exprs, dendrogram = list(Col=list(status="no"), Row=list(status="yes")), labels=list(Row=list(labels=use.row.labels)), legend=T)
    
    if (missing(plot.path) || is.null(plot.path))
    {
        plot(use.hm)
    }else{
        png(width=800, height=800, file=plot.path)
        plot(use.hm)
        dev.off()
    }
    
}

#from genefilter:::rowIQRs, just bringing it into this namespace
rowIQRs <- function (eSet) 
{
    numSamp <- ncol(eSet)
    lowQ <- rowQ(eSet, floor(0.25 * numSamp))
    upQ <- rowQ(eSet, ceiling(0.75 * numSamp))
    upQ - lowQ
}