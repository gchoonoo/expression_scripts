#setwd("/Users/choonoo/expression_scripts")

source("affy_gene_qa_qc_GC.R")


library(plyr)

# prompt user to type in the folder name where analysis will be saved
reports.dir = gsub("/","_",gsub(" ","_",(readline("Please enter the name and date and of your study: "))))
cat(paste("Your analysis will be saved in the folder:", reports.dir))

#reports.dir <- "wnv_11_14_2014"

#this file proceeds in order of the functions.  The functions themselves are not executed per say as data.preparation() but instead line by line.
#they are defined as functions simply for organizational purposes.

data.preparation <- function()
{
    #Read in the data as a GeneFeatureSet and populate the pData slot with covariates to be utilized later
  
    # prompt user to enter path of expression data
    cel.file.dir = readline("Please enter the path to your expression data: ")
    
    cat(paste("The path of your expression data is:", cel.file.dir))
    
    #cel.file.dir="Gale_WNV_plate_11_14/"
    name.func=function(x) paste(x[1:7], collapse="_")
    name.col.char="_"
    
    affy.cels <- list.files(path=cel.file.dir, pattern=".CEL$", full.names=TRUE)
    
    #macq.cels <- affy.cels[grepl("MAQC", basename(affy.cels))]
    
  
    affy.cels <- affy.cels[! grepl("MAQC", basename(affy.cels))]
    
    
    cur.samples <- sub(".CEL", "", sapply(strsplit(basename(affy.cels), name.col.char), name.func))
    
################
# Chunk of code if you need to update your annotation file, (error would be thrown at line 94) so you can fix this part and then re-run
    
#     print(cur.samples)
#     
#     ## Fix sample annotation - missing lab info for some files
#     lab_annot = read.xls('./U19 Expression Archive/Line7_LabAnnotation.xlsx', sheet=1, header=T, as.is=T)
#     lab_annot = rbind(lab_annot, c('16441x8005', '121', 'L'))
#     lab_annot = rbind(lab_annot, c('16441x8005', '138', 'L'))
#     lab_annot = rbind(lab_annot, c('16441x8005', '139', 'L'))
#     
#     cur.pdta <- data.frame(do.call("rbind", strsplit(cur.samples, "_")), stringsAsFactors=FALSE)
#     cur.pdta[,1] = gsub('X', 'x', cur.pdta[,1])
#     for (i in 1:dim(cur.pdta)[1]) {
#       if (!cur.pdta[i,5] %in% c('L', 'G')) {
#         cur.pdta[i, 6:7] = cur.pdta[i,5:6]
#         ## Find Lab annotation for RIX_ID and update cur.pdta
#         cur.pdta[i, 5] = strsplit(lab_annot[lab_annot$Mating==cur.pdta[i,1] & lab_annot$RIX_ID==cur.pdta[i,2],3], "")[[1]][1]
#       }
#     }
#     
#     cur.samples = apply(cur.pdta, 1, function(x){paste(x[1:7], sep='_', collapse='_')})
#     print(cur.samples)

################

    #cur.samples_macq <- sub(".CEL", "", sapply(strsplit(basename(macq.cels), name.col.char), name.func))
    
    #raw.exprs <- read.celfiles(affy.cels, pkgname="pd.mogene.2.1.st", sampleNames=cur.samples_macq)
    
    raw.exprs <- read.celfiles(affy.cels, pkgname="pd.mogene.2.1.st", sampleNames=cur.samples)
      
    cur.pdta <- data.frame(do.call("rbind", strsplit(sampleNames(raw.exprs), "_")), stringsAsFactors=FALSE)
    
    
    names(cur.pdta) <- c("Mating", "Number", "Sex", "Tissue", "Lab", "Inf_Status", "Timepoint")
    rownames(cur.pdta) <- sampleNames(raw.exprs)
    cur.pdta$Sex <- toupper(cur.pdta$Sex)
    
    ####An example of how to populate the RIN data in the pData slot if available....
    #rna.dta <- read.delim("/Users/bottomly/Desktop/baric_heise_u19/Sample_data/datadumps/u19_datadump_20141103/RNA.txt",  sep="\t", header=TRUE, stringsAsFactors=FALSE)
    #stopifnot(all(rownames(cur.pdta) %in% rna.dta$Source.ID))
    #rownames(rna.dta) <- rna.dta$Source.ID
    #cur.pdta <- cbind(cur.pdta, rna.dta[rownames(cur.pdta),c("RIN..RIN...or.RQS.Score", "Instrument.for.QC", "Date.of.QC"), drop=FALSE])
    #names(cur.pdta)[(ncol(cur.pdta)-2):(ncol(cur.pdta))] <- c("RIN", "Instrument", "QC_Date")
    #cur.pdta$RIN <- as.numeric(cur.pdta$RIN)
    #
    ###end example
    
    pData(raw.exprs) <- cur.pdta[sampleNames(raw.exprs),]
    
    print(pData(raw.exprs))
    
    pdata_correct = readline("Is the annotation formated correctly? ")
        
    if (toupper(substring(pdata_correct,1,1)) != "Y"){
      stop("ERROR: Please check your annotation file is formatted correctly")
    }
    
    # save data from saved reports directory (GC)
    save(raw.exprs, file=paste0("full_raw_exprs_",reports.dir,".RData"))
    
    # commenting out remove work space and sourcing the file again since it loops through the prompts again if we do this

    #rm(list=ls())
    #source("wnv_analysis.R")
    #onto qc.plot
}



#Creates the QA/QC report.  The functions not defined in the ReportingTools package are defined in the file sourced at the top of the file.
qc.plot <- function()
{
    # load data from saved reports directory (GC)
    load(paste0("full_raw_exprs_",reports.dir,".RData"))
    #load("full_raw_exprs.RData")
    base.name <- reports.dir # save base names of pdfs as reports directory which should describe the study name and date
    
    if (file.exists(reports.dir) == F)
    {
      dir.create(reports.dir)
    }
    
    if (file.exists(paste(reports.dir,"QA_QC",sep="/")) == F)
    {
      dir.create(paste(reports.dir,"QA_QC",sep="/"))
    }
    
    make.boxplot(raw.exprs, base.name=base.name, type="raw", make.pdf=T)
    make.boxplot(raw.exprs, base.name=base.name, type="norm", make.pdf=T)
    
    make.heatmap(raw.exprs, base.name=base.name, make.pdf=T, include.tissue=T)
    #tissue kind of dominates all here, so will look at just the brain and spleen subsets
    
    #spleen
    make.heatmap(raw.exprs[,pData(raw.exprs)$Tissue == "Sp"], base.name=paste0(base.name, "_spleen"), make.pdf=T, include.tissue=T)
    
    #brain
    make.heatmap(raw.exprs[,pData(raw.exprs)$Tissue == "Br"], base.name=paste0(base.name, "_brain"), make.pdf=T, include.tissue=T)
    
    plot.bac.spikes(raw.exprs, base.name=base.name, make.pdf=T)
    plot.polya.spikes(raw.exprs, base.name=base.name, make.pdf=T)
    
    sum.df <- make.summary.table(raw.exprs)
    
    #then make the qa/qc report.
    
   
    
    # Prompt user to enter title of QA/QC (GC)
    title_sum = readline("Please enter the title for your QA/QC results or hit enter if same as name of the folder for your analysis: ")
    
    if(title_sum == ""){
      qa.qc.rep <- HTMLReport(shortName="qa_qc", title_sum=paste("QA/QC Results for Plate", reports.dir, sep=" "), reportDirectory=file.path(reports.dir, "QA_QC"))
      cat("QA/QC Title: ",paste("QA/QC Results for Plate", reports.dir, sep=" "))
    }else{
      qa.qc.rep <- HTMLReport(shortName="qa_qc", title_sum=paste("QA/QC Results for Plate", title_sum, sep=" "), reportDirectory=file.path(reports.dir, "QA_QC"))
      cat("QA/QC Title: ",paste("QA/QC Results for Plate", title_sum, sep=" "))  
    }
    
    # Prompt user to enter summary for QA/QC (GC)
    summary.text = readline("Please enter a summary for your QA/QC report, you can view the plots as PDF's in your folder: ")
    
    #summary.text <- ""
    
    publish(hwrite("Summary", br=T, heading=4), qa.qc.rep)
    publish(hwrite(summary.text), qa.qc.rep)
    
    #MAQC samples
    
    system(paste("convert /Users/bottomly/Desktop/baric_heise_u19/analyses/MAQC_raw_11_14_2014.pdf", file.path(dirname(path(qa.qc.rep)), "maqc.png")))
    himg <- hwriteImage("maqc.png", link="maqc.png", width=400)
    publish(hwrite(himg, br=TRUE), qa.qc.rep)
    publish(hwrite(paste0('Figure ', report.entities(qa.qc.rep, node.type="img"), '. ', "Raw intensity distributions for the MAQC control samples."), heading=5), qa.qc.rep)
    
    #image of polya spikes
    
    qa.qc.rep <- add.image.to.report(qa.qc.rep, base.name=base.name, plot.type="polya_spikes",reports.dir=reports.dir)
    
    #image of bac spikes
    
    qa.qc.rep <- add.image.to.report(qa.qc.rep, base.name=base.name, plot.type="bac_spikes",reports.dir=reports.dir)
    
    #raw boxplot
    
    qa.qc.rep <- add.image.to.report(qa.qc.rep, base.name=base.name, plot.type="raw", plot.text="Boxplot summaries of the expression intensity distribution of each of the samples colored by RIX cross and ordered by infection status, tissue and time within each cross.",reports.dir=reports.dir)
    
    #normalized boxplot
    
    qa.qc.rep <- add.image.to.report(qa.qc.rep, base.name=base.name, plot.type="norm",reports.dir=reports.dir)
    
    #heatmap
    
    heat.text <- "Shown is a heatmap of the 1000 most variable genes.  On the bottom is a legend that indicates which timepoint (Timepoint=), RIX cross (Mating=), tissue (Tissue=) and infection status (Inf_status=)
                                    a given sample belongs to (black bar)."
    
    qa.qc.rep <- add.image.to.report(qa.qc.rep, base.name=base.name, plot.type="heatmap", plot.text=heat.text,reports.dir=reports.dir)
    
    #spleen
    
    qa.qc.rep <- add.image.to.report(qa.qc.rep, paste0(base.name=base.name, "_spleen"), plot.type="heatmap", plot.title=paste("Clustering of samples and genes based on experimental variables limited to only the spleen samples"), plot.text="",reports.dir=reports.dir)
    
    #brain
    
    qa.qc.rep <- add.image.to.report(qa.qc.rep, paste0(base.name=base.name, "_brain"), plot.type="heatmap", plot.title=paste("Clustering of samples and genes based on experimental variables limited to only the brain samples"), plot.text="",reports.dir=reports.dir)
    
    #table 1
    publish(hwrite("Table 1.  Experimental Design after QA/QC procedure", heading=5), qa.qc.rep)
    
    publish(make.pretty.df(sum.df, rec.limit=1), qa.qc.rep)
    
    system(paste("open",finish(qa.qc.rep) ))
    
}

#This is how the analysis was performed. 
analysis <- function()
{   
    load(paste0("full_raw_exprs_",reports.dir,".RData"))
    
    use.exprs <- rma(raw.exprs, target="core")
    
    #it might be better to do this seperately for each tissue in the future, but for now...
    use.exprs <- featureFilter(use.exprs)
    
    #First we carry out the differential expression analysis seperately for Spleen and Brain returning a list of 
    
    # Prompt to get mock day for comparisons
    
    mock = readline("Would you like to make DE comparisons to closest mock or specified day? (type closest or use mock day) ")
    
    if(mock == 'closest'){
    
    html.res <- lapply(c("Sp", "Br"), function(x)
           {
                all.res.list <- de.bycross.to.mock(use.exprs, contrast.correction="global", subset.contrasts.by.tissue=x)
    
                return(contrast.summary.tc(all.res.list$results))
           })
    
    names(html.res) <- c("Sp", "Br")
    
    #run it again, but this time keep only the initial results not just the summary
    tissue.res.list <- lapply(c("Sp", "Br"), function(x)
                              {
                                return(de.bycross.to.mock(use.exprs, contrast.correction="global", subset.contrasts.by.tissue=x))
                              })
    }else if(mock == 'use mock day'){
      mock_use = readline("Please enter your mock day (Ex. D12): ")
      
      html.res <- lapply(c("Sp", "Br"), function(x)
      {
        all.res.list <- de.bycross.to.mock(use.exprs, contrast.correction="global", subset.contrasts.by.tissue=x, use.mock.time=as.character(mock_use))
        
        return(contrast.summary.tc(all.res.list$results))
      })
      
      names(html.res) <- c("Sp", "Br")
      
      #run it again, but this time keep only the initial results not just the summary
      tissue.res.list <- lapply(c("Sp", "Br"), function(x)
      {
        return(de.bycross.to.mock(use.exprs, contrast.correction="global", subset.contrasts.by.tissue=x, use.mock.time=as.character(mock_use)))
      })
    }else{
      stop("ERROR: Please re-enter your mock comparisons")
    }
    
    names(tissue.res.list) <- c("Sp", "Br")
    
    if (file.exists(reports.dir) == F)
    {
        #unlink(file.path(reports.dir, "GO_Result"), recursive=T)
        
        dir.create(reports.dir)
    }
    
    #GO results
    
    univ <- select(mogene21sttranscriptcluster.db, keys=featureNames(use.exprs), columns=c("ENTREZID"), keytype="PROBEID")
    
    #summarize per tissue with the end result being: list(tissue=list(comparison=go.result))
    go.prof <- lapply(html.res, function(x)
                      {
                        #here x is the summarized result from the mock comparisions see 'contrast.summary.tc'.
                        #we first figure out what all the significance columns are based on those from the first Mating
                        #and assuming they are the same, which should be safe here.
                        signif.cols <- colnames(x$de.res[[1]])[grep("Signif", colnames(x$de.res[[1]]))]
                        
                        #We then iterate over all of the mock comparisons and for each one extract the significant results for each Mating in a list for basic.go.enrich.tc
                        de.res <- lapply(signif.cols, function(y)
                                    {
                                        m.list <- lapply(x$de.res, function(z)
                                            {
                                                return(z[is.na(z[,y]) == F & z[,y] != 0,])
                                            })
                                        
                                        return(basic.go.enrich.tc(m.list, univ))
                                    })
                        
                        names(de.res) <- signif.cols
                        return(de.res)
                      })
    
    #will go with the first categories that give signifiant GO terms:
    print("Significant GO comparisons for Brain")
    print(sapply(go.prof$Br, function(x) nrow(summary(x))))
    #D7.W.D2.M.Signif  D4.W.D2.M.Signif D12.W.D2.M.Signif 
    #               0                 0                72
    
    print("Significant GO comparisons for Spleen")
    print(sapply(go.prof$Sp, function(x) nrow(summary(x))))
    #D7.W.D2.M.Signif  D4.W.D2.M.Signif  D2.W.D2.M.Signif D12.W.D2.M.Signif 
     #         324               610               526               189 
    
    #so Day 2 for Spleen and Day 12 for Brain...
 
    # prompt user to enter sig GO comparisons to use (GC)
    use.tissue.brain = readline("Please enter the GO category with the most significant terms for brain: ")
 
    use.tissue.spleen = readline("Please enter the GO category with the most significant terms for spleen: ")
 
    use.tissue <- c(Br=as.character(use.tissue.brain), Sp=as.character(use.tissue.spleen))
    
    #use.tissue <- c(Br="D12.W.D2.M.Signif", Sp="D2.W.D2.M.Signif")
    
    #make the report pages for these two tissues
    go.tissue <- lapply(names(go.prof), function(tissue)
                        {
                            go.sum <- summary(go.prof[[tissue]][[use.tissue[tissue]]])
                            split.go <- split(go.sum, go.sum$Cluster)
                            
                            pretty.time <- sub("D", "Day ", strsplit(use.tissue[tissue], "\\.")[[1]][1])
                            pretty.tissue <- switch(tissue, Br="Brain", Sp="Spleen")
                            
                            message(paste(pretty.time, pretty.tissue))
                            
                            go.list <- mapply(function(cur.df, cur.name)
                                              {
                                                    message(cur.name)
                                                    
                                                    go.rep <- HTMLReport(shortName=paste0("GO_", cur.name, "_", tissue),
                                                                         title=paste(pretty.tissue, pretty.time, "Infected-Mock Gene Ontology Enrichment Results for Cross", cur.name),
                                                                         reportDirectory=file.path(reports.dir, "GO_Results"))
                                                    
                                                    temp.df <- go.df.compareCluster(cur.df,
                                                                                    report.path=file.path(dirname(path(go.rep)),paste0("GOPages_", cur.name, "_", tissue)))
                                                    publish(temp.df, go.rep)
                                                    finish(go.rep)
                                                    return(go.rep)
                                              }, split.go, names(split.go))
                            
                            return(go.list)
                        })
    
    names(go.tissue) <- names(go.prof)
    
    #Text summaries
    
    summary.rep <- HTMLReport(shortName="summary", title="Overall Differential Expression Summary", reportDirectory=file.path(reports.dir, "DE_Summary"))
    publish(hwrite('Differential expression was assessed per timepoint/mock comparison for each tissue (BH FDR < .05; adjusted globally within cross)'), summary.rep)
   
    tab.sum.br <- tabular.summary.from.list(tissue.res.list$Br$results, direction.as.col=F)
    tab.sum.sp <- tabular.summary.from.list(tissue.res.list$Sp$results, direction.as.col=F)
    
    all.tab.sum <- rbind(cbind(tab.sum.br, Tissue="Brain"), cbind(tab.sum.sp, Tissue="Spleen"))
    
    de.sum.df <- dcast(Cross+Comparison~Tissue+Direction, data=all.tab.sum, value.var="Total")
    
    de.sum.df.ord <- ddply(.data=de.sum.df, .variables="Cross", .fun=function(df)
                           {
                                base.times <- df$Comparison
                                day.only <- sapply(strsplit(df$Comparison, "\\."), "[", 1)
                                names(base.times) <- day.only
                                df$Comparison <- factor(df$Comparison, levels=base.times[order.times(day.only)], ordered=T)
                                
                                df[order(df$Comparison),]
                           })
    
    publish(hwrite('Table of DE Results', heading=5), summary.rep)
    publish(make.pretty.df(de.sum.df.ord, sep.header=T), summary.rep)
    
    
    #then a heatmap summary of DE genes
    summary.mat.br <- html.res$Br$summary.dta
    colnames(summary.mat.br)[-1] <- paste("Brain", colnames(summary.mat.br)[-1], sep="_")
    summary.mat.sp <- html.res$Sp$summary.dta
    colnames(summary.mat.sp)[-1] <- paste("Spleen", colnames(summary.mat.sp)[-1], sep="_")
    
    summary.mat <- merge(summary.mat.br, summary.mat.sp, by="Symbol", all=T, incomparables=NA)
    summary.mat[is.na(summary.mat)] <- 0
    
    rownames(summary.mat) <- summary.mat$Symbol
    summary.mat <- as.matrix(summary.mat[,-1])
    
    binary.dist <- function(x) dist(x, method="binary")
    
    use.summary.mat <- summary.mat[order(rowSums(summary.mat), decreasing=T),][1:min(1500, nrow(summary.mat)),]
    
    binary.hm <- regHeatmap(use.summary.mat, dendrogram = list(distfun=binary.dist),
                            labels=list(Col=list(nrow=15), Row=list(labels=rep("", nrow(use.summary.mat)))), scale="none", legend=3, breaks=2)#,
    
    png(file=file.path(dirname(path(summary.rep)),"de_summary_heat.png"), width=10, height=9, units="in", res=300)
    plot(binary.hm)
    dev.off()
    himg <- hwriteImage('de_summary_heat.png', link='de_summary_heat.png', width=400)
    publish(hwrite(paste('Clustering of DE Results by Cross and Gene (top',nrow(use.summary.mat),'genes ranked by number of crosses)'), heading=5), summary.rep)
    publish(hwrite(himg, br=TRUE), summary.rep)
    
    #then make a link to a searchable summary matrix
    summary.rep.mat <- HTMLReport(shortName="summary_matrix", title="Explore Heatmap", reportDirectory=dirname(path(summary.rep)))
    publish(summary.mat, summary.rep.mat)
    finish(summary.rep.mat)
    
    publish(Link(summary.rep.mat, report=summary.rep), summary.rep)
    
    #GO Summary
    #use.tissue looks like c(Br="D12.W.D2.M.Signif", Sp="D2.W.D2.M.Signif")
    
    clust.res <- data.frame()
    gene.clusts <- list()
    
    for(i in names(use.tissue))
    {
        temp.comp.dta <- go.prof[[i]][[use.tissue[i]]]@compareClusterResult
        temp.comp.dta$Cluster <- paste(temp.comp.dta$Cluster, switch(i, Br="Brain", Sp="Spleen"), sapply(strsplit(use.tissue[i], "\\."), "[", 1),sep="_")
        clust.res <- rbind(clust.res, temp.comp.dta)
        
        temp.gene.clusts <- go.prof[[i]][[use.tissue[i]]]@geneClusters
        names(temp.gene.clusts) <- paste(names(temp.gene.clusts), switch(i, Br="Brain", Sp="Spleen"), sapply(strsplit(use.tissue[i], "\\."), "[", 1) ,sep="_")
        
        gene.clusts <- append(gene.clusts, temp.gene.clusts)
    }
    
    #a hack to improve how the categories are displayed..., plot the top 15 categories in terms of adjusted p-value instead of ranking them by count...
    
    new.clust.res <- ddply(.data=clust.res, .variables=.(Cluster), .fun=function(df, showCategory=15)
                {
                    return(df[order(df$p.adjust, decreasing=F),][1:min(showCategory, nrow(df)),])
                })
    
    test.class <- new("compareClusterResult", compareClusterResult=new.clust.res, fun=enrichGO, geneClusters=gene.clusts )
    
    validObject(test.class)
    
    go.obj.plot <- plot(test.class, showCategory=15)
    ggsave(go.obj.plot, file=file.path(dirname(path(summary.rep)), 'go_summary.png'), width=16, height=8)
    himg <- hwriteImage('go_summary.png', link='go_summary.png', width=400)
    publish(hwrite('Summary of Biological Process GO Categories for Infected-Mock comparison', heading=5), summary.rep)
    publish(hwrite('GeneRatio: Number of DE genes in a given category / Total DE genes assigned to a GO term'), summary.rep)
    publish(hwrite(himg, br=TRUE), summary.rep)
    
    finish(summary.rep)
    
    #DE reports
    res.list <- mapply(function(cur.tissue, de.res)
                       {
                            message(cur.tissue)
                            mapply(function(cur.mating, result)
                            {
                                message(cur.mating)
                               
                                max.fc.row <- apply(result[,grep("logFC", names(result))], 1, function(x) max(abs(x)))
                                
                                result <- result[order(max.fc.row, decreasing=T),]
                                
                                pretty.tissue <- switch(cur.tissue, Br="Brain", Sp="Spleen")
                                
                                cur.rep <- HTMLReport(shortName=paste(cur.mating, pretty.tissue, sep="_"), title=paste("Differential Expression Results for Cross", sub("M", "", cur.mating), pretty.tissue), reportDirectory=file.path(reports.dir, "DE_Results"))
                                publish(hwrite('Gene Expression figures are shown for the 500 genes ranked by absolute foldchange'), cur.rep)
                                
                                pub.df <- publish(result, cur.rep, cur.mating=cur.mating, eSet=use.exprs,limit=500, Tissue=cur.tissue, .modifyDF=list(timecourse.images, order.signif.lfc))
                                finish(cur.rep)
                                return(cur.rep)
                            }, names(de.res$de.res), de.res$de.res)
                            
                       }, names(html.res), html.res)
    
    ##make the index page
    index.page <- HTMLReport(shortName = "index",title = paste(reports.dir,"Expression Results", sep=" "), reportDirectory = reports.dir)
    publish(Link("Quality Control Summary", target=file.path(reports.dir, "QA_QC" ,"qa_qc.html"), report=index.page), index.page)
    
    publish(hwrite(hmakeTag('br')), index.page)
    publish(Link(summary.rep, report=index.page), index.page)
    
    publish(hwrite(hmakeTag('br')), index.page)
    publish(Link(res.list$Sp, report = index.page), index.page)
    publish(Link(res.list$Br, report = index.page), index.page)
    
    publish(hwrite(hmakeTag('br')), index.page)
    publish(Link(go.tissue$Sp, report=index.page), index.page)
    publish(Link(go.tissue$Br, report=index.page), index.page)
    
    finish(index.page)
    
}

data.preparation()
qc.plot()
analysis()

if(interactive()) x11()
