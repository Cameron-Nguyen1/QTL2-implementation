#QTL_Mapping_R
deps = c("datapasta","jsonlite","qtl2","qtl","qtl2convert","paletteer","grid","gridExtra","ggplot2","dplyr","qtl2ggplot", "biomaRt")
lapply(deps,library,character.only=TRUE)
set.seed(8675309)

#Define QTL plot function
plot_qtl = function(scan1,perm_sum,map,pdf_name){
    myl = list()
    z = 1
    print(perm_sum)
    for (i in 1:length(colnames(scan1))){
        p1 = ggplot_scan1(scan1,map=map,lodcolumn=i,col=paletteer_d("ggthemes::calc")[i]) + 
        geom_hline(yintercept=perm_sum[z],linetype="solid",color='yellow3',size=1.1)+
        geom_hline(yintercept=perm_sum[z+1],linetype="solid",color='orange1',size=1.1)+
        geom_hline(yintercept=perm_sum[z+2],linetype="solid",color='red3',size=1.1)+
        ggtitle(colnames(scan1)[i])+theme(plot.title = element_text(hjust = 0.5))
        #p1 <- plot(dat.out,lodcolumn = i, marker.map, col=paletteer_d("ggthemes::calc")[i]) + title(colnames(dat.out)[i])
        myl[[colnames(scan1)[i]]] = p1
        z = z+3
    }
    g=grid.arrange(arrangeGrob(grobs=myl,ncol=2))
    ggsave(paste0(pdf_name,".pdf"),plot=g,device='pdf',width=24,height=24)
}

#Read cross and run analysis
cross = read_cross2("MYCONTROL.yaml")  # <-- Replace MYCONTROL with your file.

map = insert_pseudomarkers(cross$gmap, step=1)
pr = calc_genoprob(cross, map, error_prob=0.002)
apr = genoprob_to_alleleprob(pr)
kinship = calc_kinship(pr)
Xcovar = get_x_covar(cross)
out = scan1(pr, cross$pheno, Xcovar=Xcovar, cores=4)
operm = scan1perm(pr, cross$pheno, Xcovar=Xcovar, n_perm=10000, cores=10)
operm_sum = summary(operm,alpha=c(0.2, 0.1, 0.05))
out_peaks = find_peaks(out, cross$gmap, threshold=min(operm_sum))

#Plot QTL, merge and format data out
plot_qtl(out,operm_sum,cross$gmap,"MYBACKCROSS") # <-- Name the output.
z = data.frame(t(operm_sum))
colnames(z) = c("Alpha_0.2","Alpha_0.1","Alpha_0.05")
z$rownames = rownames(z)
z_merge = merge(out_peaks, z, by.x = "lodcolumn", by.y = "rownames", all.x = TRUE) %>% mutate("Significance"= lod >= Alpha_0.05)
z_merge = z_merge[,-2]
colnames(z_merge) = c("Phenotype","Chromosome","Position","LOD",colnames(z_merge[5:length(colnames(z_merge))]))
write.csv(z_merge,file="QTL_Peaks.csv",row.names=FALSE)

#Query QTL peak bandwidths
to_query = unlist(z_merge[z_merge$Significance==TRUE,][,c(1,2)])
for (i in 1:(length(to_query)/2)){
    if (length(to_query)/2==1){
        lodc = to_query[["Phenotype"]]
        CHR = to_query[["Chromosome"]]
    }else{
        lodc = to_query[[paste0("Phenotype",as.character(i))]]
        CHR = to_query[[paste0("Chromosome",as.character(i))]]
    }
    x = qtl2::bayes_int(out,cross$gmap,
        lodcolumn=lodc,  ### Significant Phenotype LOD peak i
        chr=as.numeric(CHR),  ### Significant Chr belonging to SP_LOD peak i
        prob=0.95)
    POS1 = x[1,1]*1000000
    POS2 = x[1,3]*1000000

    ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    genes = getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'start_position', 'end_position', 'strand', 'description'), 
        filters = c('chromosome_name', 'start', 'end'), 
        values = list(CHR, POS1, POS2), mart = ensembl) %>% arrange(external_gene_name=="",external_gene_name)
    genes = genes %>% mutate("Phenotype"=rep(lodc,nrow(genes)),"Chromosome"=rep(CHR,nrow(genes))) %>% relocate("Phenotype",.before=1) %>% relocate("Chromosome",.before=2)
    if (file.exists("QTL_Gene_Info.csv") == TRUE){
        write.table(genes,file="QTL_Gene_Info.csv",row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }else{
        write.csv(genes,file="QTL_Gene_Info.csv",row.names=FALSE)
    }
}