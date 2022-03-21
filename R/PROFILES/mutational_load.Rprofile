.First<-function(){
    #Initiating profile
    cat("Successfully loaded profile for germline analysis \n")
}

#Setting paths to libraries, rm if not using uppmax.
#    .libPaths(c("/home/sesagger/R/packages", .libPaths()))

#Make new environment.
    .env<-new.env()

#Loading necessary packages
    .env$packages <- c("tidyverse","gdata","utils","data.table","ggrepel","RColorBrewer")

    .env$libraries<-function(){
    sapply(.env$packages,require,character=TRUE)
    }   

    .env$libraries()


#Div

#Input stuff
    .env$lavfil <- function(df,filnavn,separator){
        if(missing(separator)){
            separator="\t"
        }
        fwrite(df, file = filnavn, append = FALSE, quote = "auto", sep = separator,row.names = FALSE,)
    }

    .env$hentfil<-function(filnavn,varnavn,separator="\t",head=TRUE){
        if(separator=="plink"){
            separator=""
        }
        data<-fread(filnavn,header=head,stringsAsFactors=FALSE,sep=separator)
        assign(paste(deparse(substitute(varnavn))),data,envir=.GlobalEnv)
        }


#Data handling stuff
    .env$clean_plink <- function(df){
        df1<-df
            # Makes new cleaned-up df based on df.
            # chooses all cases with non-missing values (should be removed before R, but just in case)
        df1$Z_FST<-scale(df1$FST)
        com_df<-df1[complete.cases(df1),];
        newdf<-data.frame(com_df)
        newdf$Location<-str_replace(newdf$SNP,"_",":")
        assign(paste(deparse(substitute(df)),"clean",sep='_'),newdf,envir = .GlobalEnv)
    }

    .env$clean_vt <- function(df){
        df<-rename(df,"FST"="WEIR_AND_COCKERHAM_FST")
        df<-subset(df,!is.na(df$FST))
        df_sig<-df[complete.cases(df),];
        SNP<-c(1:(nrow(df_sig)))
        newdf<-data.frame(SNP,df_sig)
        newdf<-unite(newdf, Location, c(CHROM, POS), sep=":",remove=FALSE)
        newdf$Location<-as.factor(newdf$Location)
        # levels(newdfsig)[levels(newdfsig)=="chrX"]<-'39'
        # levels(newdfsig)[levels(newdfsig)=="chrM"]<-'40'
        # newdfsig$CHR<-as.factor(str_replace(newdfsig$CHROM,'chr',''))
        assign(paste(deparse(substitute(df)),"clean",sep='_'),newdf,envir = .GlobalEnv)
    }


#Statistics stuff
    .env$hp_calc<-function (df, bin_size, bin_step){
    df3 <- c()
    j <- 1
    print("Binning starting, Do you have $C1,$C2 and $G0? Should be done on plink.frq.counts files")
    while (j < 39) {
        a <- subset(df, CHR == j)
        if (dim(a)[1] != 0) {
            i <- 0
            while (i < max(a$POS)) {
                bin_start <- i + 1
                bin_end <- i + bin_size
                b <- subset(a, a$POS >= bin_start & a$POS <= bin_end)
                C1_sum <- sum(b$C1, na.rm = TRUE)
                print(C1_sum)
                C2_sum <- sum(b$C2, na.rm = TRUE)
                G0_sum <- sum(b$G0, na.rm = TRUE)
                hp_1 <- 2 * C1_sum * C2_sum/(C1_sum + C2_sum)^2
                if (C1_sum != 0 | C2_sum != 0 | G0_sum != 0) {
                  df2 <- cbind(j, bin_start, bin_end, C1_sum, 
                    C2_sum, G0_sum, hp_1)
                  df3 <- bind_rows(df3, as.data.frame(df2))
                }
                i <- i + bin_step
                print("Finished",j,bin_start,bin_end)
            }
        }
        print(paste0("chr ", j, " done"))
        j <- j + 1
    }
    print("Summing beginning")
    df3$hp <- 2 * df3$C1_sum * df3$C2_sum/(df3$C1_sum + df3$C2_sum)^2
    print("Z-transforming beginning")
    df3$Z_hp <- scale(df3$hp_1, center = TRUE, scale = TRUE)
    print("Naming beginning")
    binsi <- bin_size/1000
    binst <- bin_step/1000
    assign(paste(deparse(substitute(df)), "w", paste(deparse(substitute(binsi))), 
        "s", paste(deparse(substitute(binst))), "summed", sep = "_"), 
        df3, envir = .GlobalEnv)
    summary(df3)
}
#Summation stuff
    .env$pull_tables<-function(df){
        print("Deleterious mutations")
        subset(df,grepl("elet",Extra)&FCR!=0) %>% dplyr::select("CHR","SNP","Gene","Location","POS","NMISS","FST","Z_FST","#Uploaded_variation","Allele","Consequence","Amino_acids","Codons","Extra","external_gene_name","A1","FCR","Eriks") %>% distinct(Location,Gene,.keep_all=TRUE) %>% print(n=100)

        print("Missense mutations")
        subset(df,grepl("missense",Consequence)&FCR!=0) %>% dplyr::select("CHR","SNP","Gene","Location","POS",
            "NMISS","FST","Z_FST","#Uploaded_variation","Allele","Consequence","Amino_acids","Codons","Extra","external_gene_name","A1","FCR","Eriks") %>% distinct(Location,Gene,.keep_all=TRUE) %>% print(n=100)

        print("Highest Z_FST")
        df %>% arrange(desc(Z_FST)) %>% dplyr::select("CHR","SNP","Gene","Location","POS","NMISS","FST","Z_FST","#Uploaded_variation","Allele","Consequence","Amino_acids","Codons","Extra","external_gene_name","A1","FCR","Eriks")  %>% distinct(Location,Gene,.keep_all=TRUE)%>% head(n=10)

    }

#Annotation stuff
    .env$annotategenes<-function(df){
        #You need to have the object ann loaded, this should be a VEP annotation file
        df_ann<-merge(df,ann,all.x=TRUE,by="Location")
        assign(paste(deparse(substitute(df)),"ann",sep='_'),df_ann,envir = .GlobalEnv)
    }

    .env$relevant<-function(df,fst){
        #Subsetting FST to the given level.
        newdf<-subset(df,Z_FST>=fst)
        assign(paste(deparse(substitute(df)),deparse(substitute(fst)),"fst",sep='_'),newdf,envir = .GlobalEnv)
    }

    .env$addgenes<-function(df){
        ensembl<-useEnsembl(biomart = 'genes', dataset = 'cfamiliaris_gene_ensembl', version = 99)
        df1<-subset(df,Gene!="-")
        df2<-df1 %>% count(Gene)
        print("annotation starting")
        genes<-getBM(c('ensembl_gene_id','external_gene_name','start_position','end_position'),filter='ensembl_gene_id',values=df2$Gene,mart=ensembl)
        print("Start merging")
        newdf<-merge(df,genes,by.x="Gene", by.y="ensembl_gene_id",all.x=TRUE)
        newdf$external_gene_name<-as.factor(newdf$external_gene_name)
        assign(paste(deparse(substitute(df)),"genes",sep='_'),newdf,envir = .GlobalEnv)
    }

    .env$addsnp<-function(df){
        snp<-useMart("ENSEMBL_MART_SNP",dataset="clfamiliaris_snp")
        nydata <- matrix(ncol=4, nrow=dim(df)[1])
        nydata<-as.data.frame(nydata)
        names(nydata)<-c('refsnp_id','allele','chr_name','chrom_start')
        udtraek1<-c("empty","empty","empty","empty")
        for (i in dim(df)[1]){
            nydata[i,]<-getBM( c("refsnp_id", "allele", "chr_name", "chrom_start"), filters = c("chr_name", "start", "end"), values = list(df$CHR[i], df$bp[i], df$bp[i]), mart = .env$snp)

            assign(paste(deparse(substitute(df)), "snps", sep = "_"), 
            nydata, envir = .GlobalEnv)
        }
    }

#Figures
    .env$manhat<-function(df,jpgname,afgraensx,afgraensy,snpofint){
        #It works even when warning shows up, but the missing snps will not be highlighted. Quick and dirty, not usefull for publishing
        if((missing(afgraensx))){
            afgraensx<-c(0,max(df$POS))
        }
        if((missing(afgraensy))){
            afgraensy<-c(0,1)
        }   
        if((missing(snpofint))){
            snpofint<-FALSE
        }
        newdfsig <- subset(df, !grepl("Un", df$CHR) & FST >0)
        newdfsig$CHR <- as.numeric(newdfsig$CHR)
        newdfsig$POS <- as.numeric(newdfsig$POS)
        jpeg(jpgname)
            manhattan(newdfsig, chr = "CHR", bp = "POS", p = "FST", logp = FALSE, xlim=afgraensx,
            ylab = "Weir and Cockerham Fst",ylim=afgraensy,highlight=snpofint)
        dev.off()
    }

.env$gg.manhattan_pos <- function(df1,chr , threshold,  ylims, title,hlight,col){
    #https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/
    #GGplot2-Manhattan-Plot-Function fixed for FST
    #Needed packages
    # library(readr)
    # library(ggrepel)
    # library(ggplot2)
    # library(dplyr)
    # library(RColorBrewer)
    if(is.data.frame(df1)==FALSE){
        print("Input must be a data.frame")
    }
    else if(is.data.table(df1)==TRUE){
        print("Input must be data.frame, not data.table, use as.data.frame()")
    }else{
        if((missing(threshold))){
            threshold<-5
        }
        if((missing(col))){
            col<-c("#5D82BB", "#3B64A5", "#1E4F9E", "#103B7E", "#082B64"
                )
        }
        if((missing(title))){
            title<-"Unknown"
        }  
        if((missing(hlight))){
            hlight<-NULL
        } 
         
        # format df
        df<- subset(df1,Z_FST>0&CHR==chr)
        df.tmp <- df %>%

        # Compute chromosome size
        group_by(CHR) %>% 
        summarise(chr_len=max(POS)) %>%

        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
        dplyr::select(-chr_len) %>%

        # Add this info to the initial dataset
        left_join(df, ., by=c("CHR"="CHR")) %>%

        # Add a cumulative position of each SNP
        arrange(CHR, POS) %>%
        mutate( BPcum=POS+tot) %>%

        # Add highlight and annotation information
        #mutate( is_highlight=ifelse(SNP %in% hlight, "Non-coding transcript variant", "Other")) %>%
        mutate( is_annotate=ifelse(Z_FST > 5, "yes", "no")) %>%
        mutate( is_annotate1=ifelse(Z_FST > 5.5, "yes", "no"))


        # get chromosome center positions for x-axis
        axisdf <- df.tmp %>% group_by(POS,Consequence) %>% summarize( #Add POS here if you want to look single CHR
        center=( max(BPcum) + min(BPcum) ) / 2 )

        ggplot(df.tmp, aes(x=BPcum, y=Z_FST)) +
        # Show all points
        geom_point(aes(color=as.factor(CHR)),  size=1) +
        scale_color_manual(values = c("#5D82BB","red"),name="Consequence") +

        # custom X axis: 
        #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) + #Remove for single CHR (i think)
        scale_y_continuous(expand = c(0, 0), limits = ylims) +  #expand=c(0,0) #removes space between plot area and x axis 

        # add plot and axis titles
        ggtitle(title) +
        labs(x = "Position",y=expression(ZF["ST"])) +#CHange to "Position" if looking at single CHR

        # add genome-wide sig and sugg lines
        geom_hline(yintercept = 5.5,color="orange")+
        geom_hline(yintercept = 5,color="orange",linetype="dashed")+

        # Add highlighted points
        geom_point(data=subset(df.tmp, is_annotate=="yes"), color="orange", size=1) +
        geom_point(data=subset(df.tmp, is_annotate1=="yes"), color="orangered", size=2) +

        #geom_point(data=subset(df.tmp, is_highlight=="Non-coding transcript variant"), color="red",size=3) +

        # Add label using ggrepel to avoid overlapping, I usally don't use the labels so I've just commented it out
        #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +


        # Custom the theme:
        theme_bw(base_size = 22) +
        theme(                  
          plot.title = element_text(hjust = 0.5),
          legend.position="none",
          legend.text=element_text(size=6),
          legend.title=element_text(size=8),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )
    }
}

.env$gg.manhattan_chr <- function(df1, threshold,  ylims, title,hlight,col){
    #https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/
    #GGplot2-Manhattan-Plot-Function fixed for FST
    #Needed packages
        # library(readr)
        # library(ggrepel)
        # library(ggplot2)
        # library(dplyr)
        # library(RColorBrewer)
      # format df
    if(is.data.frame(df1)==FALSE){
        print("Input must be a data.frame")
    }
    else if(is.data.table(df1)==TRUE){
        print("Input must be data.frame, not data.table, use as.data.frame()")
    }else{      
        if((missing(threshold))){
            threshold<-5
        }
        if((missing(col))){
            col<-c("#5D82BB", "#3B64A5", "#1E4F9E", "#103B7E", "#082B64")
        }
        if((missing(title))){
            title<-"Unknown"
        } 

        df<- subset(df1,Z_FST>0)
        df.tmp <- df %>%

        # Compute chromosome size
        group_by(CHR) %>% 
        dplyr::summarise(chr_len=max(POS)) %>%

        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
        dplyr::select(-chr_len) %>%

        # Add this info to the initial dataset
        left_join(df, ., by=c("CHR"="CHR")) %>%

        # Add a cumulative position of each SNP
        arrange(CHR, POS) %>%
        mutate( BPcum=POS+tot) %>%

        # Add highlight and annotation information
        # mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
        mutate( is_annotate=ifelse(Z_FST > 5, "yes", "no")) %>%
        mutate( is_annotate1=ifelse(Z_FST > 5.5, "yes", "no"))

        # get chromosome center positions for x-axis
        axisdf <- df.tmp %>% group_by(CHR) %>% summarize( #Add POS here if you want to look single CHR
        center=( max(BPcum) + min(BPcum) ) / 2 )

        ggplot(df.tmp, aes(x=BPcum, y=Z_FST)) +
        # Show all points
        geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1) +
        scale_color_manual(values = rep(col, 38 )) +

        # custom X axis: 
        scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) + #Remove for single CHR (i think)
        scale_y_continuous(expand = c(0, 0), limits = ylims) +  #expand=c(0,0)removes space between plot area and x axis 

        # add plot and axis titles
        ggtitle(title) +
        labs(x = "Chromosome",y=expression(ZF["ST"])) +

        # add genome-wide sig and sugg lines
        geom_hline(yintercept = 5,color="orange",linetype="dashed")+
        geom_hline(yintercept = 5.5,color="orange")+

        # Add highlighted points
        geom_point(data=subset(df.tmp, is_annotate=="yes"), color="orange", size=1) +
        geom_point(data=subset(df.tmp, is_annotate1=="yes"), color="orangered", size=2) +


        # Add label using ggrepel to avoid overlapping, I usally don't use the labels so I've just commented it out
        #geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +


        # Custom the theme:
        theme_bw(base_size = 22) +
        theme( 
          plot.title = element_text(hjust = 0.5),
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )
    }
}

.env$gg.heatmap<-function(df,title){



    df %>% gather(key="Group",value="MAF","FCR","Eriks") %>%  
    ggplot(aes(y=Group,x=as.factor(POS),fill=MAF))+geom_tile()+
    
    # geom_rect(data=df, size=1, fill=NA, colour="black",
    #     aes(xmin=hlight - 0.5, xmax=hlight + 0.5, ymin=hlight - 0.5, ymax=hlight + 0.5)) +
    
    labs(x = "Position")+



    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0))+
    ggtitle(title) +

   
    # Custom theme:
    theme_bw(base_size = 22) +
    #theme_classic() +
    theme(
        #legend.position="none",
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),        
        axis.text.x = element_text(
            angle = 90
            #,hjust = -40
            ,margin = NULL
            ,size=8 
            ,vjust = 0.5 #styrer højre-venstre når teksten er vendt
            ),
        axis.ticks.x=element_blank(),
    )
}



#Add functions to environment. When doing it like this you get alle the functions directly loaded, but without them counting as objects in your normal environment, to see objects from this profile write objects(.env). If it says that one of the functions is unknown, write attach(.env) in commandline. This typically happens if there's keyboard input while it's loading. Or if there's a messed up parentheses somewhere.
attach(.env)

#Message when shutting R down.
.Last<-function(){
    cat("Goodbye \n")
}

