################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
###                                     ########################################################################################################
### File_S1_-_R_script_Nätt_et_al_2019  ########################################################################################################
###                                     ########################################################################################################
################################################################################################################################################
##--------------------------------------------------------------------------------##                                                          ##
## Introduction:                                                                                                                              ##  
##--------------------------------------------------------------------------------##                                                          ##
## This R script was supplemented with the paper Nätt et al 2019 (PLoS Biology).                                                              ##  
## Following the instructions below will generate some of the results published in the article.                                               ##
## Script usage comes with no warranties. The most up to date version on the script can be found on github (https://github.com/Danis102).     ## 
## If encountering any problems please post a comment on the github repository or email Dr Daniel Nätt (daniel.natt@liu.se).                  ##
##                                                                                                                                            ##
## Setup:                                                                                                                                     ##
## The script uses data from the Supplemantery Tables file. More particularly:                                                                ##
##    Table S2 (Phenotype data)                                                                                                               ##
##    Data S1 (Sports output)                                                                                                                 ##
##    Data S2 (Mintmap output)                                                                                                                ##
##    Table S4 (Mintmap results)                                                                                                              ##
##                                                                                                                                            ##
## Note: It is important that the original Supplemantery Tables file, provided as an Excel (.xlsx), is unpacked, and that the correct path    ##    
## to this file is provided to the script.                                                                                                    ##  
##                                                                                                                                            ##
## Copyright (R)                                                                                                                              ##
## The script can be used under the Creative Commons Attribution (CC BY) license. This means tha any one may copy, distribute, or reuse       ##
## the content as long as the author and original source are properly cited.                                                                  ##
##                                                                                                                                            ##
## Citation:                                                                                                                                  ##
## Nätt et al (2019) "Human Sperm Displays Rapid Responses to Diet", Plos Biology                                                             ##
##                                                                                                                                            ##
##                                                                                                                                            ##
##                                                                                                                                            ##
##                                                                                                                                            ##
################################################################################################################################################
##--------------------------------------------------------------------------------##
####       Step 1: Install dependencies 
##--------------------------------------------------------------------------------##         

## Note: 
## This script was developed on a Linux (Mint 19.1 Tessa) system.
## using R version 3.4.4.
## OBS! It has been tested on Windows 10 but not on Mac/iOS.

## The following r-packages are needed:
# System wide: "ggplot2", "readxl", "reshape2", "RColorBrewer", "scales", "ggthemes", "cowplot", "dendextend", "ape"
# Linux specific: "extrafont"

## Check dependancies 
 # Change this to FALSE if you are using windows
check.depend <- function(){
                      dep <- c("ggplot2", "readxl", "reshape2", "RColorBrewer", "scales", "ggthemes",  "cowplot", "dendextend", "ape")
                      inst <- dep %in% rownames(installed.packages())
                      if(any(inst==FALSE)){warning("The following package was not installed: ", paste(dep[!inst], collapse ="; "))}
                      else{cat("Good! All packages appear to be installed.")}
                      }
check.depend()



###########################################################################################################################################
##--------------------------------------------------------------------------------##
####       Step 2: Download and unpack the input file
##--------------------------------------------------------------------------------## 

# The Supplemenatary tables file published with the article is used as input.
# The file is availble at github (https://github.com/Danis102) or on the journal homepage. 
# The file names must be: 
#   S2_Table.xlsx (phenotype data)
#   S1_Data.xlsx (sports output)
#   S2_Data.xlsx (mintmap output)
#   S4_Table.xlsx (tsRNA annotation etc)

###########################################################################################################################################
##--------------------------------------------------------------------------------##
####       Step 3: Add the to location if the input file 
##--------------------------------------------------------------------------------## 

### Add the path were you have stored the supplementary files
path="<Your path to the folder with unpacked supplementary files>"


###########################################################################################################################################
##--------------------------------------------------------------------------------##
####       Step 4: Generate and run functions.
##--------------------------------------------------------------------------------## 

## This function will generate all data and graphs in Fig 2 and Fig 3: 
generate.results <- function(path){
                            #### Set up envirnoment ###
                            options(scipen=999)
                            require(ggplot2)
                            require(readxl)
                            require(reshape2)
                            require(RColorBrewer)
                            require(scales)
                            require(ggthemes)
                            require(cowplot)
                            require(dendextend)
                            require(ape)
                            #require(extrafont)
                            ###############################################################################################################
                            #### Set up colors colors ###
                                      colfunc_sports <- colorRampPalette(c("#094A6B", "#FFFFCC", "#9D0014"))
                                      rgb_vec_sports <- c(colfunc_sports(6), "#6E6E6E", "#BCBCBD")
                                      rgb_vec_sports <- c(colfunc_sports(6), "#6E6E6E")
                                      rgb_vec_sports[3] <- "#A5BF7E"
                                      #show_col()
                                      
                                      colfunc_mint_types <- colorRampPalette(c("#094A6B", "#FFFFFF", "#9D0014"))
                                      rgb_vec_mint_types <- rev(colfunc_mint_types(5))
                                      
                                      colfunc_mint_isodecoders <- colorRampPalette(c("white", "black"))
  
                            #### Setup functions ####
                                ##-------- Function for reading data 
                                  read.data <- function(path){
                                                    Sports_PAC <- list(Pheno=NA, Anno=NA, countTable=NA)
                                                    Mint_PAC <- list(Pheno=NA, Anno=NA, countTable=NA, Result=NA)
                                                   
                                                  ## Sports S1 Data ##
                                                          Sports <- as.data.frame(read_xlsx(paste0(path, "/S1_Data.xlsx"), sheet=1, skip=4, col_names = TRUE))
                                                          Sports_PAC$Anno <- Sports[,1:5]
                                                          Sports_PAC$countTable <- Sports[,-(1:5)]
                                                          rownames(Sports_PAC$Anno) <-  Sports_PAC$Anno$"RNA_sequence"
                                                          rownames(Sports_PAC$countTable) <- Sports_PAC$Anno$"RNA_sequence"

                                                    ## Mintmap S2 Data ## 
                                                          Mintmap <- as.data.frame(read_xlsx(paste0(path, "/S2_Data.xlsx"), sheet=1, skip=4, col_names = TRUE))
                                                          Mint_PAC$Anno <- Mintmap[, c(1:7, 53)] 
                                                          Mint_PAC$countTable <- Mintmap[,-c(1:7,53)]
                                                          rownames(Mint_PAC$Anno) <-  Mint_PAC$Anno$"tsRNA sequence"
                                                          rownames(Mint_PAC$countTable) <-  Mint_PAC$Anno$"tsRNA sequence"
                                                          colnames(Mint_PAC$Anno)[3] <-  "Length"
                                                    
                                                    ## Results statistics Minmap S4 Table ##
                                                          Mint_res <- suppressMessages(as.data.frame(read_xlsx(paste0(path, "/S4_Table.xlsx"), sheet=1, skip=3, col_names = TRUE)))
                                                          Mint_PAC$Result <- Mint_res[1:1747, c(-7, -12, -16, -17, -19)] 
        
                                                    ## Pheno Table S2 ##
                                                          Pheno_1 <- suppressMessages(as.data.frame(read_xlsx(paste0(path, "/S2_Table.xlsx"), sheet=1, skip=3, col_names = TRUE)))
                                                          Pheno_2 <- Pheno_1[-c(8:10, 19:21, nrow(Pheno_1)), c(1, 19:ncol(Pheno_1))]         
                                                          Pheno_3 <- t(Pheno_2[,-1])
                                                          colnames(Pheno_3) <- Pheno_2[,1]
                                                          Pheno_4 <- cbind(data.frame(Participant = do.call("rbind", strsplit(rownames(Pheno_3), "_"))[,1], Diet= do.call("rbind", strsplit(rownames(Pheno_3), "_"))[,2]), Pheno_3)
                                                          Sports_PAC$Pheno <- Pheno_4
                                                          Mint_PAC$Pheno <- Pheno_4[!grepl("_Start", rownames(Pheno_3)),]
                                                          
                                                    ## Match order ##
                                                          Mint_PAC$countTable <- Mint_PAC$countTable[,match(rownames(Mint_PAC$Pheno), colnames(Mint_PAC$countTable))]
                                                          Sports_PAC$countTable <- Sports_PAC$countTable[,match(rownames(Sports_PAC$Pheno), colnames(Sports_PAC$countTable))]
                                                          Mint_PAC$countTable <- Mint_PAC$countTable[match(rownames(Mint_PAC$Anno), rownames(Mint_PAC$countTable)),]
                                                          Sports_PAC$countTable <- Sports_PAC$countTable[match(rownames(Sports_PAC$Anno), rownames(Sports_PAC$countTable)),]
                                                          
                                                          stopifnot(identical(rownames(Mint_PAC$Pheno), colnames(Mint_PAC$countTable)))
                                                          stopifnot(identical(rownames(Sports_PAC$Pheno), colnames(Sports_PAC$countTable)))
                                                          stopifnot(identical(rownames(Mint_PAC$Anno), rownames(Mint_PAC$countTable)))
                                                          stopifnot(identical(rownames(Sports_PAC$Anno), rownames(Sports_PAC$countTable)))
                                                          
                                                    return(list(Sports=Sports_PAC, Mint=Mint_PAC))
                                }
                                ##-------- Function for generating RPM (Sports)
                                  generate.RPM <- function(PAC){
                                                      lib_sizes <- colSums(PAC$countTable)
                                                      counts_RPM <- data.frame(matrix(NA, nrow=nrow(PAC$countTable), ncol=ncol(PAC$countTable)))
                                                      colnames(counts_RPM) <- colnames(PAC$countTable)
                                                      rownames(counts_RPM) <- rownames(PAC$countTable)
                                                      for (i in 1:length(lib_sizes)){
                                                                        counts_RPM[,i] <- (PAC$countTable[,i]/(lib_sizes[i]/1000000))
                                                                        }
                                                      return(counts_RPM)
                                }
                                ##-------- Function for filtering based on RPM and size (both Mint and Sports)
                                  RPM.filt <- function(PAC, size, rpm, pheno_targets){
                                                          ### Check col and row names
                                                          if(!identical(rownames(PAC$Anno), rownames(PAC$RPM))){stop("Error: Not matching rownames in input files! (Anno vs RPM)")}
                                                          if(!identical(rownames(PAC$Pheno), colnames(PAC$RPM))){stop("Error: Not matching rownames in input files!\nRun the order.res function.")}
                                                          ### Subset data by groups
                                                              sub <- as.character((PAC$Pheno[pheno_targets[[1]]])[,1]) %in% pheno_targets[[2]]
                                                              PAC$RPM  <- PAC$RPM[,sub] 
                                                              PAC$Pheno  <- PAC$Pheno[sub,]
                                                          ### RPM filter
                                                              start <- nrow(PAC$Anno)
                                                              idx_rpm <- rowSums(PAC$RPM >= rpm[1]) >= round(ncol(PAC$RPM)*(rpm[2]*0.01))
                                                          ### Size filter 
                                                              idx_size <- PAC$Anno$Length >= size[1] & PAC$Anno$Length <= size[2]
                                                          ### In total
                                                              idx <- rowSums(cbind(idx_rpm, idx_size)) == 2
                                                          ### Apply filter
                                                              PAC$Anno  <- PAC$Anno[idx,]
                                                              PAC$countTable  <- PAC$RPM[idx,]
                                                              PAC$RPM  <- PAC$RPM[idx,]
                                                              end <- nrow(PAC$Anno)
                                                              #cat("\n\nIn total:\n", end, "was retained.\n", start-end, "was removed.")
                                                          ## Check
                                                          if(!identical(rownames(PAC$Anno), rownames(PAC$countTable))){stop("Error: Not matching rownames in input files! (Anno vs countTable)")}
                                                          if(!identical(rownames(PAC$Pheno), colnames(PAC$countTable))){stop("Error: Not matching rownames/colnames in input files! (Pheno vs countTable)")}
                                                          if(!identical(rownames(PAC$Anno), rownames(PAC$RPM))){stop("Error: Not matching rownames in input files! (Anno vs RPM)")}
                                                          if(!identical(rownames(PAC$Pheno), colnames(PAC$RPM))){stop("Error: Not matching rownames/colnames in input files! (Pheno vs RPM)")}
                                                          #cat("\nGood! Row and col names identical after filtering. Ready to proceed.") 
                                                          return(PAC)
                                  }
                                ##-------- Function for calculation for data (Sports)
                                  FC.bio.sports <- function(PAC_RPM, group_column, groups, paired=TRUE, noAnno.rm=TRUE){
                                        ### log2FC Paired
                                        if(paired==TRUE){
                                                    PAC_RPM$mean_Log2FC <- data.frame(Biotype=PAC_RPM$Anno$Biotype, 
                                                                                        mean_FC = rowMeans(log2((PAC_RPM$RPM[,grepl(groups[1], PAC_RPM$Pheno[,group_column])]) / (PAC_RPM$RPM[,grepl(groups[2], PAC_RPM$Pheno[,group_column])]))))    
                                        ### log2FC Independent
                                        }else{
                                                    PAC_RPM$mean_Log2FC <- data.frame(Biotype=PAC_RPM$Anno$Biotype, 
                                                                                        mean_FC = log2(rowMeans(PAC_RPM$RPM[,grepl(groups[1], PAC_RPM$Pheno[,group_column])]) / rowMeans(PAC_RPM$RPM[,grepl(groups[2], PAC_RPM$Pheno[,group_column])])))    
                                        }
                                        ### Sample sum biotype calculation
                                        PAC_RPM$Sample_sumRPM_biotype <- aggregate(PAC_RPM$RPM, list(as.character(PAC_RPM$Anno$Biotype)), sum)
                                        colnames(PAC_RPM$Sample_sumRPM_biotype)[1] <-"Biotype" 
                                        PAC_RPM$Experiment_sumRPM_biotype <- aggregate(rowMeans(PAC_RPM$RPM), list(as.character(PAC_RPM$Anno$Biotype)), sum)
                                        colnames(PAC_RPM$Experiment_sumRPM_biotype) <- c("Biotype", "meanRPM") 
                                        return(PAC_RPM)
                                        }
                                ##-------- Function for generating jitterplot (Sports)
                                  plot.jitter <- function(input, FC_col, feat_col, limits, Ypos_n, colors="Black"){
                                        					require(ggplot2)
                                                  if(class(input)=="data.frame"){
                                        						input <- list(input=input)}
                                        					plot_lst <- list(NA)
                                        					for(i in 1:length(input)){
                                        							exp <- names(input)[i]
                                        							input[[i]] <- input[[i]][!is.na(input[[i]][, colnames(input[[i]])==FC_col]),]
                                        							if(i==length(input)){
                                        							  perc_up_agg <- aggregate(input[[i]][,FC_col], list(as.character(input[[i]][,feat_col])), function(x){ sum(as.numeric(x > 0))/length(x)})
                                        							  perc_up <- round(perc_up_agg$x, digits=3)*100
                                        							  perc_up <-  perc_up[match(as.character(levels(input[[i]][,feat_col])), as.character(perc_up_agg$Group.1))]
                                        								plot_lst[[i]] <- ggplot(input[[i]], aes_string(x=feat_col, y=FC_col, col=feat_col, fill=feat_col))+
                                        								        geom_hline(yintercept=0, col="#707177", cex=0.6) +
                                        								        geom_jitter(position=position_jitter(0.2), cex=1.5)+
                                        	                      stat_summary(geom = "crossbar", fun.y=median, fun.ymax = median, fun.ymin = median, width=0.7, cex=0.4, position = "identity", col="Black") +								  			
                                        								        geom_boxplot(width=0.3, fill="white", col="black", alpha=0.7,  outlier.shape = NA)+
                                        												geom_text(stat="count", aes(label=paste0("n=",..count.., "\nup:", perc_up, "%")), size=5, y=Ypos_n, col="Black") +
                                        												labs(title=paste0(exp) , x="Biotype" , y = paste0(FC_col)) +
                                        												theme_classic()+
                                        												scale_y_continuous(limits =limits) +
                                        												theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 0.95), axis.text.y = element_text(size=15), axis.title.y= element_blank())+
                                        												scale_color_manual(values=colors)
                                        								        #coord_flip()
                                        								names(plot_lst)[i] <- exp
                                        								} else {
                                        								plot_lst[[i]] <- ggplot(input[[i]], aes_string(x=feat_col, y=FC_col, fill=feat_col))+
                                        												geom_violin (width=0.9, trim=FALSE, scale="width")+
                                        												geom_boxplot(width=0.2, fill="white", outlier.shape = NA)+
                                        												geom_text(stat="count", aes(label=paste0("n=",..count..)), size=5, y=Ypos_n) +
                                        												labs(title=paste0(exp) , x="Biotype" , y = paste0(FC_col)) +
                                        												geom_hline(yintercept=0)+
                                        												theme_classic()+
                                        												scale_y_continuous(limits = limits) +
                                        												theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 0.95), axis.text.y = element_text(size=15), axis.title.y= element_blank())+
                                        												coord_flip()
                                        								names(plot_lst)[i] <- exp
                                        								}
                                        					}
                                        				return(plot_lst)
                                }
                                ##-------- Function for generating stacked bargraph (Sports)
                                  plot.stacked.bar <- function(input, target.biotypes, color_vec){
                                                data <- input[input$Biotype %in% as.character(target.biotypes), ]
                                                data_perc <- apply(data[,-1], 1, function(x){x/colSums(data[,-1])})
                                                colnames(data_perc) <- data[,1]
                                                data_perc_long <- melt(data_perc)
                                                colnames(data_perc_long) <- c("Sample", "Biotype", "Percent")
                                                data_perc_long$Biotype <- factor(as.character(data_perc_long$Biotype), levels=rev(target.biotypes))
                                                data_perc_long <- cbind(data_perc_long, data.frame(Participant=do.call("rbind", strsplit(as.character(data_perc_long$Sample), split="_"))[,1],
                                                                                     Diet=do.call("rbind", strsplit(as.character(data_perc_long$Sample), split="_"))[,2]))
                                                data_perc_long <- data_perc_long[order(as.numeric(gsub("S", "", data_perc_long$Participant))),]
                                                data_perc_long$Sample <- factor(as.character(data_perc_long$Sample), levels=as.character(unique(data_perc_long$Sample))) 
                                                p1<- ggplot(data_perc_long, aes(x=Sample, y=Percent, fill=Biotype)) +
                                                                              geom_bar(stat="identity", width=1.0, size=1.0) + 
                                                								              geom_hline(yintercept=0, col="#707177", cex=0.6) +
                                                                              geom_rangeframe(aes(y=c(0, rep(1, length(Percent)-1))))+
                                                  												    scale_fill_manual(values=rev(color_vec))+
                                                                              theme_tufte()+
                                                                              theme(
                                                                  					      plot.caption =  element_text(size=12, face= "bold"),
                                                                  					      axis.title.y=element_text(size=16, face= "bold"), 
                                                                  					      axis.title.x= element_blank(), 
                                                                  					      axis.text=element_text(size=14),
                                                                  					      axis.text.x = element_text(angle=45, hjust=1),
                                                                  					      panel.background = element_blank()) 
                                                return(p1)
                                                }
                                ##-------- Function for generating pie plot biotype (Sports)
                                  plot.pie <- function(data, target.biotypes, color_vec, angle=-25){
                                                    data <- data[match(rev(target.biotypes), as.character(data$Biotype)), ]
                                                    data$percent <- (data$meanRPM/sum(data$meanRPM))*100 
                                                    p1  <- suppressMessages(pie(data$meanRPM, labels=paste(data$Biotype, round(data$percent , digits=0), "%"), col=rev(color_vec), init.angle = angle))
                                                    print(p1)
                                                    rp <- recordPlot()
                                                    return(p1)
                                  }
                                ##-------- Function for generating pie plot percent up or down (Sports)
                                  plot.perc.pie <- function(input,  target.biotypes, colors, angle){
                                                                input$Biotype <- factor(input$Biotype, levels=target.biotypes)
                                                							  perc_up_agg <- aggregate(input$mean_FC, list(input$Biotype), function(x){ sum(as.numeric(x > 0))/length(x)})
                                                							  perc_up <- round(perc_up_agg$x, digits=4)*100
                                                							  df_1 <- data.frame(biotype=perc_up_agg[,1], up=perc_up, down=100-perc_up)
                                                                r_lst <- setNames(split(df_1, seq(nrow(df_1))), df_1$biotype)
                                                                pie_lst <- lapply(r_lst, function(x){
                                                                            df <- melt(x)
                                                                            p <- ggplot(df, aes(x="", y=value, fill=variable)) +
                                                                                    geom_bar(width = 1, stat = "identity") +
                                                                                    coord_polar("y", start=angle) +
                                                                                    scale_fill_manual(values=colors) +
                                                                                    ggtitle(paste0("%up: ", x$up, "\n%down: ", x$down)) +
                                                                                    #geom_text(aes(y=-1, poistion="dodge", label = paste0(round(value, digits=0),"%")), size=6, color=c("white","black"))+
                                                                                    theme_minimal()+
                                                                                    theme(
                                                                                    axis.text=element_blank(),
                                                                                    legend.position="none",
                                                                                    axis.title = element_blank(),
                                                                                    panel.border = element_blank(),
                                                                                    panel.background = element_blank(),
                                                                                    panel.grid=element_blank(),
                                                                                    axis.ticks = element_blank(),
                                                                                    plot.title=element_text(size=12, hjust = 0.5))
                                                                            return(p)
                                                                            })
                                                               p <- plot_grid(pie_lst[[1]],
                                                                                pie_lst[[2]],
                                                                                pie_lst[[3]],
                                                                                pie_lst[[4]],
                                                                                pie_lst[[5]],
                                                                                pie_lst[[6]],
                                                                                pie_lst[[7]],
                                                                                labels = names(pie_lst), hjust = 0, vjust = 1, ncol = 7, nrow = 1,  rel_widths = c(1, 1, 1), rel_heights = c(1, 1, 1))
         
                                                                
                                                          return(p)
                                  }
                                ##-------- Function for calculations for tsRNA types (Mint)  
                                  prep.mint.types <- function(PAC, paired=FALSE, pheno_target, anno_target){
                                                                              prep_lst <- list(NULL)
                                                                              data_types <- aggregate(PAC$RPM, list(factor(PAC$Anno[,anno_target])), sum, simplify=TRUE)
                                                                              if(paired==TRUE){
                                                                                        data_types_mat <- data_types[,-1]
                                                                                        prep_lst$Info$Description <-"
                                                    Experimental design: Paired samples.
                                                    Fold changes were calculated on the ratios between repeated measures of the same individual.
                                                    'Anno_aggregated' = Sums were calculated over the factor in anno_targets.
                                                    'Types' = Sums were calculated over anno_targets before ratios were calculated.
                                                    'Frag' =  No aggregation, ratios were calculated on individual fragments."
                                                                                        prep_lst$Info$Paired <- TRUE
                                                                                        prep_lst$Info$Pheno_target <- pheno_target
                                                                                        prep_lst$Info$Anno_target <- anno_target
                                                                                        prep_lst$Pheno <- PAC$Pheno[PAC$Pheno[,pheno_target[[1]]] %in% pheno_target[[2]],]
                                                                                        prep_lst$Anno_aggregated <- melt(data_types, id.vars="Group.1")
                                                                                        paired_types_log2FC <- cbind(data.frame(Group.1=data_types[,1]),  log2(data_types_mat[, PAC$Pheno[pheno_target[[1]]] == pheno_target[[2]][1]] / data_types_mat[,PAC$Pheno[pheno_target[[1]]] == pheno_target[[2]][2]]))
                                                                                        prep_lst$Types_log2FC <- melt(paired_types_log2FC, id.vars="Group.1")
                                                                                        paired_frag_log2FC <- cbind(data.frame(Group.1=PAC$Anno[,anno_target]), log2(PAC$RPM[,PAC$Pheno[pheno_target[[1]]] == pheno_target[[2]][1]] / PAC$RPM[,PAC$Pheno[pheno_target[[1]]] == pheno_target[[2]][2]]))
                                                                                        prep_lst$Frag_log2FC <- melt(paired_frag_log2FC, id.vars="Group.1")
                                                                                        prep_lst$Frag_log2FC$MINTbase_ID <- rep(PAC$Anno$MINTbase_ID, times=nrow(prep_lst$Pheno)/2)
                                                                                        }
                                                                              if(paired==FALSE){
                                                                                        prep_lst$Info$Description <-"
                                                    Experimental design: Independent samples.
                                                    Fold change calculations were based on ratios between independent groups.
                                                    'Anno_aggregated' = Sums were calculated over the factor in anno_targets.
                                                    'Types' = Sums of anno_targets were calculated before ratios between groups were calculated.
                                                    'Frag' = No aggregation, ratios between groups were calculated on individual fragments."
                                                                                        prep_lst$Info$Paired <- FALSE
                                                                                        prep_lst$Info$Pheno_target <- pheno_target
                                                                                        prep_lst$Info$Anno_target <- anno_target
                                                                                        prep_lst$Pheno <- PAC$Pheno[PAC$Pheno[,pheno_target[[1]]] %in% pheno_target[[2]],]
                                                                                        prep_lst$Anno_aggregated <- melt(data_types, id.vars="Group.1")
                                                                                        data_types_mat <- data_types[,-1]
                                                                                        prep_lst$Types_log2FC <- data.frame(Group.1=data_types[,1],  value= log2(rowMeans(data_types_mat[,PAC$Pheno[pheno_target[[1]]] == pheno_target[[2]][1]]) / rowMeans(data_types_mat[,PAC$Pheno[pheno_target[[1]]] == pheno_target[[2]][2]])))
                                                                                        prep_lst$Frag_log2FC <- data.frame(Group.1=PAC$Anno[,anno_target],  value= log2(rowMeans(PAC$RPM[,PAC$Pheno[pheno_target[[1]]] == pheno_target[[2]][1]]) / rowMeans(PAC$RPM[,PAC$Pheno[pheno_target[[1]]] == pheno_target[[2]][2]])))
                                                                                        prep_lst$Frag_log2FC$MINTbase_ID <- PAC$Anno$MINTbase_ID
                                                                                        }
                                                                              return(prep_lst[-1])
                                  }
                                ##-------- Function for plotting tsRNA types (Mint)                                   
                                  plots.mint.types <- function(prep, color_vec){
                                                                    plot_lst <- list(NULL)
                                                                    lvls <-   c("5'-tRF", "5'-half", "i-tRF", "3'-tRF", "3'-half")
                                                                    ## Bars for mean RPM per type - log10
                                                                    ymax <- max(prep$Anno_aggregated$value)*0.6
                                                                    prep$Anno_aggregated$Group.1 <- factor(prep$Anno_aggregated$Group.1, levels=lvls)
                                                                    plot_lst$Types_MeanRPM_log10 <- ggplot(prep$Anno_aggregated, aes(x=Group.1, y=value, fill=Group.1)) +
                                                                                                 				stat_summary(geom = "errorbar",  width=0.3, size=1.0, fun.data = mean_se, position = "identity", na.rm =TRUE) +
                                                                                              					stat_summary(geom = "bar", colour="black", width=0.6, size=1.0, fun.y = mean, na.rm =TRUE) +
                                                                              												  geom_text(stat="count", aes(label=paste0("n=",..count..)), size=5, y=ymax, col="Black", na.rm =TRUE) +
                                                                                            					  labs(title=paste0(prep$Info$Pheno_target[[2]], collapse=" vs ")) + 
                                                                                               					ylab("Mean log10 RPM +/- SE") + 
                                                                                              					scale_fill_manual(values=color_vec) +
                                                                                              					scale_x_discrete(labels=levels(prep$Anno_aggregated$Group.1)) +
                                                                                                        scale_y_continuous(trans = 'log10', breaks=c(1,10,100,1000,10000,100000)) +
                                                                                                        geom_rangeframe(aes(x=Group.1, y=c(1, rep(100000, length(value)-1))))+   
                                                                                                        theme_tufte()+
                                                                                                        theme(legend.position="none",
                                                                                              					      axis.title.y=element_text(size=16, face= "bold"), 
                                                                                              					      axis.title.x= element_blank(), 
                                                                                              					      axis.text=element_text(size=14),
                                                                                              					      axis.text.x = element_text(angle=45, hjust=1),
                                                                                              					      panel.background = element_blank())
                                                                    #Paired#############################################################################################################
                                                                    if(prep$Info$Paired==TRUE){
                                                                            ## Error bars for log2_FC types - paired (with jitter)
                                                                            prep$Types_log2FC$Group.1 <- factor(prep$Types_log2FC$Group.1, levels=lvls)
                                                                            ylims = c(-1.5, 1.5)
                                                                            plot_lst$Types_log2FC <- ggplot(prep$Types_log2FC, aes(x= Group.1, y=value, fill= Group.1)) +
                                                                                              					geom_hline(yintercept = 0, linetype="dotted", size=1, color="azure4")+
                                                                                                        geom_jitter(aes(color= Group.1), position=position_jitter(0.1), cex=2.5) +
                                                                                                 					    stat_summary(geom = "errorbar", width=0.5, size=1, fun.data = mean_se, position = "identity") +
                                                                                              					      stat_summary(geom = "point", stroke=1, shape=21, size = 5.0, fun.y = mean, position = "identity") +
                                                                              												  geom_text(stat="count", aes(label=paste0("n=",..count..)), size=5, y=max(ylims), col="Black") +
                                                                                            					  labs(title=paste0(prep$Info$Pheno_target[[2]], collapse=" vs ")) + 
                                                                                               					ylab("Mean log2FC (RPM) +/- SE") +
                                                                                              					scale_fill_manual(values=color_vec) +
                                                                                                        scale_color_manual(values=c(color_vec[1:2], "grey",  color_vec[4:5])) +
                                                                                              					scale_x_discrete(labels=levels(prep$Types_log2FC$Group.1)) +
                                                                                                        coord_cartesian(ylim = ylims) +
                                                                                                        geom_rangeframe(aes(x=Group.1, y=c(ylims[1], rep(ylims[2], length(value)-1))))+
                                                                                                        theme_tufte()+
                                                                                              					theme(legend.position="none",
                                                                                              					      axis.title.y=element_text(size=16, face= "bold"),
                                                                                              					      axis.title.x= element_blank(),
                                                                                              					      axis.text=element_text(size=14),
                                                                                              					      axis.text.x = element_text(angle=45, hjust=1))
                                                                          
                                                                              }
                                                                    #Independent#############################################################################################################
                                                                    if(prep$Info$Paired==FALSE){
                                                                            ## Error bars for log2_FC types - independent
                                                                            prep$Types_log2FC$Group.1 <- factor(prep$Types_log2FC$Group.1, levels=lvls)
                                                                            ylims = c(-0.5, 0.5)
                                                                            plot_lst$Types_log2FC <- ggplot(prep$Types_log2FC, aes(x=Group.1, y=value, fill=Group.1)) +
                                                                                              					geom_hline(yintercept = 0, linetype="dashed", size=1, color="azure4")+
                                                                                              					stat_summary(geom = "point", stroke=0.5, shape=21, size = 5.0, fun.y = mean, position = "identity") +
                                                                              												  geom_text(stat="count", aes(label=paste0("n=",..count..)), size=5, y=max(ylims), col="Black") +
                                                                                            					  labs(title=paste0(prep$Info$Pheno_target[[2]], collapse=" vs ")) + 
                                                                                               					ylab("Mean Log2FC between groups (RPM) +/- SE") + 
                                                                                              					scale_fill_manual(values=color_vec) +
                                                                                              					scale_x_discrete(labels=levels(prep$Types_log2FC$Group.1)) +
                                                                                                        coord_cartesian(ylim = ylims) +
                                                                                                        geom_rangeframe()+                				  
                                                                                                        theme_tufte()+
                                                                                              					theme(legend.position="none",
                                                                                              					      axis.title.y=element_text(size=16, face= "bold"), 
                                                                                              					      axis.title.x= element_blank(), 
                                                                                              					      axis.text=element_text(size=14),
                                                                                              					      axis.text.x = element_text(angle=45, hjust=1))
                                                                            
                                                                            
                                                                            ## Error bars for differencs in RPM - independent
                                                                            ymax <- max(prep$Anno_aggregated$value)*0.7
                                                                            prep$Anno_aggregated$Group.1 <- factor(prep$Anno_aggregated$Group.1, levels=lvls)
                                                                            prep$Anno_aggregated$Group.2 <- rep(prep$Pheno[,prep$Info$Pheno_target[[1]]], each=5)
                                                                            plot_lst$Types_Diff_Errorbar <- ggplot(prep$Anno_aggregated, aes(x=Group.1, y=value, group=interaction(Group.2, Group.1), fill=Group.1)) +
                                                                                                        stat_summary(geom = "errorbar",  width=0.8, size=0.5, fun.data = mean_se, position = "dodge") +
                                                                                              					stat_summary(geom = "point", stroke=0.5, shape=21, size = 5.0, fun.y = mean, position = position_dodge(width=0.8)) +
                                                                              	                        stat_summary(geom = "text", aes(label=paste0("n=",..y..)), size = 4.0, fun.y = length, vjust = -23,  position = position_dodge(width=1)) +
                                                                                              					labs(title=paste0(prep$Info$Pheno_target[[2]], collapse=" vs ")) + 
                                                                                               					ylab("Mean RPM +/- SE") + 
                                                                                              					scale_fill_manual(values=color_vec) +
                                                                                              					scale_x_discrete(labels=levels(prep$Types_log2FC$Group.1)) +
                                                                                                        coord_cartesian(ylim = c(0, ymax)) +
                                                                                                        #scale_y_continuous(limits=c(0, ymax)) +
                                                                                                        theme_tufte()+
                                                                                              					theme(legend.position="none",
                                                                                              					      axis.title.y=element_text(size=16, face= "bold"), 
                                                                                              					      axis.title.x= element_blank(), 
                                                                                              					      axis.text=element_text(size=14),
                                                                                              					      axis.text.x = element_text(angle=45, hjust=1))+
                                                                                                        geom_rangeframe()
                                                                            
                                                                            ## Error bars for differencs in RPM log10 - independent
                                                                            ymax <- max(prep$Anno_aggregated$value)*0.7
                                                                            prep$Anno_aggregated$Group.1 <- factor(prep$Anno_aggregated$Group.1, levels=lvls)
                                                                            prep$Anno_aggregated$Group.2 <- rep(prep$Pheno[,prep$Info$Pheno_target[[1]]], each=5)
                                                                            ns <- as.numeric(table(prep$Pheno[,prep$Info$Pheno_target[[1]]]))
                                                                            plot_lst$Types_Diff_Errorbar <- ggplot(prep$Anno_aggregated, aes(x=Group.1, y=value, group=interaction(Group.2, Group.1), fill=Group.1)) +
                                                                                                        stat_summary(geom = "errorbar",  width=0.8, size=0.5, fun.data = mean_se, position = "dodge") +
                                                                                              					stat_summary(geom = "point", stroke=0.5, shape=21, size = 5.0, fun.y = mean, position = position_dodge(width=0.8)) +
                                                                                              					labs(title=paste0(prep$Info$Pheno_target[[2]], collapse=" vs "), caption=paste0("n= ", ns[1], "/", ns[2])) + 
                                                                                               					ylab("Mean RPM log10 +/- SE") + 
                                                                                              					scale_fill_manual(values=color_vec) +
                                                                                              					scale_x_discrete(labels=levels(prep$Types_log2FC$Group.1)) +
                                                                                                        scale_y_continuous(trans = 'log10', breaks=c(1,10,100,1000,10000,100000,1000000)) +
                                                                                                        geom_rangeframe(aes(x=Group.1, y=c(1, rep(1000000, length(value)-1))))+   
                                                                                                        theme_tufte()+
                                                                                                        theme(legend.position="none",
                                                                                              					      plot.caption =  element_text(size=12, face= "bold"),
                                                                                              					      axis.title.y=element_text(size=16, face= "bold"), 
                                                                                              					      axis.title.x= element_blank(), 
                                                                                              					      axis.text=element_text(size=14),
                                                                                              					      axis.text.x = element_text(angle=45, hjust=1),
                                                                                              					      panel.background = element_blank())                          
                                      
                                                                            ## Error bars for log2_FC types - independent
                                                                            plot_lst$Types_Diff_Boxplot <-ggplot(prep$Anno_aggregated, aes(x=Group.1, y=value, group=interaction(Group.2, Group.1), fill=Group.1))+
                                                                                                      	geom_hline(yintercept=0, col="azure4")+
                                                                                                      	geom_boxplot(width=0.7, aes(fill=Group.1, group=interaction(Group.2, Group.1)), outlier.shape = NA)+
                                                                                                      	geom_point(cex=2.0, position = position_jitterdodge(jitter.width=0.2), aes(color=Group.2, group=interaction(Group.2, Group.1)))+
                                                                                                      	expand_limits(y=0.8*max(prep$Anno_aggregated$value)) +
                                                                                                      	xlab(NULL)+
                                                                                                        ylab("RPM") + 
                                                                                                      	scale_fill_manual(values=c("#094A6B", "#3F7280", "#769A96", "#ACC2AB", "#E3EAC1"))+
                                                                                                      	scale_color_manual(values=c("grey", "black"))+
                                                                                                      	scale_shape_discrete(solid=TRUE)+
                                                                                                        theme_tufte()+
                                                                                              					theme(legend.position="none",
                                                                                              					      axis.title.y=element_text(size=16, face= "bold"), 
                                                                                              					      axis.title.x= element_blank(), 
                                                                                              					      axis.text=element_text(size=14),
                                                                                              					      axis.text.x = element_text(angle=45, hjust=1))+
                                                                                                        geom_rangeframe()
                                                                                                 				
                                                                            }
                                                                return(plot_lst[-1])
                                                                  
                                  }
                                ##-------- Function for calculations of isodecoders (Mint)
                                  prep.mint.isodec <- function(PAC, paired=FALSE, filter=100, pheno_targets = list(target_col="Diet", target_groups=c( "Sugar", "Healthy")), col.IDs="Participant") {
                                                  PAC$Pheno[,pheno_targets[[1]]] <- factor(PAC$Pheno[,pheno_targets[[1]]], levels=pheno_targets[[2]])
                                                  tRNA_sub <- as.character(unique(do.call("c", strsplit(PAC$Anno$Isodecoder, ", "))))
                                        					tRNA_sub_sum_df <- data.frame(matrix(NA, nrow=length(tRNA_sub), ncol=ncol(PAC$RPM))) # Make an empty dataframe
                                        					colnames(tRNA_sub_sum_df) <- colnames(PAC$RPM)
                                        					tRNA_sub_perc <- data.frame(matrix(NA, nrow=length(tRNA_sub), ncol=5)) # Make an empty dataframe
                                        					colnames(tRNA_sub_perc) <- c("tsRNA_5'", "half_5'", "tsRNA_i'", "half_3'", "tsRNA_3'")
                                        					# Loop for making ind means over tRNA isotype
                                        					# And for making % tRF type destributions
                                        						for(i in 1:length(tRNA_sub)){
                                        							nam <- tRNA_sub[i]
                                        							sub <- PAC$RPM[grepl(paste0("\\<", tRNA_sub[i], "\\>"), PAC$Anno$Isodecoder),] # search for only whole words
                                        							tRNA_sub_sum_df[i,] <- colSums(sub)
                                        							rownames(tRNA_sub_sum_df)[i] <- nam
                                        							sub_anno <- PAC$Anno[grepl(paste0("\\<", tRNA_sub[i], "\\>"), PAC$Anno$Isodecoder),]
                                        								stopifnot(identical(rownames(sub), rownames(sub_anno))) 
                                        							sub_means <- data.frame(tsRNA_type=sub_anno$"tsRNA type", mean_RPM=rowMeans(sub), row.names = rownames(sub))
                                        							total <- sum(sub_means$mean_RPM)
                                        							tRNA_sub_perc$"tsRNA_5'"[i] <- (sum(sub_means$mean_RPM[sub_means$tsRNA_type == "5'-tRF"])/total)*100
                                        							tRNA_sub_perc$"half_5'"[i] <- (sum(sub_means$mean_RPM[sub_means$tsRNA_type == "5'-half"])/total)*100
                                        							tRNA_sub_perc$"tsRNA_i'"[i] <- (sum(sub_means$mean_RPM[sub_means$tsRNA_type == "i-tRF"])/total)*100
                                        							tRNA_sub_perc$"tsRNA_3'"[i] <- (sum(sub_means$mean_RPM[sub_means$tsRNA_type == "3'-tRF"])/total)*100
                                        							tRNA_sub_perc$"half_3'"[i] <- (sum(sub_means$mean_RPM[sub_means$tsRNA_type == "3'-half"])/total)*100							
                                        							rownames(tRNA_sub_perc)[i] <- nam
                                        							}
                                        					idx <- rowMeans(tRNA_sub_sum_df) > filter
                                        						#Sys.sleep(0.1)
                                        						#message(noquote(paste("tRNA isotypes passed filter ", filter, ":", sep="")))
                                        						#Sys.sleep(0.1)
                                        						#print(c(table(idx)[1],table(idx)[2]))
                                        					tRNA_sub_sum_df <- tRNA_sub_sum_df[idx,]
                                        					tRNA_sub_perc <- tRNA_sub_perc[idx,]
                                        					tRNA_sub_sum_df_ord <- tRNA_sub_sum_df[order(rowMeans(tRNA_sub_sum_df), decreasing = TRUE),]
                                        					tRNA_sub_perc_ord <- tRNA_sub_perc[order(rowMeans(tRNA_sub_sum_df), decreasing = TRUE),]
                                        						stopifnot(identical(as.character(rownames(PAC$Pheno)), as.character(colnames(tRNA_sub_sum_df)))==TRUE)
                                        						stopifnot(identical(as.character(rownames(tRNA_sub_perc_ord)), as.character(rownames(tRNA_sub_sum_df_ord)))==TRUE)
                                        					RPM_means_df <- data.frame(variable=factor(rownames(tRNA_sub_sum_df_ord), levels=rownames(tRNA_sub_sum_df_ord)), means=rowMeans(tRNA_sub_sum_df_ord), SE=(apply(tRNA_sub_sum_df_ord,1,sd))/(sqrt(ncol(tRNA_sub_sum_df_ord))))
                                        					rownames(tRNA_sub_sum_df_ord) <- paste0("RPM_", rownames(tRNA_sub_sum_df_ord)) 
                                        					RPM_pheno <- cbind(PAC$Pheno[,c(1:2)], t(tRNA_sub_sum_df_ord)) 
                                        					### No diff and log2FC with independent samples 
                                        					if(paired==FALSE){
                                        					      list_wide <- list(RPM=RPM_pheno)
                                        					      list_long <- list(melt(RPM_pheno, id.vars= colnames(RPM_pheno)[!(grepl("RPM_", colnames(RPM_pheno)))]))
                                        					}
                                        					### Diff and log2FC with paired samples
                                        					if(paired==TRUE){
                                                  					log2FC <- log2(tRNA_sub_sum_df_ord[, PAC$Pheno[,pheno_targets[[1]]] == pheno_targets[[2]][1]] 
                                                  					                    / tRNA_sub_sum_df_ord[, PAC$Pheno[,pheno_targets[[1]]] == pheno_targets[[2]][2]])
                                                            rownames(log2FC) <-  gsub("RPM_", "log2FC_", rownames(log2FC))
                                                  						stopifnot(identical(colnames(log2FC), rownames(PAC$Pheno)[PAC$Pheno[, pheno_targets[[1]]] == pheno_targets[[2]][1]] ))
                                                  					
                                                  						log2FC_pheno <- cbind(PAC$Pheno[as.character(PAC$Pheno[,pheno_targets[[1]]]) %in% pheno_targets[[2]][1], c(1:2)], t(log2FC))
                                                  			list_wide <- list(RPM=RPM_pheno, log2FC=log2FC_pheno) # Save wide format in wide list
                                                  					RPM_pheno_long <- melt(RPM_pheno, id.vars= colnames(RPM_pheno)[!(grepl("RPM_", colnames(RPM_pheno)))])
                                                  					log2FC_pheno_long <- melt(log2FC_pheno, id.vars= colnames(log2FC_pheno)[!(grepl("log2FC_", colnames(log2FC_pheno)))])
                                                  			list_long <- list(RPM=RPM_pheno_long, log2FC=log2FC_pheno_long)
                                                  					}
                                        		list_final <- list(Long_formats=list_long, Wide_formats=list_wide, percent_tRF_type=tRNA_sub_perc_ord, Summary=RPM_means_df)
                                        		return(list_final)
                                        	}
                                ##-------- Function for plotting isodecoders (Mint)  
                                  plot.mint.isodec <- function(preped_list, paired, group_col, colfunc2, col_tRFtype){
                                                                            plot_lst <- list(NA)
                                                                            if(max(preped_list[[4]]$means)<100000){breaks <- c(1,10,100,1000,10000,100000)
                                                                            }else{breaks <- c(1,10,100,1000,10000,100000,1000000)}
                                                                            col_isotype <- colfunc2(length(unique(preped_list[[1]][[1]]$variable)))
                                                                         ## Log10 mean RPM (All)
                                                                            plot_df <- preped_list[[4]]
                                                                            plot_df$variable <- factor(plot_df$variable , levels=rev(levels(plot_df$variable)))
                                                                            plot_lst$log10_meanRPM <- ggplot(plot_df, aes(x=variable, y=means, fill=variable,
                                                                                              						  ymax = means + (means > 0)*SE,
                                                                                              						  ymin = means - (means < 0)*SE)) +
                                                                                              					geom_errorbar(width=0.5, size=1.0, colour="black") +
                                                                                              					geom_col(width = 0.9, cex=0.2, colour="black")+
                                                                                              					labs(title="Mean RPM")+
                                                                                              					ylab("Log10 RPM +/- SE") +
                                                                                              					scale_fill_manual(values=c(col_isotype))+
                                                                                              					theme_classic()+
                                                                                              					theme(legend.position="none", axis.title.y= element_blank(), panel.grid.major.y =  element_line(linetype="dashed", colour="grey", size=0.5), panel.grid.major.x = element_line(colour="grey", size=0.5), axis.text.x = element_text(angle = 0, hjust = 0), axis.text.y = element_text(angle = 0, hjust = 0), axis.line.x =element_blank())+
                                                                                              					scale_y_log10(limits = c(min(breaks),max(breaks)), breaks=breaks)+
                                                                                              					coord_flip()
                                                                        ## Percent filled bar (All)
                                                                            plot_df <- preped_list[[3]]
                                                                            plot_df$variable <- factor(rownames(plot_df), levels=rev(rownames(plot_df)))
                                                                            plot_df <- melt(plot_df, id.vars = "variable", variable.name = "type", value.name = "value")
                                                                            plot_df$type <- factor(plot_df$type , levels = c("tsRNA_5'", "half_5'", "tsRNA_i'", "tsRNA_3'", "half_3'"))
                                                                            plot_lst$Percent <- ggplot(plot_df, aes(x=variable, y=value, fill=type)) +
                                                                                                					#geom_hline(yintercept = 0, size=1, color="azure4")+
                                                                                                					geom_col(width = 0.9, cex=0.2, colour="black", position="fill")+
                                                                                                					labs(title="Mean percent tRF content")+
                                                                                                					ylab("%") +
                                                                                                					scale_fill_manual(values=col_tRFtype)+
                                                                                                					theme_classic()+
                                                                                                					theme(axis.title.y= element_blank(), axis.text.x = element_text(angle = 0, hjust = 0), axis.text.y = element_blank(), axis.line.x =element_blank(), axis.line.y =element_blank())+
                                                                                                					coord_flip(ylim=c(0, 1))
                                                                        if(paired==FALSE){
                                                                                    ## Difference independent groups 
                                                                                    plot_df <- preped_list[[1]][[1]]
                                                                                    if(is.null(group_col)){groups <- 1}
                                                                                    if(!is.null(group_col)){groups <- plot_df[,group_col]}
                                                                                    plot_df$variable <- factor(plot_df$variable , levels=rev(levels(plot_df$variable)))
                                                                                    plot_lst$Diff_RPM <- ggplot(plot_df, aes(x=variable, y=value, group=groups, fill=variable)) +
                                                                                                    					geom_hline(yintercept = 0, size=1, color="azure4")+
                                                                                                       					stat_summary(geom = "errorbar",  width=0.8, size=0.5, fun.data = mean_se, position = "dodge") +
                                                                                                    					  stat_summary(geom = "point", stroke=0.5, shape=21, size = 2.0, aes(fill=plot_df$variable), fun.y = mean, position = position_dodge(width=0.8)) +
                                                                                                    					labs(title=paste0("RPM ", paste0(rev(levels(plot_df[,group_col])), collapse = " then ")))+
                                                                                                    					ylab("Difference RPM +/- SE") + 
                                                                                                    					scale_fill_manual(values=rep(c("black", "white"), length(col_isotype)))+
                                                                                                    					scale_x_discrete(labels=gsub("difference_", "", levels(plot_df$variable)))+
                                                                                                    					theme_classic()+
                                                                                                    					theme(legend.position="none", axis.title.y= element_blank(), panel.grid.major.y =  element_line(linetype="dashed", colour="grey", size=0.5), panel.grid.major.x = element_line(colour="grey", size=0.5), axis.text.x = element_text(angle = 0, hjust = 0), axis.text.y = element_blank(), axis.line.x =element_blank(), axis.line.y =element_blank())+
                                                                                                              coord_flip()
                                                                                        
                                                                                    ## Difference independent groups log10
                                                                                    plot_df <- preped_list[[1]][[1]]
                                                                                    if(is.null(group_col)){groups <- 1}
                                                                                    if(!is.null(group_col)){groups <- plot_df[,group_col]}
                                                                                    plot_df$variable <- factor(plot_df$variable , levels=rev(levels(plot_df$variable)))
                                                                                    plot_lst$Diff_RPM_log10 <- ggplot(plot_df, aes(x=variable, y=value, group=groups, fill=variable)) +
                                                                                                    					geom_hline(yintercept = 0, size=1, color="azure4")+
                                                                                                       					stat_summary(geom = "errorbar",  width=0.8, size=0.5, fun.data = mean_se, position = "dodge") +
                                                                                                    					  stat_summary(geom = "point", stroke=0.5, shape=21, size = 2.0, aes(fill=plot_df$variable), fun.y = mean, position = position_dodge(width=0.8)) +
                                                                                                    					labs(title=paste0("RPM ", paste0(rev(levels(plot_df[,group_col])), collapse = " then ")))+
                                                                                                    					ylab("Difference RPM +/- SE") + 
                                                                                                    					scale_fill_manual(values=rep(c("black", "white"), length(col_isotype)))+
                                                                                                    					scale_x_discrete(labels=gsub("difference_", "", levels(plot_df$variable)))+
                                                                                                    					theme_classic()+
                                                                                                    					theme(legend.position="none", axis.title.y= element_blank(), panel.grid.major.y =  element_line(linetype="dashed", colour="grey", size=0.5), panel.grid.major.x = element_line(colour="grey", size=0.5), axis.text.x = element_text(angle = 0, hjust = 0), axis.text.y = element_blank(), axis.line.x =element_blank(), axis.line.y =element_blank())+
                                                                                              					      scale_y_log10(limits = c(min(breaks),max(breaks)), breaks=breaks)+
                                                                                                              coord_flip()
                                                                                                         }
                                                                        if(paired==TRUE){
                                                                                    ## Log2 fold change C vs B with jitter  
                                                                                    plot_df <- preped_list[[1]][[2]]
                                                                                    plot_df$variable <- factor(plot_df$variable , levels=rev(levels(plot_df$variable)))
                                                                                    plot_lst$Log2FC <- ggplot(plot_df, aes(x=variable, y=value, fill=variable)) +
                                                                                                      					geom_hline(yintercept = 0, size=1.5, color="azure4")+
                                                                                                                geom_jitter(aes(color= "grey"), position=position_jitter(0.1), cex=1.3) +
                                                                                                           				stat_summary(geom = "errorbar",  width=0.5, size=1.0, fun.data = mean_se, position = "identity") +
                                                                                                        					stat_summary(geom = "point", colour="black", stroke=1.5, shape=21, size = 3.5, fun.y = mean, position = "identity") +
                                                                                                      					labs(title="Log2 Fold change C vs B")+
                                                                                                      					ylab("Log2 Fold change +/- SE") +
                                                                                                      					scale_fill_manual(values=c(col_isotype))+
                                                                                                                scale_color_manual(values="bisque3") +
                                                                                                      					scale_x_discrete(labels=gsub("log2FC_", "", levels(plot_df$variable)))+
                                                                                                      					theme_classic()+
                                                                                                      					theme(legend.position="none", axis.title.y= element_blank(), panel.grid.major.y =  element_line(linetype="dashed", colour="grey", size=0.5), panel.grid.major.x = element_line(colour="grey", size=0.5), axis.text.x = element_text(angle = 0, hjust = 0), axis.text.y = element_blank(), axis.line.x =element_blank(), axis.line.y =element_blank())+
                                                                                                      					coord_flip(ylim=c(-1.3, 2.25))
                                                                                                       			
                                            
                                                                                    
                                                                                    
                                                                                            }
                                                                            return(plot_lst[-1])
                                  }
                                ##-------- Function for generating scaled differences of top candidates (Mint)  
                                  top.diff.scaling<- function(Mint_PAC,  center=TRUE, scale=TRUE){
                                                      ### Extract top candidates from PAC using info from Table S6
                                                      top_fragments <- Mint_PAC$Result[Mint_PAC$Result$"p value4" < 0.05,]
                                                      top_fragments <-top_fragments[rowMeans(cbind(top_fragments$Sugar, top_fragments$Healthy)) >=5,]
                                                      top_fragments <- top_fragments[grepl("\\<MT_SerTGA|\\<MT_ThrTGT|\\<MT_ValTAC|\\<MT_HisGTG|\\<LysCTT|\\<ArgCCG|\\<ArgCCT|\\<LeuCAA", top_fragments$Isodecoder2),]
                                                      idx <- as.character(Mint_PAC$Anno$"MINTbase ID") %in% top_fragments$ID
                                                      master_top <- Mint_PAC[c("Pheno","Anno", "RPM")]
                                                      master_top$Anno <- master_top$Anno[idx,]
                                                      master_top$RPM <- master_top$RPM[idx,]
                                                      stopifnot(identical(rownames(master_top$Anno), rownames(master_top$RPM)))
                                                      stopifnot(identical(rownames(master_top$Pheno), colnames(master_top$RPM)))              
                                                      ### Generate scaled (comparable) differences between Sugar and Healthy   
                                                      diff_ind <- (master_top$RPM[,master_top$Pheno$Diet=="Sugar"]) - (master_top$RPM[,master_top$Pheno$Diet=="Healthy"])
                                                      Scaled_differences <- t(as.matrix(scale(t(diff_ind), center=center, scale=scale)))
                                                      return(list(Info=data.frame(Center=center, Scaled=scale), Top_PAC=master_top, Top_results=top_fragments, Scaled_difference=Scaled_differences))
                                  }
                                ##-------- Function for plotting circular cladogram
                                  plot.mint.fan <- function(cluster_data, nclus=4, id_labels=FALSE) {
                                                            scaled_data <- cluster_data$Scaled_difference
                                                            labels <- gsub("SerTGA, |LeuCAG, ", "",  cluster_data$Top_PAC$Anno$Isodecoder)
                                                            color_iso <- data.frame(iso=c("LysCTT", "ArgCCG", "ArgCCT", "LeuCAA","MT_SerTGA",  "MT_ThrTGT", "MT_ValTAC", "MT_HisGTG"), col=colfunc_sports(8))
                                                            color_iso_lab <- color_iso[match(labels, color_iso$iso),]
                                                            HS_col <- color_iso_lab$col
                                                            if(id_labels==TRUE){
                                                                labels  <- paste(as.character(cluster_data$Top_PAC$Anno$'MINTbase ID'), labels)
                                                                }
                                                            dists <- dist(scaled_data)
                                                            hc <- hclust(dists)        # First prepare df with orignal data
                                                            clus_df<- cutree(hc, nclus)
                                                            #hc$labels <- color_iso$iso  		# Then change labels 
                                                            hc$labels <- labels 		# Then change labels 
                                                            palette <- c('#9D0014','#094A6B','#9D0014','#094A6B','#9D0014','#094A6B','#094A6B')[1:nclus]
                                                            clus  <- cutree(hc, nclus)
                                                            X <- as.phylo(hc)
                                                            edge.clus <- sapply(1:nclus,function(i)max(which(X$edge[,2] %in% which(clus==i))))
                                                            order     <- order(edge.clus)
                                                            edge.clus <- c(min(edge.clus), diff(sort(edge.clus)))
                                                            edge.clus <- rep(order, edge.clus)
                                                            plot(X, type='fan',
                                                                 tip.color=as.character(HS_col), edge.color=palette[edge.clus],
                                                            	 	edge.width= 4, label.offset=0.05, no.margin=TRUE, cex=1.4,
                                                            	 	show.tip.label=TRUE)
                                                            Sys.sleep(2)
                                                            p <- recordPlot()
                                                            pal <- palette[clus_df]
                                                            df <- data.frame(SEquence=names(clus_df), cluster=clus_df, col=pal)
                                                            return(list(plot=p, data=df))
                                  }
                            cat("\n-----------------------------------------------------------------------------------------\n")   
                            cat("Step 1: Setting up the environoment with packages and functions.\n")

                            ###############################################################################################################3
                            #### Reading data  ####
                            cat("\n-----------------------------------------------------------------------------------------\n")   
                            cat("Step 2: Reading and organizing data files from supplementary tables.\n") 
                                import_lst <- read.data(path)
                            ###############################################################################################################                            
                            #### Generate variables and plots for Sports ####
                            cat("\n-----------------------------------------------------------------------------------------\n")                                  
                            cat("Step 3: Sports analysis.\n")
                                  cat("---Step 3.1: Calculating RPM from Sports. Applying a low abundance filter to remove the worst noise (size= 16-45 bp, at least 0.01 PRM in all samples.\n")
                                    import_lst$Sports$RPM  <- generate.RPM(import_lst$Sports)
                                    import_lst$Sports_filtered  <- suppressMessages(RPM.filt(import_lst$Sports, size = c(16, 45), rpm = c(0.01, 100), pheno_targets= list(target_col="Diet", target_groups=c("Start", "Healthy", "Sugar"))))
                                  cat("---Step 3.2: Using RPM to calculate log2FC from Sports.\n")
                                    import_lst$Sports_filtered$Calc <- FC.bio.sports(import_lst$Sports_filtered, group_column="Diet", groups=c("Sugar", "Healthy"), paired=TRUE, noAnno.rm=TRUE)[-c(1:4)]
                                  cat("---Step 3.3: Plotting Sports data.\n")
                                      target.biotypes = c("rRNA", "tRNA", "Mt_tRNA", "miRNA", "piRNA", "lincRNA", "other")
                                      ##-------- Stacked barplot
                                      cat("-------Step 3.3.1: Fig 2A - Percent biotypes in stacked barplot of each participants.\n")
                                        import_lst$Plots$Sports$Fig2A_stacked_bar_plot <- plot.stacked.bar(import_lst$Sports_filtered$Calc$Sample_sumRPM_biotype, target.biotypes[c(7,1:6)], color_vec=rgb_vec_sports[c(7,1:6)])
                                      ##-------- Pie chart biotypes
                                      cat("-------Step 3.3.2: Fig 2B - Percent biotypes in pie chart of experiment means.\n")
                                        suppressMessages(plot.pie(import_lst$Sports_filtered$Calc$Experiment_sumRPM_biotype, target.biotypes,  color_vec=rgb_vec_sports, angle=-10))   
                                        Sys.sleep(2)
                                        import_lst$Plots$Sports$Fig2B_pie_chart_biotype <- recordPlot()
                                      ##-------- Jitterplot
                                      cat("-------Step 3.3.3: Fig 2D - Log2 fold changes as jitterplot of each small RNA.\n")
                                        import_lst$Sports_filtered$Calc$mean_Log2FC <- import_lst$Sports_filtered$Calc$mean_Log2FC[import_lst$Sports_filtered$Calc$mean_Log2FC$Biotype %in% as.character(target.biotypes), ]
                                        import_lst$Sports_filtered$Calc$mean_Log2FC$Biotype <- factor(as.character(import_lst$Sports_filtered$Calc$mean_Log2FC$Biotype), levels=as.character(target.biotypes))
                                        import_lst$Plots$Sports$Fig2D_jitterplot <- plot.jitter(import_lst$Sports_filtered$Calc$mean_Log2FC, FC_col="mean_FC", feat_col="Biotype", limits=c(-1.1, 1.05), Ypos_n=1.05, colors=rgb_vec_sports)                                   
                                      ##-------- Pie chart direction 
                                      cat("-------Step 3.3.4: Fig 2D - Percent up or down as pie plots of each biotype.\n")
                                         import_lst$Plots$Sports$Fig2D_pie_chart_direction <- suppressMessages(plot.perc.pie(import_lst$Sports_filtered$Calc$mean_Log2FC,  target.biotypes, colors=c("cornsilk4", "cornsilk3"), angle=0))
                            ###############################################################################################################                              
                            #### Generate variables and plots for Mint #####
                            cat("\n-----------------------------------------------------------------------------------------\n")           
                            cat("Step 4: Mintmap analysis.\n")           
                                  cat("---Step 4.1: countTable is already provided as RPM from Mintmap. A stricter filter was applied to remove more noise than in Sports (size= 16-45 bp, at least 1 PRM in 50% of samples.\n")
                                      import_lst$Mint$RPM <- import_lst$Mint$countTable 
                                      import_lst$Mint_filt <- suppressMessages(RPM.filt(import_lst$Mint, size = c(16, 45), rpm = c(1, 50), pheno_targets = list(target_col="Diet", target_groups=c( "Sugar", "Healthy"))))
                                  cat("---Step 4.2: tsRNA type analysis.\n")
                                    cat("------Step 4.2.1: Calculating values over tsRNA types.\n")
                                          import_lst$Mint_filt$Types$Data <- suppressMessages(prep.mint.types(import_lst$Mint_filt, paired=TRUE, pheno_target=list(column="Diet", groups=c("Sugar", "Healthy")), anno_target="tsRNA type"))
                                    cat("------Step 4.2.1: Fig 2F/G - Plotting tsRNA types.\n")
                                          import_lst$Plots$Mint$Fig2FG_types <- suppressMessages(plots.mint.types(import_lst$Mint_filt$Types$Data, color_vec=rgb_vec_mint_types))
                                  cat("---Step 4.3: tsRNA isodecoder analysis.\n")
                                      cat("------Step 4.3.1: Calculating values over tsRNA isodecoders, only considering isodecoders with >100 RPM.\n")
                                          import_lst$Mint_filt$Isodecoder$Data <- suppressMessages(prep.mint.isodec(import_lst$Mint_filt, filter=100, paired=TRUE, pheno_targets = list(target_col="Diet", target_groups=c("Sugar", "Healthy")), col.IDs="Participant"))
                                      cat("------Step 4.3.2: Fig 2H - Plotting tsRNA isodecoders.\n")
                                          import_lst$Plots$Mint$Fig2H_isodecoders <-  suppressMessages(plot.mint.isodec(import_lst$Mint_filt$Isodecoder$Data, paired=TRUE, group_col="Diet", colfunc_mint_isodecoders, rgb_vec_mint_types))
                                  cat("---Step 4.4: tsRNA cluster analysis of siginficant fragments in significant isodecoders.\n")
                                          cat("------Step 4.4.1: Extracting and scaling the differences between Sugar and Healthy of the top tsRNA in Table S6 (p<0.05; mean RPM >=5).\n")
                                          import_lst$Mint$Cluster <- top.diff.scaling(import_lst$Mint)
                                          cat("------Step 4.4.2: Fig 3 - Plotting circular cladogram.\n")
                                          import_lst$Plots$Mint$Fig3_cluster_cladogram <- plot.mint.fan(import_lst$Mint$Cluster, nclus=4)
                              return(import_lst[c(4,1,3,2,5)])     
                              }

### Run function on the Supplementary_Tables_S1-S6.xlsx file:                                  
results <- generate.results(path)


###########################################################################################################################################
##--------------------------------------------------------------------------------##
####       Step 5: Getting to know the structure in results 
##--------------------------------------------------------------------------------##
## Show results structure results list:
str(results, max.level =2, give.attr = FALSE)

## Notes:
## All plots have been saved in "Plots".
## Raw, filted and calculated data are stored in Sport and Mint respective "folders".

## Show structure of the Plots folder:
str(results$Plots, max.level =2, give.attr = FALSE)

## You can either use the plot_grid function in the cowplot package to plot side by side.
## Or you can plot individually graphs.  




###########################################################################################################################################
##--------------------------------------------------------------------------------##
####       Step 6: Plotting the graphs 
##--------------------------------------------------------------------------------##

require(cowplot)
###-------------------------------------------------------------------
## Sports
#
# Fig 2A:
results$Plots$Sports[[1]]

# Fig 2B:
results$Plots$Sports[[2]] 

# Fig 2D:
plot_grid(results$Plots$Sports[[4]], results$Plots$Sports[[3]]$input, labels = NULL, hjust = 0, vjust = 1, ncol = 1, nrow = 2,  rel_widths = c(1, 1), rel_heights = c(0.3, 1)) 
          

###-------------------------------------------------------------------
## Mint
# Fig 2E and Fig 2F:
plot_grid(results$Plots$Mint[[1]][[1]], results$Plots$Mint[[1]][[2]], labels = NULL, hjust = 0, vjust = 1, ncol = 1, nrow = 2,  rel_widths = c(1, 1), rel_heights = c(1, 1)) 

# Fig 2H:
plot_grid(results$Plots$Mint[[2]][[1]], results$Plots$Mint[[2]][[3]], results$Plots$Mint[[2]][[2]], labels = NULL, hjust = 0, vjust = 1, ncol = 3, nrow = 1,  rel_widths = c(1, 1), rel_heights = c(1, 1)) 

# Fig 3
results$Plots$Mint[[3]][[1]] 





###########################################################################################################################################
##--------------------------------------------------------------------------------##
####       SessionInfo
##--------------------------------------------------------------------------------##
#
# sessionInfo()
# 
# R version 3.4.4 (2018-03-15)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 19.1
# 
# Matrix products: default
# BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] ape_5.3            dendextend_1.12.0  cowplot_0.9.4      extrafont_0.17     ggthemes_4.2.0     scales_1.0.0       RColorBrewer_1.1-2 reshape2_1.4.3     readxl_1.3.1       ggplot2_3.1.1     
# 
# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.1        cellranger_1.1.0  pillar_1.4.1      compiler_3.4.4    plyr_1.8.4        viridis_0.5.1     tools_3.4.4       digest_0.6.19     nlme_3.1-131      lattice_0.20-35   tibble_2.1.2      gtable_0.3.0      viridisLite_0.3.0 pkgconfig_2.0.2  
# [15] rlang_0.3.4       cli_1.1.0         rstudioapi_0.10   parallel_3.4.4    gridExtra_2.3     Rttf2pt1_1.3.7    withr_2.1.2       dplyr_0.8.1       stringr_1.4.0     grid_3.4.4        tidyselect_0.2.5  glue_1.3.1        R6_2.4.0          purrr_0.3.2      
# [29] extrafontdb_1.0   magrittr_1.5      assertthat_0.2.1  colorspace_1.4-1  labeling_0.3      stringi_1.4.3     lazyeval_0.2.2    munsell_0.5.0     crayon_1.3.4  
#             
                                
                                