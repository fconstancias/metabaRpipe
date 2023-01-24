#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'
#'check_primers(path = "~/Documents/GitHub/interfaces/data/raw/agro/raw/ERR5194893/", FWD = reverseComplement(DNAString("GTGCCAGCMGCCGCGGTAA")), REV =  (DNAString("GGACTACHVGGGTWTCTAAT")))

check_primers <- function(path_dir,
                          n_samples = 1,
                          FWD,
                          REV,
                          sep = "[^_]+",
                          export = FALSE){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(ShortRead); require(Biostrings)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using ShortRead version ", packageVersion('ShortRead'),'\n\n'))
  cat(paste0('\n##',"You are using Biostrings version ", packageVersion('Biostrings'),'\n\n'))
  
  cat('################################\n\n')
  
  
  ## ------------------------------------------------------------------------
  # generate com rev comp primer sequences
  
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  FWD.orients
  
  ## ------------------------------------------------------------------------
  # define function to get promer hits
  primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  ## ------------------------------------------------------------------------
  # get run directories
  
  ## get run directories (under the raw_files_path)
  list.dirs(path = path_dir, 
            full.names = FALSE, 
            recursive = FALSE) %>% as.vector() -> run_list  
  
  out <- vector("list", length(run_list)) 
  names(out) <- run_list
  
  ## ------------------------------------------------------------------------
  
  for(i in seq_along(run_list)) {
    
    Fs <- sort(list.files(file.path(path_dir, 
                                    run_list[i]), 
                          pattern = glob2rx(as.character(file_pattern[1])), 
                          full.names = TRUE))
    
    Rs <- sort(list.files(file.path(path_dir, 
                                    run_list[i]), 
                          pattern = glob2rx(as.character(file_pattern[2])), 
                          full.names = TRUE))
    
    
    exists <- file.exists(Fs) & file.exists(Rs)
    
    Fs <- Fs[exists]
    Rs <- Rs[exists]
    
    # readFastq(fnRs) %>% #idea to filter based on length
    #   ShortRead::sread()
    
    if(length(Rs) != length(Fs)) stop ("Forward and reverse files do not match.")
    
    sample.names <- basename(Fs) %>%
      str_extract(paste0(sep))
    
    cat(paste0('\n# sample names list starts with : '))
    head(sample.names)
    
    ## ------------------------------------------------------------------------
    
    # set.seed(seed_value) #random  generator necessary for reproducibility
    # 
    # ii <- sample(length(sample.names),
    #              round(length(sample.names) * (prop.sample/100),0) + 1 ) 
    
    
    ## ------------------------------------------------------------------------
    
    rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = Fs[[n_samples]]), 
          FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = Rs[[n_samples]]), 
          REV.ForwardReads = sapply(REV.orients, primerHits, fn = Fs[[n_samples]]), 
          REV.ReverseReads = sapply(REV.orients, primerHits, fn = Rs[[n_samples]])) -> out[[i]]
    
  }
  
  if(export != FALSE){
    
    dir.create(export, recursive = TRUE)
    
    out %>%
      do.call(rbind.data.frame, .)  %>%
      data.frame() %>%
      rownames_to_column("id") %>%
      write_tsv(file = paste0(export, "primers_check.tsv"))
  }
  return(out)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#' run_atropos()
#'

run_atropos <- function(raw_files_path,
                        atropos = "atropos",
                        cut_dir = "dada2/00_atropos_primer_removed",
                        PRIMER_F = "CCTAYGGGRBGCASCAG",
                        PRIMER_R = "GGACTACNNGGGTATCTAAT",
                        NSLOTS = 6,
                        raw_file_pattern = c("*_R1_*.gz","*_R2_*.gz"),
                        cut_file_pattern = c("_primersout_R1_.fastq.gz","_primersout_R2_.fastq.gz"),
                        MIN_L = 100,
                        sep = "[^_]+"){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(ShortRead); require(Biostrings)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using ShortRead version ", packageVersion('ShortRead'),'\n\n'))
  cat(paste0('\n##',"You are using Biostrings version ", packageVersion('Biostrings'),'\n\n'))
  cat(paste0('\n##',"You are using atropos version "));  system2(atropos, args = "trim --version")
  
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  list.dirs(path = raw_files_path, 
            full.names = FALSE, 
            recursive = FALSE) %>% as.vector() -> run_list  
  
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  
  # setwd(raw_files_path)
  # setwd("./..")
  
  for(i in seq_along(run_list)) {
    
    
    cut_path <- file.path(cut_dir, 
                          run_list[i])
    
    dir.create(cut_path, showWarnings = TRUE, recursive = TRUE)
    
    cat(paste0('\n# output dir :  ',cut_path,'\n'))
    
    
    if(length(fnRs) != length(fnFs)) stop("Forward and reverse files do not match.")
    
    if(length(fnRs) == 0) stop("No raw fastq files detected, please check...")
    
    
    sample.names <- basename(fnFs) %>%
      str_extract(paste0(sep))
    
    cat(paste0('\n# sample names list starts with : \n'))
    head(sample.names)
    
    fnFs_cut <- file.path(cut_path, paste0(sample.names, cut_file_pattern[1]))
    fnRs_cut <- file.path(cut_path, paste0(sample.names, cut_file_pattern[2]))
    
    sum(nchar(PRIMER_F),nchar(PRIMER_R))/2 * 2/3 -> MIN_F
    
    ## ------------------------------------------------------------------------
    
  }
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'


run_dada2_qplot <- function(prop.sample = 20,
                            aggregate = TRUE,
                            cut_dir = "dada2/00_atropos_primer_removed",
                            qplot_dir = "dada2/01_dada2_quality_profiles",
                            cut_file_pattern = c("*_R1_*.fastq.gz","*_R2_*.fastq.gz"),
                            seed_value = 123,
                            sep = "[^_]+",
                            export = TRUE){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------
  ## get run directories (under the raw_files_path)
  list.dirs(path = cut_dir, 
            full.names = FALSE, 
            recursive = FALSE) %>% as.vector() -> run_list  
  
  out_fwd <- vector("list", length(run_list)) # Fwd and Rev
  names(out_fwd) = c(run_list)
  out_rev <- vector("list", length(run_list)) # Fwd and Rev
  names(out_rev) = c(run_list)
  
  ## ------------------------------------------------------------------------
  
  
  for(i in seq_along(run_list)) {
    
    cut_path <- file.path(cut_dir, 
                          run_list[i])
    
    out_qplot <- file.path(qplot_dir, 
                           run_list[i])
    
    cat(paste0('\n# output dir :  ',out_qplot,'\n'))
    
    
    dir.create(out_qplot, showWarnings = TRUE, recursive = TRUE)
    
    
    fnFs_cut <- sort(list.files(cut_path, 
                                pattern = glob2rx(paste0("*", cut_file_pattern[1])), 
                                full.names = TRUE))
    fnRs_cut <- sort(list.files(cut_path, 
                                pattern = glob2rx(paste0("*", cut_file_pattern[2])), 
                                full.names = TRUE))
    
    
    exists <- file.exists(fnFs_cut) & file.exists(fnRs_cut)
    
    fnFs <- fnFs_cut[exists]
    fnRs <- fnRs_cut[exists]
    
    # readFastq(fnRs) %>% #idea to filter based on length
    #   ShortRead::sread()
    
    if(length(fnRs) != length(fnFs)) stop ("Forward and reverse files do not match.")
    
    sample.names <- basename(fnFs) %>%
      str_extract(paste0(sep))
    
    cat(paste0('\n# sample names list starts with : \n'))
    cat(sample.names)
    
    ## ------------------------------------------------------------------------
    
    
    ## ------------------------------------------------------------------------
    ### Forward
    set.seed(seed_value) #random  generator necessary for reproducibility
    
    ii <- sample(length(sample.names),
                 round(length(sample.names) * (prop.sample/100),0)+ 1 ) 
    
    qplot_f <- plotQualityProfile(fnFs[ii], aggregate = as.logical(aggregate))
    
    out_fwd[[i]] <- qplot_f +
      geom_vline(xintercept = lines_seq(max(qplot_f$data$Cycle)), alpha = 0.3, size = 0.3, linetype = 'dashed') +
      improve_headers +
      ggtitle(paste0(run_list[i]," - ","Forward reads"))
    
    ## ------------------------------------------------------------------------
    ### Reverse
    
    qplot_r <- plotQualityProfile(fnRs[ii], aggregate = aggregate)
    
    out_rev[[i]] <- qplot_r +
      geom_vline(xintercept = lines_seq(max(qplot_r$data$Cycle)), alpha = 0.3, size = 0.3, linetype = 'dashed') +
      improve_headers +
      ggtitle(paste0(run_list[i]," - ","Reverse reads"))
    
    
    if (export == TRUE){
      ggsave(plot = out_fwd[[i]],
             path= out_qplot,
             device="pdf",
             filename = paste0(run_list[i],"_forward.pdf"),
             width = 410,
             height = 300,
             units = 'mm')
      
      ggsave(plot = out_rev[[i]] ,
             path= out_qplot,
             device="pdf",
             filename = paste0(run_list[i],"_reverse.pdf"),
             width = 410,
             height = 300,
             units = 'mm')
    }
    
  }
  out <- list("fwd_plot" = out_fwd,
              "rev_plot" = out_rev)
  return(out)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'

run_dada2_filter_denoise_merge_reads <- function(trunclen,
                                                 maxee,
                                                 truncQ = 6,
                                                 minLen = 100,
                                                 nthreads = 6,
                                                 maxLen = Inf,
                                                 nbases = 20000000,
                                                 pool = "pseudo",
                                                 priors = "FALSE",
                                                 minover = 12,
                                                 cut_dir = "dada2/00_atropos_primer_removed",
                                                 filt_dir = "dada2/02_dada2_filtered_denoised_merged",
                                                 cut_file_pattern = c("_primersout_R1_.fastq.gz","_primersout_R2_.fastq.gz"),
                                                 filt_pattern = c("_R1_filtered.fastq.gz","_R2_filtered.fastq.gz"),
                                                 sep = "[^_]+",
                                                 seed_value = 123,
                                                 remove_input_fastq = TRUE,
                                                 export = TRUE)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------
  ## get run directories (under the cut_dir)
  list.dirs(path = cut_dir, 
            full.names = FALSE, 
            recursive = FALSE) %>% as.vector() -> run_list  
  
  seqtab <- vector("list", length(run_list)) # Fwd and Rev
  names(seqtab) = c(run_list)
  
  plot <- vector("list", length(run_list)) # Fwd and Rev
  names(plot) = c(run_list)
  
  track <- vector("list", length(run_list)) # Fwd and Rev
  names(track) = c(run_list)
  
  out_fwd <- vector("list", length(run_list)) # Fwd and Rev
  names(out_fwd) = c(run_list)
  
  out_rev <- vector("list", length(run_list)) # Fwd and Rev
  names(out_rev) = c(run_list)
  
  ## ------------------------------------------------------------------------
  
  for(i in seq_along(run_list)) {
    
    
    cut_path <- file.path(cut_dir,
                          run_list[i])
    
    filt_path <- file.path(filt_dir,
                           run_list[i])
    
    
    dir.create(filt_path, showWarnings = TRUE, recursive = TRUE)
    
    cat(paste0('\n# output dir :  ',filt_path,'\n'))
    
    fnFs_cut <- sort(list.files(cut_path, pattern = glob2rx(paste0("*", cut_file_pattern[1])), full.names = TRUE))
    fnRs_cut <- sort(list.files(cut_path, pattern = glob2rx(paste0("*", cut_file_pattern[2])), full.names = TRUE))
    
    exists <- file.exists(fnFs_cut) & file.exists(fnRs_cut)
    
    fnFs <- fnFs_cut[exists]
    fnRs <- fnRs_cut[exists]
    
    # readFastq(fnRs) %>% #idea to filter based on length
    #   ShortRead::sread()
    
    if(length(fnRs) != length(fnFs)) stop ("Forward and reverse files do not match.")
    
    sample.names <- basename(fnFs) %>%
      str_extract(paste0(sep))
    
    cat(paste0('\n# sample names list starts with : \n'))
    cat(sample.names)
    
    if(sample.names %>% length()  <= 1) stop("There is something wrong with the provided samples or you just provided one sample.")
    
    ## ------------------------------------------------------------------------
    filtFs <- file.path(filt_path, paste0(sample.names, filt_pattern[1]))
    filtRs <- file.path(filt_path, paste0(sample.names, filt_pattern[2]))
    
    ## ------------------------------------------------------------------------
    
    cat(paste0('\n# filterAndTrim \n'))
    
    # error when multi = T and loads of samples ?
    # https://github.com/benjjneb/dada2/issues/273
    # ?filterAndTrim
    # If memory is an issue, execute in a clean environment and reduce the chunk size n and/or the number of threads.
    # now OMP = TRUE and multi = FALSE
    
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=trunclen, trimLeft = 0, trimRight = 0,
                         maxEE=maxee, maxLen = maxLen,  rm.phix=TRUE, maxN=0, minLen=minLen, verbose = T,
                         compress=TRUE, multithread= nthreads, truncQ = truncQ, n = 1e5, OMP = FALSE)
    
    out2 <- data.frame(row.names = sample.names,
                       out)
    
    cat('\n# Filtering and trimming done with the following parameters:')
    cat(str_c('\n# Forward pair: trimming at ',trunclen[1],' nts and max expected error ',maxee[1]))
    cat(str_c('\n# Reverse pair: trimming at ',trunclen[2],' nts and max expected error ',maxee[2],'\n\n'))
    
    ## ------------------------------------------------------------------------
    
    cat(str_c('\n# Filtered fastq files were generated in : "', filt_path,'" \n'))
    
    ## ------------------------------------------------------------------------
    
    # get filtered reads if they still exist after trimming
    filtFs <- file.path(filt_path, paste0(sample.names, filt_pattern[1]))
    filtRs <- file.path(filt_path, paste0(sample.names, filt_pattern[2]))
    
    exists <- file.exists(filtFs) & file.exists(filtRs)
    filtFs <- filtFs[exists]
    filtRs <- filtRs[exists]
    
    if(length(filtFs) != length(filtRs)) stop("Forward and reverse filtered files do not match.")
    
    ## ------------------------------------------------------------------------
    cat(str_c('\n# derepFastq \n'))
    
    derepFs <- derepFastq(filtFs, verbose=FALSE)
    derepRs <- derepFastq(filtRs, verbose=FALSE)
    
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names[exists]
    names(derepRs) <- sample.names[exists]
    
    cat('\n# Dereplication done\n')
    
    ## ------------------------------------------------------------------------
    set.seed(seed_value) #random  generator necessary for reproducibility
    cat(str_c('\n# learnErrors  \n'))
    
    
    
    ## ------------------------------------------------------------------------
    cat(paste0('\n# dada \n'))
    set.seed(seed_value) #random  generator necessary for reproducibility
    
    if(identical(names(derepFs), names(derepRs)) != TRUE ) stop("Samples names are not consistent between Forward and Reverse dereplicated samples")
    
    if(priors != "FALSE"){
      priors %>% 
        Biostrings::readDNAStringSet() %>%  data.frame() %>%  pull(".") -> priors_seq}
    
    if(priors == "FALSE"){
      priors_seq = character(0)}
    
    dadaFs <- dada(derepFs, 
                   err=errF, 
                   multithread= nthreads, 
                   pool=ifelse(pool == "FALSE", as.logical(pool), pool),
                   priors = priors_seq)
    
    dadaRs <- dada(derepRs, 
                   err=errR, 
                   multithread= nthreads, 
                   pool=ifelse(pool == "FALSE", as.logical(pool), pool),
                   priors = priors_seq)
    
    
    cat('\n# DADA2 algorithm performed \n')
    
    # https://github.com/benjjneb/dada2/issues/77
    # dada2:::checkConvergence(dadaFs[[1]])
    # dada2:::checkConvergence(dadaRs[[1]])
    
    ## ------------------------------------------------------------------------
    cat(str_c('\n# mergePairs with defined ',minover ,' nt overlap \n'))
    set.seed(seed_value) #random  generator necessary for reproducibility
    
    if(sample.names %>% length()  > 1){
      if(identical(names(dadaFs), names(dadaRs)) != TRUE ) stop("Samples names are not consistent between Forward and Reverse dadas")
      if(identical(names(derepFs), names(derepRs)) != TRUE ) stop("Samples names are not consistent between Forward and Reverse dereplicated sequences")
    }
    
    
    mergers <- mergePairs(dadaFs, derepFs,
                          dadaRs, derepRs,
                          minOverlap = minover,
                          justConcatenate = FALSE,
                          maxMismatch = 0)
    
    cat('# Pairs were merged\n')
    
    ## ------------------------------------------------------------------------
    seqtab[[i]] <- makeSequenceTable(mergers)
    
    cat(paste0('# Number of samples: ',dim(seqtab[[i]])[1], '\n'))
    cat(paste0('# Number of detected variants (ASVs): ',dim(seqtab[[i]])[2]))
    cat("# The variants (ASVs) have the following length distribution:")
    table(nchar(getSequences(seqtab[[i]]))) # plot in the future ?
    
    cat(str_c('\n# saving seqtab as ',str_c(filt_path,"/",run_list[i],"_seqtab.rds") ,'\n'))
    
    if( export ==TRUE){
      saveRDS(seqtab[[i]], str_c(filt_path,"/",run_list[i],"_seqtab.rds"))
      
    }
    
    ## ------------------------------------------------------------------------
    cat(str_c('\n# Plotting Sequences/ASV distribution to ',str_c(filt_path, "/", "seq_distrib_",run_list[i],".pdf") ,'\n\n'))
    
    plotLengthDistro <- function(st) {
      tot.svs <- table(nchar(colnames(st)))
      tot.reads <- tapply(colSums(st), nchar(colnames(st)), sum)
      df <- data.frame(Length=as.integer(c(names(tot.svs), names(tot.reads))),
                       Count=c(tot.svs, tot.reads),
                       Type=rep(c("ASVs", "Reads"), times=c(length(tot.svs), length(tot.reads))))
      p <- ggplot(data=df, aes(x=Length, y=Count, color=Type)) + geom_point() + 
        facet_wrap(~Type, scales="free_y") + theme_bw() + xlab("Amplicon Length")
      
      return(p)
    }
    
    plotLengthDistro(seqtab[[i]]) + scale_y_log10() + 
      ggtitle(str_c("Sequence / ASV length distribution : ",run_list[i], " Run")) -> plot[[i]]
    
    if( export ==TRUE){
      ggsave(str_c(filt_path, "/","seq_distrib_",run_list[i],".pdf"), plot=plot[[i]], width = 9, height = 8)
    }
    ## ------------------------------------------------------------------------
    cat(str_c('\n# Generating summary \n'))
    
    
    
    if( export ==TRUE){
      write_tsv(as.data.frame(track[[i]]),str_c(filt_path,"/",run_list[i],"_track_analysis.tsv"))
      
      save(track, 
           seqtab,
           file=paste0(filt_path,"/",run_list[i],".RData"))
    }
    
    cat("\n# The distribution of merged kept is the following:\n")
    summary(track[[i]]$merged_pc)
    
    # assign(paste0(name.run,"_track"), track)
    # assign(paste0(name.run,"_seqtab"), seqtab)
    
    
    #load(paste0(output,"/",name.run,".RData"))
    
    file.remove(filtFs, filtRs)
    
    if( remove_input_fastq ==TRUE){
      file.remove(fnFs_cut, fnRs_cut)
      
    }
    
  }
  
  out <- list("track" = track,
              "seqtab" = seqtab,
              "plot" = plot,
              "out_fwd" = out_fwd,
              "out_rev" = out_rev)
  
  return(out)
}



#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note TODO: directly export phyloseq object at this point...
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'

run_dada2_mergeRuns_removeBimeraDenovo <- function(seqtab = NULL,
                                                   track = NULL,
                                                   merged_run_dir = "dada2/03_dada2_merged_runs_chimera_removed",
                                                   chimera_method = "consensus",
                                                   trim_length,
                                                   nthreads = 6,
                                                   collapseNoMis = FALSE,
                                                   minOverlap = 20, # https://github.com/benjjneb/dada2/issues/518
                                                   filt_dir = "dada2/02_dada2_filtered_denoised_merged",
                                                   export = TRUE,
                                                   seed_value = 123,
                                                   return = TRUE){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(phyloseq)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))
  
  cat('################################\n\n')
  
  
  ## ------------------------------------------------------------------------
  # Distribution of variants
  cat("\n# The variants (ASVs) have the following length distribution:\n")
  table(nchar(getSequences(seqtab.raw)))
  
  cat(paste0('\n# Reads shorter than ',trim_length[1],'bp and longer than ',trim_length[2], 'bp are going to be removed.\n'))
  
  
  plotLengthDistro <- function(st) {
    tot.svs <- table(nchar(colnames(st)))
    tot.reads <- tapply(colSums(st), nchar(colnames(st)), sum)
    df <- data.frame(Length=as.integer(c(names(tot.svs), names(tot.reads))),
                     Count=c(tot.svs, tot.reads),
                     Type=rep(c("ASVs", "Reads"), times=c(length(tot.svs), length(tot.reads))))
    
    p <- ggplot(data=df, aes(x=Length, y=Count, color=Type)) + geom_point() + facet_wrap(~Type, scales="free_y") + theme_bw() + xlab("Amplicon Length")
    
    return(p)
    
  }
  
  plotLengthDistro(seqtab.raw) + scale_y_log10() + 
    ggtitle(str_c("Overall Sequence / ASV length distribution ")) +
    geom_vline(xintercept = trim_length[1], size = 0.1, colour = "red", alpha = 0.8, linetype = 2, show.legend = "min") +
    geom_vline(xintercept = trim_length[2], size = 0.1, colour = "red", alpha = 0.8, linetype = 2 ) -> plot
  
  
  if (export == TRUE){
    ggsave(str_c(merged_run_dir,"/","seqtab_distrib",".pdf"),plot=plot, width = 9, height = 8)
  }
  ## ------------------------------------------------------------------------
  
  
  # Trim the unespecific amplifications from our dataset
  seqtab <- seqtab.raw[,nchar(colnames(seqtab.raw)) %in% seq(trim_length[1],
                                                             trim_length[2])]
  
  
  cat(paste0('\n# Reads shorter than ',trim_length[1],'bp and longer than ',trim_length[2], 'bp were removed.\n'))
  cat("\n# The variants (ASVs) after length filtering have the following length distribution:\n")
  table(nchar(getSequences(seqtab)))
  
  
  cat(paste0("\n# A total of ", round((sum(colSums(seqtab)) * 100) / sum(colSums(seqtab.raw)), digits = 2), "% reads were kept after length filtering.\n\n"))
  cat(paste0("\n# A total of ", round((dim(seqtab)[2] * 100) / dim(seqtab.raw)[2], digits = 2), "% ASVs were kept after length filtering.\n\n"))
  
  
  ## ------------------------------------------------------------------------
  full_join(summary, 
            data.frame(sample = rownames(st.all),
                       tabled_joined = rowSums(st.all),
                       chimera_out = rowSums(seqtab.raw),
                       length_filtered = rowSums(seqtab)), by='sample') %>% 
    mutate(tabled_pc = round(tabled_joined /merged, 2)) %>%
    mutate(chimera_out_pc = round(chimera_out/tabled, 2)) %>% 
    mutate(length_filtered_pc = round(length_filtered/chimera_out, 2)) -> track
  
  #write_tsv(data.frame(track),str_c(output,"/",name,"_track_analysis.tsv"))
  
  cat("\n# The distribution of chimera reads kept is the following:\n")
  summary(track$chimera_out_pc)
  
  
  ## ------------------------------------------------------------------------
  if(collapseNoMis==TRUE){
    
    if (export == TRUE){
      cat(str_c('\n# Saving uncollapsed .rds and fasta files as well as summary .tsv \n'))
      
      
      saveRDS(seqtab, str_c(merged_run_dir,"/uncollapsed_no-chim-seqtab.rds"))
      
      uniquesToFasta(seqtab,
                     str_c(merged_run_dir,"/uncollapsed_no-chim-seqtab.fasta"),
                     ids= str_c("asv",c(1:ncol(seqtab)), ";size=", colSums(seqtab)))
      
      write_tsv(track, str_c(merged_run_dir,"/uncollapsed_track_analysis.tsv"))
    }
    cat('\n# You have decided to run collapseNoMismatch on your dataset. Please note that it is only helpful IF you are working with several sequencing runs and it might take a long time to run. You might want to go further (taxonomy, ...) on your uncollapsed seqtable while it is running \n')
    
    collapsed_100 <- collapseNoMismatch(seqtab,  
                                        minOverlap = minOverlap,
                                        identicalOnly = FALSE)
    
    cat(str_c('\n# Saving collapsed .rds and fasta files as well as summary .tsv  \n'))
    
    physeq <- merge_phyloseq(otu_table(t(collapsed_100), taxa_are_rows=TRUE))
    
    ASV_seq <- Biostrings::DNAStringSet(taxa_names(physeq))
    names(ASV_seq) <- taxa_names(physeq)
    
    physeq <- merge_phyloseq(physeq, 
                             ASV_seq)
    
    taxa_names(physeq) <- paste0("ASV", str_pad(seq(ntaxa(physeq)),
                                                nchar(ntaxa(physeq)),
                                                pad = "0"))
    
    if (export == TRUE){
      saveRDS(collapsed_100, str_c(merged_run_dir,"/minOverlap_",minOverlap,"_collapse_no_mismatch_no-chim-seqtab.rds"))
      
      saveRDS(physeq, str_c(merged_run_dir,"/physeq.rds"))
      
      uniquesToFasta(collapsed_100,
                     str_c(merged_run_dir,"/minOverlap_",minOverlap,"_collapse_no_mismatch_no-chim-seqtab.fasta"),
                     ids= str_c("asv",c(1:ncol(collapsed_100)), ";size=", colSums(collapsed_100)))
    }
    
    track %>% 
      left_join(data.frame(sample = rownames(collapsed_100),
                           collapsed_100 = rowSums(collapsed_100)),
                by='sample') %>% 
      mutate(collapsed_100_pc = round(collapsed_100 / length_filtered, digits = 10)) -> track.final
    
    
    if (export == TRUE){
      
      write_tsv(track.final, str_c(merged_run_dir,"/track_analysis.tsv"))
      
      cat(str_c('# Your final 100% clustered ASV table can be found in "', paste0(merged_run_dir,"/seqtab.rds"),'"\n'))
      cat(str_c('# A FASTA file with your final ASVs was written in "',paste0(merged_run_dir,"/seqtab.fasta"), '"\n'))
      
      cat(str_c('# In "',paste0(merged_run_dir,"/track_analysis.tsv"),"\" you will find a table where you can check the loss of reads in each step. Check it out to see if everything's correct!",'\n'))
      # cat(str_c('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',str_c(output,"/",name,"_track_analysis.tsv "),'." and go further.\n'))
    }
  }
  ## ------------------------------------------------------------------------
  if(collapseNoMis==FALSE){
    
    physeq <- merge_phyloseq(otu_table(t(seqtab), taxa_are_rows=TRUE))
    ASV_seq <- Biostrings::DNAStringSet(taxa_names(physeq))
    names(ASV_seq) <- taxa_names(physeq)
    
    physeq <- merge_phyloseq(physeq, 
                             ASV_seq)
    
    taxa_names(physeq) <- paste0("ASV", str_pad(seq(ntaxa(physeq)),
                                                nchar(ntaxa(physeq)),
                                                pad = "0"))
    
    if (export == TRUE){
      cat(str_c('\n# Saving .rds and fasta files as well as summary .tsv \n'))
      
      saveRDS(seqtab, str_c(merged_run_dir,"/no-chim-seqtab.rds"))
      
      uniquesToFasta(seqtab,
                     str_c(merged_run_dir,"/no-chim-seqtab.fasta"),
                     ids= str_c("asv",c(1:ncol(seqtab)), ";size=", colSums(seqtab)))
      
      write_tsv(track, str_c(merged_run_dir,"/track_analysis.tsv"))
      
      saveRDS(physeq, str_c(merged_run_dir,"/physeq.rds"))
      
      cat(str_c('# Your final ASV table can be found in "', paste0(merged_run_dir,"/seqtab.rds"),'"\n'))
      cat(str_c('# A FASTA file with your final ASVs was written in "',paste0(merged_run_dir,"/seqtab.fasta"), '"\n'))
      
      cat(str_c('# In "',paste0(merged_run_dir,"/track_analysis.tsv"),"\" you will find a table where you can check the loss of reads in each step. Check it out to see if everything's correct!",'\n'))
      # cat(str_c('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',str_c(output,"/",name,"_track_analysis.tsv "),'." and go further.\n'))
    }
    cat('\n# chimera removal step is done. You can go further into the analysis (taxonomy, phylogeny) or explore collapseNoMismatch IF you are dealing with multiple runs ... it might take very long \n\n')
  }
  ## ------------------------------------------------------------------------
  
  if (return == TRUE){
    if(collapseNoMis==TRUE){
      
      return(list("seqtab" = collapsed_100,
                  "track" = track,
                  "plot" = plot,
                  "physeq" = physeq))
    }else{
      return(list("seqtab" = seqtab,
                  "track" = track,
                  "plot" = plot,
                  "physeq" = physeq))
    }
  }
  
  
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'



run_dada_taxonomy <- function(seqtab = NULL,
                              taxa_dir = "dada2/04_dada2_taxonomy",
                              threshold = 60,  # used for DECIPHER and dada2 if outputBootstraps = FALSE
                              tryRC = FALSE,
                              collapseNoMis = FALSE,
                              db,
                              db_species,
                              nthreads = 6,
                              merged_run_dir = "dada2/03_dada2_merged_runs_chimera_removed",
                              export = TRUE,
                              seed_value = 123)
{
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  
  if(is.null(seqtab)){
    
    taxa_path <- file.path(taxa_dir)
    
    dir.create(taxa_path, showWarnings = TRUE, recursive = TRUE)
    
    cat(paste0('\n# output dir :  ',taxa_path,'\n'))
    
    ## ------------------------------------------------------------------------
    
    merged_run_path <- file.path(merged_run_dir)
    
    if(collapseNoMis==FALSE){
      list.files(merged_run_path,
                 pattern = glob2rx("*uncollapsed_no-chim-seqtab.rds"),
                 full.names = TRUE,
                 recursive = TRUE) -> seqtab.nochim
    }
    if(collapseNoMis==FALSE){
      list.files(merged_run_path,
                 pattern = glob2rx("*no-chim-seqtab.rds"),
                 full.names = TRUE,
                 recursive = TRUE) -> seqtab.nochim
      
    }
    
    seqtab.nochim <- readRDS(seqtab.nochim)
  }
  if(!is.null(seqtab)){
    
    if(export==TRUE){
      taxa_path <- file.path(taxa_dir)
      
      dir.create(taxa_path, showWarnings = TRUE, recursive = TRUE)
    }
    seqtab.nochim <- seqtab
  }
  ## ------------------------------------------------------------------------
  
  dbname <- str_extract(basename(db), "[^.]+")
  
  set.seed(seed_value) #random  generator necessary for reproducibility
  
  
  # seqtab.nochim <- seqtab.nochim[,1:10] # for testing purpose only
  ## ------------------------------------------------------------------------
  cat(paste0('\n# You have decided to use: ', dbname), ' database \n')
  
  ## ------------------------------------------------------------------------
  
  taxa <- assignTaxonomy(seqtab.nochim, db,
                         outputBootstraps = TRUE,
                         multithread = nthreads,
                         verbose = TRUE,
                         minBoot = threshold,
                         tryRC = as.logical(tryRC),
                         taxLevels = c("Kingdom", "Phylum", "Class",
                                       "Order", "Family", "Genus", "Species"))
  
  
  if(file.exists(db_species))
  {
    boot_taxa <- taxa$boot
    
    taxa_Species <- addSpecies(taxa$tax, db_species,
                               verbose = TRUE,
                               allowMultiple = TRUE,
                               tryRC = as.logical(tryRC))
    
    # Create a merged table with counts and tax
    taxa_full <- left_join(as_tibble(taxa_Species, rownames = 'ASV') %>% 
                             unite("Species",Species:tail(names(.), 1), na.rm = TRUE, sep = "|"), # because / sometimes already Pseudomonas_koreensis(AF468452)/koreensis
                           (as_tibble(taxa$boot, rownames = 'ASV')),
                           by = 'ASV', suffix = c("", "_Boot"))
    
    merged_table <- as_tibble(t(seqtab.nochim), rownames = 'ASV') %>%
      left_join(taxa_full, by = 'ASV') %>%
      mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
      select(ASV_id, everything())
    
    if (export == TRUE){
      
      write_tsv(x = merged_table,
                file = paste0(taxa_path,"/", dbname,"_table.tsv"))
      
      # saveRDS(as_tibble(boot_taxa, rownames = 'ASV'), 
      #         paste0(output,"/", name,"_", dbname,"_boot.rds"))
      
      # saveRDS(as_tibble(taxa_Species, rownames = 'ASV'), 
      #         paste0(output,"/", name,"_", dbname,"_assignation.rds"))
      
      saveRDS(list(as_tibble(taxa_Species, rownames = 'ASV'),
                   as_tibble(boot_taxa, rownames = 'ASV')),
              paste0(taxa_path,"/", dbname,"_assignation.rds"))
      
      cat(paste0('# The obtained taxonomy file can be found in "', paste0(taxa_path,"/", dbname,"_assignation.rds"), '"\n'))
      # cat(paste0('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',paste0(output,"/", name,"_", dbname,"*"),'." and go further ... \n\n'))
    }
  }
  
  if(!file.exists(db_species))
  {
    if(str_extract(basename(db), "[^.]+") == "hitdb_v1")
    {  
      taxa$tax %>% mutate(Kingdom = "Bacteria")
      colnames(taxa$tax) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      
      boot_taxa %>% mutate(Kingdom_Boot = 100)
    }
    # Create a merged table with counts and tax
    taxaid <- left_join(as_tibble(taxa$tax, rownames = 'ASV'),(as_tibble(taxa$boot, rownames = 'ASV')),by = 'ASV', suffix = c("", "_Boot"))
    
    merged_table <- as_tibble(t(seqtab.nochim), rownames = 'ASV') %>%
      left_join(taxaid, by = 'ASV') %>%
      mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
      select(ASV_id, everything())
    
    if (export == TRUE){
      
      write_tsv(x = merged_table,
                file = paste0(taxa_path,"/", dbname,"_table.tsv"))
      
      saveRDS(as_tibble(taxa, rownames = 'ASV'), 
              paste0(taxa_path,"/", dbname,"_assignation.rds"))
      
      cat(paste0('# The obtained taxonomy file can be found in "', paste0(taxa_path,"/", dbname,"_assignation.rds"), '"\n'))
    } # cat(paste0('# You have to copy them to your local computer using "scp [your.user.id]@euler.ethz.ch:',paste0(output,"/", name,"_", dbname,"*"),'." and go further ... \n\n'))
  }
  
  
  
  # merge_phyloseq()
  return("merged_table" = merged_table)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'

run_DECIPHER_phangorn_phylogeny <- function(raw_files_path,
                                            nthreads = 6,
                                            output = "dada2",
                                            phylo_dir = "05_phylo",
                                            merged_run_dir = "03_dada2_merged_runs_chimera_removed",
                                            collapseNoMis = TRUE)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(DECIPHER); require(phangorn)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using DECIPHER version ", packageVersion('DECIPHER'),'\n\n'))
  cat(paste0('\n##',"You are using phangorn version ", packageVersion('phangorn'),'\n\n'))
  
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  setwd(raw_files_path)
  setwd("./..")
  
  phylo_path <- file.path(output, 
                          phylo_dir)
  
  dir.create(phylo_path, showWarnings = TRUE, recursive = TRUE)
  
  cat(paste0('\n# output dir :  ',phylo_path,'\n'))
  
  ## ------------------------------------------------------------------------
  
  merged_run_path <- file.path(output, 
                               merged_run_dir)
  
  if(collapseNoMis==FALSE){
    list.files(merged_run_path,
               pattern = glob2rx("*uncollapsed_no-chim-seqtab.rds"),
               full.names = TRUE,
               recursive = TRUE) %>% readRDS() -> seqtab.nochim
  }else{
    list.files(merged_run_path,
               pattern = glob2rx("*collapse_no_mismatch_no-chim-seqtab.rds"),
               full.names = TRUE,
               recursive = TRUE) %>% readRDS() -> seqtab.nochim
  }
  
  ## ------------------------------------------------------------------------
  # if(method=="R"){
  
  sequences <-  Biostrings::DNAStringSet(getSequences(seqtab.nochim))
  names(sequences) <- sequences  # this propagates to the tip labels of the tree
  
  
  detach('package:phangorn', unload = TRUE)
  detach('package:DECIPHER', unload = TRUE)
  
  # }
  ## ------------------------------------------------------------------------
  # if(method=="MAFFT_Fastree"){
  #   # TODO
  # }
  # physeq@refseq = Biostrings::DNAStringSet(taxa_names(physeq))# https://github.com/benjjneb/dada2/issues/613
  # 
  # taxa_names(physeq)  <- paste0("ASV", str_pad(seq(ntaxa(physeq)), 
  #                       nchar(ntaxa(physeq)), 
  #                       pad = "0"))
  
  ## ------------------------------------------------------------------------
  saveRDS(fitGTR$tree, paste0(phylo_path,"/unrooted_tree.rds"))
  saveRDS(phangorn::midpoint(fitGTR$tree), paste0(phylo_path,"/rooted_tree.rds"))
  
  return(list("unrooted_tree" = fitGTR$tree,
              "rooted_tree" = phangorn::midpoint(fitGTR$tree))) #https://github.com/joey711/phyloseq/issues/936
  
  
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note Need to arrange when Species column is not present...also test with other possible databases/ taxonomic assignments methods.
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'

run_merge_phyloseq <- function(merged_table = NULL,
                               track = NULL,
                               metadata,
                               sample_id = "sample_name",
                               taxa_dir = "dada2/04_dada2_taxonomy",
                               merged_run_dir = "dada2/03_dada2_merged_runs_chimera_removed",
                               collapseNoMis = FALSE,
                               export = "dada2/phyloseq.RDS")
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(phyloseq)
  cat(paste0('\n##',"You are using DADA2 version ", packageVersion('dada2'),'\n'))
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))
  
  cat('################################\n\n')
  
  ## ------------------------------------------------------------------------
  if(is.null(merged_table)){
    
    list.files(file.path(merged_run_dir),
               pattern = glob2rx("track_analysis.tsv"),
               full.names = TRUE,
               recursive = TRUE) %>% read_tsv() -> track
    
    list.files(file.path(taxa_dir),
               pattern = glob2rx("*.tsv"),
               full.names = TRUE,
               recursive = TRUE) %>% read_tsv() -> tmp
  }else{
    tmp <- merged_table
  }
  
  
  tmp %>%
    column_to_rownames("ASV") %>%
    select_if(is.numeric) %>%
    select(-contains("_boot")) %>%
    as.matrix() -> table2
  
  tmp %>%
    select(Kingdom:Species, ASV) %>%
    replace(is.na(.), "unknown") %>%
    column_to_rownames("ASV") %>%
    as.matrix() -> tax_table
  # fo increase genericity need to ignore case for taxonomy and ignore if rank is not present.
  
  phyloseq(tax_table(tax_table),
           otu_table(table2,
                     taxa_are_rows = TRUE))  -> physeq
  
  # prune_samples(sample_sums(physeq) > 0, physeq) -> physeq
  
  ## ------------------------------------------------------------------------
  ## add ASV as refseq part of the phyloseq object
  
  physeq@refseq = Biostrings::DNAStringSet(taxa_names(physeq)) # https://github.com/benjjneb/dada2/issues/613
  
  ## ------------------------------------------------------------------------
  
  taxa_names(physeq) <- paste0("ASV", str_pad(seq(ntaxa(physeq)),
                                              nchar(ntaxa(physeq)),
                                              pad = "0"))
  
  # ifelse(collapseNoMis == TRUE, track_file = "*track_analysis.tsv", track_file = "uncollapsed_track_analysis.tsv")
  
  if(file.exists(metadata))
  {
    full_join(track,
              readxl::read_xlsx(metadata) %>%
                column_to_rownames(sample_id) %>%
                rownames_to_column("sample")) %>%
      column_to_rownames("sample") -> meta
    
    physeq <- merge_phyloseq(physeq,
                             meta %>% sample_data())
  }
  
  if(!file.exists(metadata))
  {
    
    physeq <- merge_phyloseq(physeq,
                             track %>% column_to_rownames("sample") %>% sample_data())
    
  }
  
  if(export != FALSE){
    saveRDS(physeq, 
            export)
  }
  
  return(physeq)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'


add_phylogeny_to_phyloseq <- function(phyloseq_path,
                                      method = "R",
                                      nthreads = 6,
                                      export = FALSE,
                                      output_phyloseq = "phyloseq_phylo.RDS"){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(DECIPHER); require(phangorn); require(phyloseq)
  
  ## ------------------------------------------------------------------------
  if(is.character(phyloseq_path)){
    phyloseq_path %>%
      readRDS() -> physeq
  }else{
    phyloseq_path -> physeq
  }
  
  
  ## ------------------------------------------------------------------------
  if(method=="R"){
    
    
    
    detach('package:phangorn', unload = TRUE)
    detach('package:DECIPHER', unload = TRUE)
    
  }
  ## ------------------------------------------------------------------------
  if(method=="MAFFT_Fastree"){
    # TODO
  }
  # physeq@refseq = Biostrings::DNAStringSet(taxa_names(physeq))# https://github.com/benjjneb/dada2/issues/613
  # 
  # taxa_names(physeq)  <- paste0("ASV", str_pad(seq(ntaxa(physeq)), 
  #                       nchar(ntaxa(physeq)), 
  #                       pad = "0"))
  
  ## ------------------------------------------------------------------------
  # physeq@phy_tree <- 
  
  physeq <- merge_phyloseq(physeq,
                           phangorn::midpoint(fitGTR$tree) %>% phyloseq::phy_tree())
  
  if(export != FALSE)
  {
    dir.create(export, recursive = TRUE)
    
    physeq %>%
      saveRDS(file = file.path(export,
                               output_phyloseq))
  }
  
  
  return(physeq)
  
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note Would be worth checking again confidence score vs taxonomy.
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'

phyloseq_DECIPHER_hitDB <- function(physeq,
                                    db,
                                    threads = 4,
                                    strand="both",
                                    threshold = 10,
                                    verbose = FALSE)
  
{
  require(DECIPHER); require(tidyverse)
  load(db)
  
  #  code below from Dr. Marco Meola:
  
  for(i in c(1:length(ids))) { # go through each row and take out the taxonomy information
    if(i == 1) {
      new_matrix <- data.frame(matrix(nrow=0, ncol=0)) # First, initially an empty dataframe
    }
    # print(i)
    new_line <- as.data.frame(t(as.data.frame(ids[[i]]$taxon))) # create dataframe from the taxon information & transpose
    new_matrix <- plyr::rbind.fill(new_matrix,new_line) # Append it to the bottom of the current dataframe
    # new_list is now an object with the taxonomy information in a simple table.
    if(i == length(ids)) {
      taxid <- as.matrix(new_matrix)
      taxid <- taxid[,-1]
    }
  }
  
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
  colnames(taxid) <- ranks
  rownames(taxid) <- taxa_names(physeq)
  
  taxid %>%
    as.data.frame() %>%
    rownames_to_column('ASV') %>%
    replace(is.na(.), "unknown") -> new_tax
  
  physeq@tax_table = NULL # remove old tax_table
  
  
  merge_phyloseq(physeq,
                 new_tax %>%
                   column_to_rownames('ASV') %>%
                   as.matrix() %>%
                   tax_table()) -> physeq_new
  
  return(physeq_new)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note Would be worth checking again confidence score vs taxonomy.
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'

phyloseq_DECIPHER_tax <- function(physeq, # readRDS("data/processed/physeq_update_11_1_21.RDS") %>% subset_taxa(taxa_sums(physeq) > 100000) -> physeq
                                  threshold = 60, # 60 (very high),  50 (high), PB = 10
                                  db, # db ="~/db/DADA2/SILVA_SSU_r132_March2018.RData" db ="~/db/DADA2/Databases_pbHITdb_pbHITdb_v1.0.0_20180919_IDTAXA.RData" db ="~/db/DADA2/GTDB_r95-mod_August2020.RData" db ="~/db/DADA2/SILVA_SSU_r138_2019.RData" db ="~/db/DADA2/Fungal_LSU_v11_March2018.RData"
                                  nthreads = 6,
                                  tryRC = FALSE,
                                  export = FALSE,
                                  seed_value = 123,
                                  bootlabel = "_Confidence",
                                  tax_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                  return = TRUE)
{
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(DECIPHER); require(phyloseq)
  
  ## ------------------------------------------------------------------------
  
  dbname <- str_extract(basename(db), "[^.]+")
  
  set.seed(seed_value) #random  generator necessary for reproducibility
  
  
  # seqtab.nochim <- seqtab.nochim[,1:10] # for testing purpose only
  ## ------------------------------------------------------------------------
  cat(paste0('\n# Running DECIPHER method against : ', dbname), ' database with threshold value = ', threshold,' \n')
  
  ## ------------------------------------------------------------------------
  if (class(physeq) != "phyloseq"){
    physeq %>% readRDS() -> physeq
  }
  ## ------------------------------------------------------------------------
  sequences <-  Biostrings::DNAStringSet(physeq@refseq)
  
  names(sequences) <- taxa_names(physeq)  # this propagates to the tip labels of the tree
  
  dna <-  Biostrings::DNAStringSet(getSequences(sequences)) # Create a DNAStringSet from the ASVs
  
  load(db)
  
  # trainingSet$ranks %>% unique() # get all tax ranks
  
  ids <- IdTaxa(dna,
                trainingSet,
                strand = if_else(tryRC == TRUE, "both", "top"),
                processors = nthreads,
                verbose = TRUE,
                threshold = threshold)
  
  # if(grepl("pbHITdb", str_extract(basename(db), "[^.]+"), fixed = TRUE) == TRUE)
  # {  
  #   for(i in c(1:length(ids))) { # go through each row and take out the taxonomy information
  #     if(i == 1) {
  #       new_matrix <- data.frame(matrix(nrow=0, ncol=0)) 
  #     }
  #     new_line <- as.data.frame(t(as.data.frame(ids[[i]]$taxon))) 
  #     new_matrix <- plyr::rbind.fill(new_matrix,
  #                                    new_line) 
  #     
  #     if(i == length(ids)) {
  #       taxid <- as.matrix(new_matrix)
  #       taxid <- taxid[,-1]
  #       
  #     }
  #   }
  #   taxid %>%
  #     data.frame() -> taxid
  #   
  #   ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
  #   
  #   colnames(taxid) <- ranks
  #   rownames(taxid) <- physeq %>% 
  #     taxa_names()
  # 
  # }else{
  
  
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
  
  # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
  taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  
  colnames(taxid) <- tax_ranks
  # }  
  taxid %>% 
    replace_na("unknown") -> taxid
  
  physeq@tax_table = NULL
  
  taxid %>% 
    as.matrix() %>%
    tax_table() ->  physeq@tax_table
  
  # Create a merged table with counts and tax and confidence
  
  as(tax_table(physeq), "matrix") %>% 
    data.frame() %>%
    rownames_to_column('ASV') -> tax_table
  
  as(otu_table(physeq), "matrix") %>% 
    as.data.frame() %>%
    rownames_to_column('ASV') -> otu_table
  
  dss2df <- function(dss) data.frame(width=BiocGenerics::width(dss), seq=as.character(dss), names=names(dss))
  
  dss2df(physeq@refseq) %>%
    rownames_to_column('ASV') %>%
    dplyr::select(-names) %>%
    dplyr::rename(ASV_length = width,
                  ASV_sequence = seq)-> refseq_df
  
  tax_score <- sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$confidence[m]
  }) %>%
    t() %>% 
    data.frame() 
  
  colnames(tax_score) <- paste0(colnames(taxid), 
                                bootlabel)
  
  otu_table %>%
    left_join(tax_table, by = 'ASV') %>%
    left_join(tax_score    %>%
                rownames_to_column('ASV')) %>%
    left_join(refseq_df, by = 'ASV') %>%
    dplyr::select(ASV, everything()) -> merged_table
  # 
  # 
  # otu_table %>%
  #   left_join(tax_table, by = 'ASV') %>%
  #   left_join(tax_score    %>%
  #               rownames_to_column('ASV') %>%
  #               mutate_at(
  #                 vars(ends_with(bootlabel)),
  #                 funs(case_when(
  #                   . < threshold ~ "unknown", #       . == 21.0 ~ 5,
  #                   TRUE ~ .
  #                 ))
  #               ))
  # # dplyr::replace(vars(contains("_Confidence")) > thresold)) %>%
  # 
  # tax_score    %>%
  #   rownames_to_column('ASV') %>%
  #   mutate_at(
  #     vars(ends_with(bootlabel)),
  #     funs(case_when(
  #       . == "1+1" | . == "1+2" ~ 1,
  #       . == "1+3" | . == "1+4" ~ 2,
  #       . == "1+5" | . == "1+6" ~ 3)))
  # 
  # dat %>%
  #   mutate(var = case_when(var == 'Candy' ~ 'Candy',
  #                          var == 'Water' ~ 'Water',
  #                          TRUE ~ 'Neither-Water-Nor-Candy'))
  
  
  # mutate(ASV_id = paste0("asv",c(1:nrow(.)))) %>%
  # mutate(threshold = threshold) %>%
  
  if(export != FALSE)
  {
    dir.create(export, recursive = TRUE)
    
    write_tsv(x = merged_table,
              file = paste0(export,"/","DECIPHER_threshold_",threshold,"_",dbname,"merged_tax_table.tsv"))
    cat(paste0("saved object and output to ", export))
    saveRDS(merged_table,
            file = paste0(export,"/","DECIPHER_threshold_",threshold,"_",dbname,".RDS"))
    saveRDS(physeq,
            file = paste0(export,"/","DECIPHER_threshold_",threshold,"_",dbname,"_physeq.RDS"))
  }
  
  if(return == TRUE){
    return(physeq)
  }
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'readRDS("data/processed/physeq_update_11_1_21.RDS") %>% subset_taxa(taxa_sums(physeq) > 100000) -> physeq
#'physeq %>%
#'  phyloseq_dada2_tax(db ="~/db/DADA2/silva_nr99_v138_train_set.fa.gz",
#'                     db_species = "none") -> p_test

phyloseq_dada2_tax <- function(physeq, # readRDS("data/processed/physeq_update_11_1_21.RDS") %>% subset_taxa(taxa_sums(physeq) > 100000) -> physeq
                               threshold = 60,
                               db, # db ="~/db/DADA2/silva_nr99_v138_train_set.fa.gz"
                               db_species, # db_species ="~/db/DADA2/silva_species_assignment_v138.fa.gz"
                               nthreads = 6,
                               tryRC = FALSE,
                               seed_value = 123,
                               bootlabel = "_Boot",
                               tax_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                               export = FALSE,
                               return = TRUE,
                               full_return = TRUE){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(dada2); require(phyloseq)
  
  ## ------------------------------------------------------------------------
  
  dbname <- str_extract(basename(db), "[^.]+")
  
  set.seed(seed_value) #random  generator necessary for reproducibility
  
  
  # seqtab.nochim <- seqtab.nochim[,1:10] # for testing purpose only
  ## ------------------------------------------------------------------------
  ## ------------------------------------------------------------------------
  if (class(physeq) != "phyloseq"){
    physeq %>% readRDS() -> physeq
  }
  ## ------------------------------------------------------------------------
  # Create a merged table with counts taxa and sequences
  
  # as(otu_table(physeq), "matrix") %>% 
  #   data.frame() %>%
  #   rownames_to_column('ASV') -> otu_table
  # 
  dss2df <- function(dss) data.frame(width=BiocGenerics::width(dss), seq=as.character(dss), names=names(dss))
  
  dss2df(physeq@refseq) -> refseq_df
  # 
  # 
  # 
  # 
  # # if(!is.null(tax_table(physeq))){
  # as(tax_table(physeq), "matrix") %>% 
  #   data.frame() %>%
  #   rownames_to_column('ASV') -> tax_table
  # colnames(tax_table)[-1] <- paste0("input_taxonomy_", 
  #                                   colnames(tax_table)[-1])
  # 
  # otu_table %>%
  #   left_join(tax_table, by = 'ASV') %>%
  #   left_join(refseq_df, by = 'ASV') %>%
  #   dplyr::select(ASV, everything()) -> merged_table_pre
  # }else{
  #   
  #   otu_table %>%
  #     left_join(refseq_df, by = 'ASV') %>%
  #     dplyr::select(ASV, everything()) -> merged_table_pre
  # }
  
  ## ------------------------------------------------------------------------
  
  cat(paste0('\n# Running dada2 method against : ', dbname), ' database with threshold value = ', threshold,' \n')
  
  ## ------------------------------------------------------------------------
  sequences <-  Biostrings::DNAStringSet(physeq@refseq)
  
  names(sequences) <- taxa_names(physeq)  # this propagates to the tip labels of the tree
  
  dna <-  Biostrings::DNAStringSet(getSequences(sequences)) # Create a DNAStringSet from the ASVs
  
  if(str_extract(basename(db), "[^.]+") == "hitdb_v1")
  {  
    tax_ranks <-  c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  }
  
  taxa <- assignTaxonomy(dna, 
                         db,
                         outputBootstraps = TRUE,
                         multithread = nthreads,
                         verbose = TRUE,
                         minBoot = threshold,
                         tryRC = as.logical(tryRC),
                         taxLevels = tax_ranks)
  
  if(!file.exists(db_species))
  {
    
    # left_join(merged_table_pre,
    #           left_join(as_tibble(taxa$tax, rownames = 'ASV'),
    #                     (as_tibble(taxa$boot, rownames = 'ASV')),
    #                     by = 'ASV', suffix = c("", "_Boot")) ,
    #           by = c('ASV_sequence' = 'ASV')) -> merged_table
    
    left_join(refseq_df,
              left_join(as_tibble(taxa$tax, rownames = 'ASV'),
                        (as_tibble(taxa$boot, rownames = 'ASV')),
                        by = 'ASV', suffix = c("", "_Boot")) ,
              by = c('ASV_sequence' = 'ASV')) -> merged_table
    
    if(str_extract(basename(db), "[^.]+") == "hitdb_v1")
    {  
      merged_table %>%
        dplyr::mutate(Kingdom = "Bacteria") -> merged_table
    }
  }  
  if(file.exists(db_species))
  {
    taxa_Species <- addSpecies(taxa$tax, 
                               db_species,
                               verbose = TRUE,
                               allowMultiple = TRUE,
                               tryRC = as.logical(tryRC))
  }
  if(file.exists(db_species) || str_extract(basename(db), "[^.]+") == "hitdb_v1")
  {
    # Create a merged table with counts and tax
    taxa_full <- as_tibble(taxa_Species, 
                           rownames = 'ASV')  %>% 
      unite("Species",Species:tail(names(.), 1), na.rm = TRUE, remove = TRUE, sep = "|")
    
    left_join(taxa_full,
              as_tibble(taxa$boot, rownames = 'ASV'),
              by = 'ASV', suffix = c("", "_Boot")) -> merged_table
    
    # if(!is.null(tax_table(physeq))){
    # merged_table %>%
    #   dplyr::select(-starts_with("input_taxonomy_")) -> merged_table
    # }
    merged_table %>%
      dplyr::select_if(is.character) %>%
      # dplyr::select(-ASV) %>%    
      column_to_rownames("ASV") -> tax_table
    
  }else{
    merged_table %>%
      dplyr::select_if(is.character) %>%
      column_to_rownames("ASV") -> tax_table
    
  }
  
  tax_table %>%
    rownames_to_column("ASV_sequence") %>%
    # as.data.frame() %>%
    left_join(refseq_df) %>%
    column_to_rownames('ASV') %>%
    dplyr::select(-ASV_length, -ASV_sequence) -> tmp
  
  physeq@tax_table = NULL
  
  merge_phyloseq(physeq,
                 tmp %>% 
                   as.matrix() %>%
                   tax_table()) -> physeq
  
  if(export != FALSE){
    dir.create(export, recursive = TRUE)
    
    write_tsv(x = merged_table,
              file = paste0(export,"/","dada2_threshold",threshold,"_",dbname,"merged_tax_table.tsv"))
    cat(paste0("saved object and output to ", export))
    
    saveRDS(merged_table,
            file = paste0(export,"/","dada2_threshold",threshold,"_",dbname,".RDS"))
    
    saveRDS(physeq,
            file = paste0(export,"/","dada2_threshold",threshold,"_",dbname,"_physeq.RDS"))
    
  }
  if(return == TRUE)
  {
    if(full_return == TRUE){
      out <- list("physeq" = physeq,
                  "full_table" = merged_table)
    }else{
      out <- physeq
    }
    return(out)
  }
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'

compare_phyloseq_taxonomy <- function(physeq_A,
                                      physeq_B){
  physeq_A %>% 
    phloseq_export_otu_tax() %>% 
    select(-ASV_sequence) %>%
    dplyr::select(-sample_names(physeq_A)) %>%
    left_join(
      physeq_B %>%
        phloseq_export_otu_tax(),
      by = "ASV" ,
      suffix = c("_A", "_B")) %>% 
    select(-contains(c("Kingdom", "Phylum", "Class", "Order","length"))) %>%
    select(-ASV_sequence) %>%
    select(contains(c("_A", "_B")), everything()) -> new_old_tax_df
  
  return(new_old_tax_df)
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
# readRDS("data/processed/physeq_update_11_1_21.RDS") %>% 
#   phyloseq_DECIPHER_cluster_ASV(threshold = 99) -> test_99
#'
# test_99$cluster_table %>%
#  arrange(-cluster_size) %>% 
# select(ASV, contains("cluster_"), c("Phylum", "Class", "Family", "Genus", "Species"))
#'
#'


phyloseq_DECIPHER_cluster_ASV <- function(physeq, # readRDS("data/processed/physeq_update_11_1_21.RDS") %>% subset_taxa(taxa_sums(physeq) > 1000) -> physeq
                                          threshold = 97,
                                          nthreads = 6,
                                          # method = "complete",
                                          # showPlot = FALSE,
                                          export = FALSE,
                                          return = TRUE,
                                          full_return = TRUE){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(speedyseq); require(DECIPHER)
  
  ## ------------------------------------------------------------------------
  if (class(physeq) != "phyloseq"){
    physeq %>% readRDS() -> physeq
  }
  ## ------------------------------------------------------------------------
  
  # Create a merged table with counts taxa and sequences
  as(tax_table(physeq), "matrix") %>% 
    data.frame() %>%
    rownames_to_column('ASV') -> tax_table
  
  
  as(otu_table(physeq), "matrix") %>% 
    as.data.frame() %>%
    rownames_to_column('ASV') -> otu_table
  
  dss2df <- function(dss) data.frame(width=BiocGenerics::width(dss), seq=as.character(dss), names=names(dss))
  
  dss2df(physeq@refseq) %>%
    rownames_to_column('ASV') %>%
    dplyr::select(-names) %>%
    dplyr::rename(ASV_length = width,
                  ASV_sequence = seq)-> refseq_df
  
  
  otu_table %>%
    left_join(tax_table, by = 'ASV') %>%
    left_join(refseq_df, by = 'ASV') %>%
    dplyr::select(ASV, everything()) -> merged_table_pre
  
  
  ## ------------------------------------------------------------------------
  dna <- refseq(physeq)
  
  # aln <- DECIPHER::AlignSeqs(dna, 
  #                            processors = nthreads)
  # 
  # d <- DECIPHER::DistanceMatrix(aln, 
  #                               processors = nthreads)
  
  ## ------------------------------------------------------------------------
  
  clusters <- DECIPHER::IdClusters(
    # d, 
    dna,
    # showPlot = showPlot,
    # method = method,
    cutoff = (100-threshold) / 100, # corresponds to 97% OTUs
    processors = nthreads,
    verbose = FALSE)
  
  ## ------------------------------------------------------------------------
  
  clusters %>%
    rownames_to_column('ASV') %>%
    group_by(cluster) %>%
    add_count() %>%
    dplyr::rename(!!paste0("cluster_", (100-threshold) / 100) := cluster,
                  cluster_size = n) %>%
    full_join(merged_table_pre) -> full_table
  
  
  ## ------------------------------------------------------------------------
  
  physeq_clustered <- speedyseq::merge_taxa_vec(
    physeq,
    group = clusters$cluster,
    tax_adjust = 2)
  
  ## ------------------------------------------------------------------------
  
  if(export != FALSE){
    
    dir.create(export, recursive = TRUE) 
    
    write_tsv(x = full_table,
              file = paste0(export,"/","physeq_ASV_clust_thresold_",(100-threshold) / 100,"_table.tsv"))
    
    cat(paste0("saved objects and outputs to ", export))
    
    saveRDS(physeq_clustered,
            file = paste0(export,"/","physeq_ASV_clust_thresold_",(100-threshold) / 100,"_physeq.RDS"))
    
  }
  ## ------------------------------------------------------------------------
  
  if(return == TRUE)
  {
    if(full_return == TRUE){
      out <- list("physeq_clustered" = physeq_clustered,
                  "cluster_table" = full_table)
    }else{
      out <- physeq_clustered
    }
    return(out)
  }
}

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'phyloseq_vsearch_lulu_cluster_ASV(readRDS("data/processed/physeq_update_11_1_21.RDS"),
#'                               vsearch = "/Users/fconstan/miniconda3/envs/metabarcodingRpipeline/bin/vsearch",
#'                                 dir = ("~/"), int_rm = TRUE) -> out

phyloseq_vsearch_lulu_cluster_ASV <- function(physeq, # readRDS("data/processed/physeq_update_11_1_21.RDS") -> physeq
                                              vsearch = "vsearch",
                                              dir = ("~/"),
                                              fasta_file = "pre_lulu.fasta",
                                              match_list_file = "pre_lulu_vsearch.list.txt",
                                              threshold = 80,
                                              nthreads = 6,
                                              int_rm = TRUE,
                                              seed_value = 123,
                                              minimum_ratio_type = "min",
                                              minimum_match = 84,
                                              minimum_relative_cooccurence = 0.95,
                                              minimum_ratio = 1, # only if minimum_ratio_type = "avg"
                                              export = FALSE,
                                              return = TRUE,
                                              full_return = TRUE){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq); require(lulu)
  
  
  dss2df <- function(dss) data.frame(width=BiocGenerics::width(dss), seq=as.character(dss), names=names(dss))
  
  ## ------------------------------------------------------------------------
  if (class(physeq) != "phyloseq"){
    physeq %>% readRDS() -> physeq
  }
  
  ## ------------------------------------------------------------------------
  dir.create(dir, recursive = TRUE)
  
  paste0(dir, "/", fasta_file) -> fasta_path
  paste0(dir, "/", match_list_file) -> match_list_file_path
  
  ## ------------------------------------------------------------------------
  
  refseq(physeq) %>%
    dss2df() %>%
    dplyr::rename(sequence = seq) %>%
    dplyr::mutate(abundance = taxa_sums(physeq)) %>%
    dplyr::select(-names) %>%
    dada2::uniquesToFasta(fasta_path,
                          ids=taxa_names(physeq))
  
  as(otu_table(physeq), "matrix") %>%
    as.data.frame()  -> otu_tab
  ## ------------------------------------------------------------------------
  system2(vsearch, args = c("--version"))
  
  ## run vsearch
  system2(vsearch, args = c("--usearch_global ", fasta_path, 
                            "--db", fasta_path,
                            "--self --id", threshold/100,
                            "--iddef 1 --userout ", match_list_file_path, 
                            "--threads ", nthreads,
                            "-userfields query+target+id --maxaccepts 0 --query_cov 0.9 --maxhits 10")
  )
  # vsearch --usearch_global OTU_sequences.fasta --db OTU_sequences.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
  
  ## ------------------------------------------------------------------------
  match_list_file_path %>%
    read_tsv(col_names = c("ASVa", 
                           "ASVb", 
                           "id")) -> match_list_file
  
  if(int_rm == TRUE){
    file.remove(fasta_path, match_list_file_path)
    
  }
  ## ------------------------------------------------------------------------
  
  if (minimum_ratio_type == "avg")
  {
    curated_result <- lulu(otu_tab, 
                           match_list_file, 
                           minimum_ratio_type = minimum_ratio_type, 
                           minimum_match = minimum_match, 
                           minimum_ratio = minimum_ratio, 
                           minimum_relative_cooccurence = minimum_relative_cooccurence)
  }
  if (minimum_ratio_type == "min")
  {
    curated_result <- lulu(otu_tab, 
                           match_list_file, 
                           minimum_ratio_type = minimum_ratio_type, 
                           minimum_match = minimum_match, 
                           minimum_relative_cooccurence = minimum_relative_cooccurence)
  }
  
  ## ------------------------------------------------------------------------
  
  # Function lulu returns a list of results based on the input OTU table and match list.
  # 
  # curated_table - a curated OTU table with daughters merged with their matching parents.
  # 
  # curated_count - number of curated (parent) OTUs.
  # 
  # curated_otus - ids of the OTUs that were accepted as valid OTUs.
  # 
  # discarded_count - number of discarded (merged with parent) OTUs.
  # 
  # discarded_otus - ids of the OTUs that were identified as errors (daughters) and merged with respective parents.
  # 
  # runtime - time used by the script.
  # 
  # minimum_match - the id threshold (minimum match % between parent and daughter) for evaluating co-occurence (set by user).
  # 
  # minimum_relative_cooccurence - minimum ratio of daughter-occurences explained by co-occurence with parent (set by user).
  # 
  # otu_map - information of which daughters were mapped to which parents.
  # 
  # original_table - original OTU table.
  
  
  ## ------------------------------------------------------------------------
  colnames(curated_result$curated_table) <- sample_names(physeq)
  colnames(curated_result$original_table) <- sample_names(physeq)
  
  curated_result$curated_table %>%
    rownames_to_column('ASV') %>%
    pull("ASV") -> curated_asv
  
  
  physeq_curated <- prune_taxa(curated_asv, physeq)
  
  curated_result$curated_table %>% 
    #   dplyr::mutate_all(funs(str_replace(., " ", "")))
    # # stringr::str_replace_all()
    as.matrix() -> new_otu_tab
  
  otu_table(physeq_curated) <- otu_table(new_otu_tab, taxa_are_rows = TRUE)
  
  ## ------------------------------------------------------------------------
  
  if(export != FALSE){
    
    # write_tsv(x = curated_result$curated_table,
    #           file = paste0(dir,"/","lulu_curated_table.tsv"))
    
    openxlsx::write.xlsx(curated_result, 
                         file.path(dir,"lulu_curated_table.xlsx"))
    
    write_tsv(x = curated_result$curated_table,
              file = file.path(dir,"lulu_curated_table.tsv"))
    
    cat(paste0("saved object and output to ", dir))
    
    saveRDS(physeq_curated,
            file = file.path(dir,"lulu_curated_physeq.RDS"))
    
  }
  ## ------------------------------------------------------------------------
  
  if(return == TRUE)
  {
    if(full_return == TRUE){
      out <- list("physeq_curated" = physeq_curated,
                  "curated_result" = curated_result)
    }else{
      out <- physeq_curated
    }
    return(out)
  }
}

# comments below
co_au=TRUE

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note Generetate MIMOSA2 compatible file- <https://borenstein-lab.github.io/MIMOSA2shiny/results.html#>
#' @note https://www.youtube.com/watch?v=xZ7yc-GKcSk
#' @note https://ycl6.github.io/16S-Demo/4_picrust2_tutorial.html#Perform_statistical_analysis
#' @note https://github.com/picrust/picrust2/wiki/Frequently-Asked-Questions#how-can-i-determine-kegg-pathway-abundances-from-the-predicted-ko-abundances
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'phyloseq_picrust2(physeq = readRDS("/Users/fconstan/Documents/GitHub/amchick/data/processed/physeq_update_11_1_21.RDS"),
#'                  input_dir = ("/Users/fconstan/Desktop/picrust_test/AP-picrust/"),
#'                  output_dir = ("/Users/fconstan/Desktop/picrust_test/AP-picrust/results"),# should not exist
#'                  in_traits = c("COG,EC,KO,PFAM,TIGRFAM"),
#'                  min_reads = 2,
#'                  min_samples = 3)
#' cd /Users/fconstan/Desktop/picrust_test/chemerin_16S/picrust2_out_pipeline/EC_metagenome_out 
#' 
#' zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz
#' zless -S EC_metagenome_out/pred_metagenome_contrib.tsv.gz 

phyloseq_picrust2 <- function(physeq = NULL, # readRDS("data/processed/physeq_update_11_1_21.RDS") -> physeq
                              picrust2 = "picrust2_pipeline.py",
                              # conda_env = "picrust2",
                              input_dir = ("~/"),
                              output_dir = ("~/picrust2/"),# should not exist
                              fasta_file = "ASV.fna",
                              count_table = "ASV_table.tsv",
                              in_traits = c("EC,KO"), #  #default: EC,KO c("COG,EC,KO,PFAM,TIGRFAM")
                              min_reads = 1,
                              min_samples = 1,
                              ref_dir = FALSE, #"/Users/fconstan/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/prokaryotic/pro_ref",
                              pathway_map = FALSE, #"/Users/fconstan/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_pro.txt",
                              regroup_map = FALSE, #"/Users/fconstan/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/pathway_mapfiles/ec_level4_to_metacyc_rxn.tsv",
                              nthreads = 6,
                              m = "mp",
                              seed = 123456,
                              no_gap_fill = FALSE,
                              add_description = FALSE,
                              load_picrust2_data = FALSE, #could generate large files...
                              return = FALSE, # true / false or path of object
                              int_rm = FALSE){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq)#; require(reticulate)
  
  
  dss2df <- function(dss) data.frame(width = BiocGenerics::width(dss), seq=as.character(dss), names=names(dss))
  
  ## ------------------------------------------------------------------------
  if (class(physeq) != "phyloseq"){
    physeq %>% readRDS() -> physeq
  }
  
  ## ------------------------------------------------------------------------
  
  # reticulate::use_condaenv(condaenv = conda_env,
  #                          conda = "auto",
  #                          required = FALSE)
  # system2("conda activate",
  #         args = c(conda_env)
  # )
  ## ------------------------------------------------------------------------
  
  dir.create(input_dir, recursive = TRUE)
  
  paste0(input_dir, "/", fasta_file) -> fasta_path
  paste0(input_dir, "/", count_table) -> count_table_path
  
  ## ------------------------------------------------------------------------
  if(!is.null(physeq)){
    
    refseq(physeq) %>%
      dss2df() %>%
      dplyr::rename(sequence = seq) %>%
      dplyr::mutate(abundance = taxa_sums(physeq)) %>%
      dplyr::select(-names) %>%
      dada2::uniquesToFasta(fasta_path,
                            ids=taxa_names(physeq))
    
    as(otu_table(physeq), "matrix") %>%
      as.data.frame() %>%
      rownames_to_column('ASV') %>%
      # column_to_rownames('ASV') %>%
      write_tsv(count_table_path)
  }
  ## ------------------------------------------------------------------------
  
  ## run picrust2_pipeline.py
  system2(picrust2, 
          args = c("-s ", fasta_path, 
                   "-i ", count_table_path,
                   "-o ", output_dir,
                   "--processes ", nthreads,
                   "--stratified --per_sequence_contrib --verbose ",
                   "--in_traits ", in_traits,
                   "--min_reads", min_reads, "--min_samples", min_samples,
                   ifelse(pathway_map != FALSE, print("--pathway_map " , pathway_map), ""),
                   ifelse(ref_dir != FALSE, print("--ref_dir ", ref_dir), ""),
                   ifelse(regroup_map != FALSE, print("--regroup_map ", regroup_map), ""),
                   "--hsp_method ", m, ifelse(no_gap_fill == TRUE, print("--no_gap_fill"), ""),
                   ifelse(int_rm == TRUE, print("--remove_intermediate"), ""))
  )
  # https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.3.0-beta)#add-functional-descriptions
  # picrust2_pipeline.py -s study_seqs.fna -i seqabun.biom -o picrust2_out --processes 10 --stratified --per_sequence_contrib 
  
  ## ------------------------------------------------------------------------
  # add descriptions <https://github.com/picrust/picrust2/wiki/Add-descriptions>
  
  if (add_description == TRUE){
    
    unlist(lapply(strsplit(in_traits, ","), as.character)) -> in_traits
    
    system2("add_descriptions.py", 
            args = c("-i ", paste0(output_dir, "/", "pathways_out", "/", "path_abun_unstrat.tsv.gz"), 
                     "-m ", "METACYC",
                     "-o ", paste0(output_dir, "/", "pathways_out", "/", "path_abun_unstrat.tsv.gz"))) #path_abun_unstrat_descrip.tsv.gz
    
    # if (c("EC") %in% in_traits){
    #   
    #   system2("add_descriptions.py", 
    #           args = c("-i ", paste0(output_dir, "/", "EC_metagenome_out", "/", "pred_metagenome_unstrat.tsv.gz"), 
    #                    "-m ", "EC",
    #                    "-o ", paste0(output_dir, "/", "EC_metagenome_out", "/", "pred_metagenome_unstrat.tsv.gz")))
    # }
    # if (c("KO") %in% in_traits){
    #   
    #   system2("add_descriptions.py", 
    #           args = c("-i ", paste0(output_dir, "/", "KO_metagenome_out", "/", "pred_metagenome_unstrat.tsv.gz"), 
    #                    "-m ", "KO",
    #                    "-o ", paste0(output_dir, "/", "KO_metagenome_out", "/", "pred_metagenome_unstrat.tsv.gz"))) #pred_metagenome_unstrat_descrip.tsv.gz
    # }
    
    for (trait in in_traits){
      system2("add_descriptions.py",
              args = c("-i ", paste0(output_dir,
                                     "/",
                                     trait,
                                     "_metagenome_out",
                                     "/",
                                     "pred_metagenome_unstrat.tsv.gz"),
                       "-m ", trait,
                       "-o ", paste0(output_dir,
                                     "/",
                                     trait,
                                     "_metagenome_out",
                                     "/",
                                     "pred_metagenome_unstrat.tsv.gz")))
      
    }
    
  }
  
  # 
  # add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
  # -o EC_metagenome_out/pred_metagenome_unstrat.tsv.gz
  # 
  # add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
  # -o KO_metagenome_out/pred_metagenome_unstrat.tsv.gz
  # 
  # add_descriptions.py -i PFAM_metagenome_out/pred_metagenome_unstrat.tsv.gz -m PFAM \
  # -o PFAM_metagenome_out/pred_metagenome_unstrat.tsv.gz
  # 
  # add_descriptions.py -i TIGRFAM_metagenome_out/pred_metagenome_unstrat.tsv.gz -m TIGRFAM \
  # -o TIGRFAM_metagenome_out/pred_metagenome_unstrat.tsv.gz
  # 
  # add_descriptions.py -i COG_metagenome_out/pred_metagenome_unstrat.tsv.gz -m COG \
  # -o COG_metagenome_out/pred_metagenome_unstrat.tsv.gz
  # 
  # add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  # -o pathways_out/path_abun_unstrat.tsv.gz
  # 
  
  ## ------------------------------------------------------------------------
  
  if (load_picrust2_data == TRUE){
    
    #marker_predicted_and_nsti.tsv.gz : 
    # The first column is the name of ASV, followed by the predicted number of 16S copies per ASVs, f
    # ollowed finally by the NSTI value per ASV. 
    # When testing several datasets we found that ASVs with a NSTI score above 2 are usually noise. 
    # It can be useful to take a look at the distribution of NSTI values for your ASVs to determine how well-characterized your community is overall and whether there are any outliers.
    
    #pred_metagenome_unstrat.tsv.gz
    #pred_metagenome_contrib.tsv.gz
    #path_abun_unstrat.tsv.gz
    #path_abun_unstrat_per_seq
    #path_abun_contrib.tsv.gz
    
    pattern = c("marker_predicted_and_nsti|pred_metagenome_unstrat|pred_metagenome_contrib|path_abun_unstrat|path_abun_unstrat_per_seq|path_abun_contrib")
    
    list.files(file.path(output_dir), 
               full.names = TRUE,
               recursive = TRUE,
               pattern = pattern) %>%
      sort() -> list_files
    
    list.files(file.path(output_dir), 
               full.names = FALSE,
               recursive = TRUE,
               pattern = pattern) %>%
      sort() %>% 
      str_replace("/", "_") %>%
      str_remove(".tsv.gz") -> list_files_names
    
    
    files_list <- lapply(list_files,function(x) {
      readr::read_tsv(file = file.path(x))
    })
    
    names(files_list) <-     
      list_files_names 
    
  }
  
  
  if (return == TRUE){
    return(files_list)
  }
  
  if (is.character(return)){
    
    dir.create(return, recursive = TRUE)
    
    files_list %>%
      saveRDS(paste0(return,"/", "picrust2_R.rds"))
  }
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'

phloseq_export_otu_tax <- function(physeq = NULL){
  
  require(phyloseq); require(tidyverse)
  
  as(tax_table(physeq), "matrix") %>% 
    as.data.frame() %>%
    rownames_to_column('ASV') -> tax_table
  
  as(otu_table(physeq), "matrix") %>% 
    as.data.frame() %>%
    rownames_to_column('ASV') -> otu_table
  
  # colnames(otu_table) <- c('ASV', sample_names(physeq))
  
  dss2df <- function(dss) data.frame(width=BiocGenerics::width(dss), seq=as.character(dss), names=names(dss))
  
  dss2df(physeq@refseq) %>%
    rownames_to_column('ASV') %>%
    dplyr::select(-names) %>%
    dplyr::rename(ASV_length = width,
                  ASV_sequence = seq)-> refseq_df
  
  otu_table %>%
    left_join(refseq_df, by = 'ASV') %>%
    left_join(tax_table, by = c("ASV" = "ASV")) %>%
    dplyr::select(ASV, everything()) -> merged_table
  
  return(merged_table)
}

# comments below:
# from previous repo
master="https://raw.githubusercontent.com/fconstancias/metabaRpipe-source/master/Rscripts/functions.R"


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'
#'
#'
#'
#'
#'

run_dada2_pipe <- function(raw_files_path,
                           raw_file_pattern = c("*_R1_*.gz","*_R2_*.gz"),
                           atropos_bin = "atropos",
                           cut_dir = "dada2/00_atropos_primer_removed",
                           cut_file_pattern = c("_primersout_R1_.fastq.gz","_primersout_R2_.fastq.gz"),
                           filt_dir = "dada2/02_dada2_filtered_denoised_merged",
                           merged_run_dir = "dada2/03_dada2_merged_runs_chimera_removed",
                           taxa_dir = "dada2/04_dada2_taxonomy",
                           sep = "[^_]+",
                           V = "V4-2PCR",
                           rm_primers = TRUE,
                           PRIMER_F,
                           PRIMER_R,
                           tax_threshold = 60,
                           nbases = 20000000,
                           pool = "pseudo",
                           priors = FALSE,
                           trim_length = c(240,400),
                           trunclen = c(260,250),
                           truncQ = 6,
                           maxee = c(3,4),
                           minLen = 100,
                           minover = 15,
                           chimera_method = "consensus",
                           collapseNoMis = FALSE,
                           tryRC = TRUE,
                           metadata = "none",
                           db = "~/db/DADA2/silva_nr_v138_train_set.fa.gz",
                           db_species = "~/db/DADA2/silva_species_assignment_v138.fa.gz",
                           merge_phyloseq_export = "dada2/phyloseq.RDS",
                           run_phylo = FALSE,
                           output_phyloseq_phylo = "dada2/phyloseq_phylo.RDS",
                           save_out = "dada2/dada2_pipe_out.RDS",
                           export = TRUE,
                           remove_input_fastq = TRUE,
                           SLOTS = 6,
                           seed_value = 123,
                           return = FALSE){
  if(V == "V4-2PCR") {
    
    PRIMER_F = "GTGCCAGCMGCCGCGGTAA"
    PRIMER_R = "GGACTACHVGGGTWTCTAAT" 
    trim_length = c(220,280)
    trunclen =  c(170,160)
    maxee = c(3,4)
    minLen = 120
    minover = 40
    nbases = 100000000000
  } 
  if(V == "V4-1PCR") {
    
    PRIMER_F = "GTGCCAGCMGCCGCGGTAA"
    PRIMER_R = "GGACTACHVGGGTWTCTAAT" 
    trim_length = c(220,280)
    trunclen =  c(170,160)
    maxee = c(3,4)
    minLen = 120
    minover = 40
    nbases = 100000000000
    rm_primers = FALSE
  } 
  if(V == "V4-NIH") {
    
    PRIMER_F = "GTGYCAGCMGCCGCGGTAA"
    PRIMER_R = "GGACTACNVGGGTWTCTAAT" 
    trim_length = c(200,280)
    truncQ = Inf
    trunclen =  c(170,160)
    maxee = c(3,4)
    minLen = 120
    minover = 40
    nbases = 100000000000
    rm_primers = FALSE
    collapseNoMis = TRUE
    
  } 
  if(V == "test") {
    
    PRIMER_F = "GTGYCAGCMGCCGCGGTAA"
    PRIMER_R = "GGACTACNVGGGTWTCTAAT" 
    rm_primers = FALSE
    trim_length = c(200,280)
    truncQ = Inf
    trunclen =  c(170,160)
    maxee = c(3,4)
    minLen = 120
    minover = 40
    nbases = 10
  } 
  if(V == "V4-Addition-PRO") {
    
    PRIMER_F = "GTGCCAGCMGCCGCGGTAA"
    PRIMER_R = "GGACTACHVHHHTWTCTAAT" 
    trim_length = c(220,280)
    trunclen =  c(170,160)
    maxee = c(3,4)
    minLen = 120
    minover = 40
  } 
  if(V == "V3V4"){
    
    PRIMER_F = "CCTAYGGGRBGCASCAG"
    PRIMER_R = "GGACTACNNGGGTATCTAAT"
    trim_length = c(240,600)
    trunclen =  c(260,250)
    maxee = c(4,5)
    minLen = 160
    minover = 10
    
  }
  if(V == "V3V4_2x250"){
    
    PRIMER_F = "CCTAYGGGRBGCASCAG"
    PRIMER_R = "GGACTACNNGGGTATCTAAT"
    trunclen  = c(230,227)
    trim_length =  c(220,480)
    maxee = c(3,4)
    minLen = 160
    minover = 10
    
  }
  if(V == "ITS2"){
    
    PRIMER_F = "GTGAATCATCGAATCTTTGAA"
    PRIMER_R = "TCCTCCGCTTATTGATATGC"
    trim_length = c(200,500)
    trunclen =  c(230,220)
    maxee = c(4,5)
    minLen = 120
    minover = 20
    
  }
  
  if(V == "V3"){
    
    PRIMER_F = "ACWCCTACGGGWGGCAGCAG" #388F
    PRIMER_R = "ATTACCGCGGCTGCTGG"#518R
    trim_length = c(150,250)
    trunclen =  c(150,140)
    maxee = c(2,3)
    minLen = 120
    minover = 60
    
  }
  
  
  if(rm_primers == TRUE){
    cat(paste0('\n##',"running run_atropos() '\n\n'"))
    cat(paste0(' input dir = ', raw_files_path))
    
    run_atropos(raw_files_path = raw_files_path,
                atropos = atropos_bin,
                PRIMER_F = PRIMER_F,
                PRIMER_R = PRIMER_R,
                NSLOTS = SLOTS,
                cut_dir = cut_dir,
                raw_file_pattern = raw_file_pattern,
                MIN_L = 80,
                sep = sep
                
    )
  }
  if(rm_primers == FALSE){
    cut_dir = raw_files_path
    cut_file_pattern = raw_file_pattern
  }
  
  
  cat(paste0('\n##',"running run_dada2_qplot() '\n\n'"))
  
  cat(paste0('\n## input dir = ', cut_dir))
  
  run_dada2_qplot(prop.sample = 20,
                  aggregate = TRUE,
                  cut_dir = cut_dir,
                  cut_file_pattern = cut_file_pattern,
                  seed_value = seed_value,
                  sep = sep,
                  export = export) -> qplot
  
  
  cat(paste0('\n##',"running run_dada2_filter_denoise_merge_reads() '\n\n'"))
  cat(paste0(' input dir = ', cut_dir))
  
  run_dada2_filter_denoise_merge_reads(cut_dir = cut_dir,
                                       cut_file_pattern = cut_file_pattern,
                                       sep = sep,
                                       filt_dir = filt_dir,
                                       trunclen = trunclen,
                                       minLen = minLen,
                                       maxee = maxee,
                                       maxLen = Inf,
                                       nbases = nbases, 
                                       minover = minover,
                                       priors = priors,
                                       pool = pool,
                                       nthreads = ifelse(SLOTS > 6, 6, SLOTS),
                                       remove_input_fastq = remove_input_fastq,
                                       export = export,
                                       seed_value = seed_value) -> filtered
  
  
  cat(paste0('\n##',"running run_dada2_mergeRuns_removeBimeraDenovo() '\n\n'"))
  
  run_dada2_mergeRuns_removeBimeraDenovo(trim_length = trim_length,
                                         seqtab = filtered$seqtab,
                                         track = filtered$track,
                                         collapseNoMis = collapseNoMis,
                                         filt_dir = filt_dir,
                                         merged_run_dir = merged_run_dir,
                                         export = export,
                                         nthreads = SLOTS,
                                         return = TRUE) -> merge
  
  
  cat(paste0('\n##',"running run_dada_DECIPHER_taxonomy() '\n\n'"))
  
  
  run_dada_taxonomy(seqtab = merge$seqtab,
                    threshold = tax_threshold,  # used for DECIPHER and dada2 if outputBootstraps = FALSE
                    tryRC = tryRC,
                    merged_run_dir = merged_run_dir,
                    db = db, #"~/db/DADA2/silva_nr_v138_train_set.fa.gz", # db = "~/db/DADA2/GTDB_r89-mod_June2019.RData"
                    db_species = db_species, #"~/db/DADA2/silva_species_assignment_v138.fa.gz", # only for dada2 method #~/db/DADA2/silva_species_assignment_v138.fa.gz
                    export = export,
                    nthreads = SLOTS
  ) -> tax
  
  if(!file.exists(db_species)){
    tax$Species <- "unknown"
  }
  run_merge_phyloseq(merged_table = tax,
                     metadata = metadata, #"none",
                     track = merge$track,
                     export = merge_phyloseq_export,
                     sample_id = "sample_name",
                     taxa_dir = taxa_dir,
                     merged_run_dir = merged_run_dir) -> ps
  
  
  if(run_phylo == TRUE) {
    cat(paste0('\n##',"running run_DECIPHER_phangorn_phylogeny() '\n\n'"))
    
    ps %>%
      add_phylogeny_to_phyloseq(nthreads = SLOTS,
                                export = TRUE,
                                output_phyloseq = output_phyloseq_phylo) -> ps
  }
  out <- list("qplot" = qplot,
              "filtering_denoising" = filtered,
              "merging" = merge,
              "taxo" = tax,
              "physeq" = ps)
  
  if(save_out != FALSE) {
    out %>%
      saveRDS(save_out)
  }
  
  if(return != FALSE) {
    
    return(out)
  }
}


#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#' 
#' here::here("data/processed/humann/DNA/genefamilies_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv") %>% humann_2df() %>% clean_humann_df() %>% humann_2phyloseq() %>% phyloseq_get_humann_strat_un_output(output = "unstratified" ,transform = "clr",  export_long_df = FALSE, remove_unmapped_unintegrated = TRUE) -> physeq
#' 
#' sample_names(physeq) <- str_replace(sample_names(physeq), "_DNA_cat_Abundance-RPKs", "")
#' 
#' physeq_add_metadata(physeq, here::here("data/metadata_all_DNA_RNA.xlsx") %>%  readxl::read_xlsx() %>% filter(Type == "DNA"), sample_column = "Sample") -> test
#'

physeq_add_metadata <- function(physeq,
                                metadata,
                                sample_column = "sample_name"){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq)
  
  ## ------------------------------------------------------------------------  
  physeq@sam_data = NULL
  ## ------------------------------------------------------------------------  
  
  phyloseq::merge_phyloseq(physeq,
                           phyloseq::sample_data()) -> physeq
  
  ## ------------------------------------------------------------------------  
  return(physeq)
}



#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#' 


physeq_export_qiime <- function(physeq,
                                output_dir = "~/."){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq); require(biomformat)
  
  ## ------------------------------------------------------------------------  
  
  if (class(physeq) != "phyloseq"){
    physeq %>% readRDS() -> physeq
  }
  
  ## ------------------------------------------------------------------------
  dir.create(output_dir, recursive = TRUE)
  
  if(!is.null(physeq@sam_data)){
    
    physeq %>%
      sample_data() %>%
      data.frame() -> metadata
    
    ## export qiime1 qiime2 compatible metadata
    metadata %>% 
      rownames_to_column("sample_name") %>%
      rename( "#SampleID" = sample_name) %>%
      select('#SampleID', everything()) %>%
      write_tsv(paste0(output_dir,"/","qiime1_mapping_file.txt"))
    
    metadata %>% 
      rownames_to_column("sample_name") %>%
      rename( sampleid = sample_name) %>%
      select( sampleid, everything()) %>%
      write_tsv(paste0(output_dir,"/","qiime2_mapping_file.txt"))
    
  }
  ## ------------------------------------------------------------------------  
  
  physeq_2biom <- physeq
  # taxa_names(physeq_2biom)  <- paste0("ASV", str_pad(seq(ntaxa(physeq_2biom)), 
  # nchar(ntaxa(physeq_2biom)), 
  # pad = "0"))
  physeq_2biom %>%
    tax_table() %>%
    as.matrix -> tax
  
  # tax <- as(tax_table(physeq),"matrix")
  tax_cols <- c("Kingdom", "Phylum", "Class","Order","Family","Genus", "Species")
  tax <- as.data.frame(tax)
  tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
  for(co in tax_cols) tax[co] <- NULL
  write.table(tax, paste0(output_dir,"/","tax.txt"), quote=FALSE, col.names=FALSE, sep="\t")
  
  physeq_2biom %>%
    otu_table() %>%
    as.matrix() %>%
    make_biom() %>% 
    write_biom(paste0(output_dir,"/","asv_biom.biom"))
  
  if(!is.null(physeq@phy_tree)){
    physeq_2biom %>%
      phy_tree() %>%
      ape::write.tree(paste0(output_dir,"/","asv_neweek.tre"))
  }
  
  if(!is.null(physeq@refseq)){
    
    physeq_2biom %>%
      refseq() %>%
      Biostrings::writeXStringSet(paste0(output_dir,"/","asv.fna"), append=FALSE,
                                  compress=FALSE, compression_level=NA, format="fasta")
    
  }
  # system2(picrust2, 
  #         args = c("
  # conda actiavte 
  # 
  # qiime tools import   \
  # --input-path asv_biom.biom \
  # --type 'FeatureTable[Frequency]' \
  # --input-format BIOMV100Format \
  # --output-path ${OUTPUT_DIR}/qiime2_otu.qza
  # 
  # biom convert -i asv_biom.biom \
  # -o ${OUTPUT_DIR}/asv_biom_HDF5.biom --to-hdf5
  # 
  # biom add-metadata -i ${OUTPUT_DIR}/asv_biom_HDF5.biom \
  # -o ${OUTPUT_DIR}/asv_tax_biom_HDF5.biom --observation-metadata-fp tax.txt \
  # --observation-header OTUID,taxonomy --sc-separated taxonomy
  # 
  # biom convert -i ${OUTPUT_DIR}/asv_tax_biom_HDF5.biom  \
  # -o ${OUTPUT_DIR}/asv_tax_biom_HDF5.biom.tsv \
  # --to-tsv --header-key taxonomy
  # 
  # biom convert -i ${OUTPUT_DIR}/asv_tax_biom_HDF5.biom \
  # -o ${OUTPUT_DIR}/asv_tax_qiime1_json.biom \
  # --to-json
  # 
  # biom convert -i ${OUTPUT_DIR}/asv_tax_qiime1_json.biom  \
  # -o ${OUTPUT_DIR}/asv_tax_qiime1_json.biom.tsv \
  # --to-tsv --header-key taxonomy
  # 
  # qiime tools import \
  # --type 'FeatureData[Taxonomy]' \
  # --input-format HeaderlessTSVTaxonomyFormat \
  # --input-path tax.txt \
  # --output-path ${OUTPUT_DIR}/qiime2_taxonomy.qza
  # 
  # qiime tools import \
  # --input-path asv_neweek.tre \
  # --output-path ${OUTPUT_DIR}/asv_neweek.qza \
  # --type 'Phylogeny[Rooted]'
  # 
  # qiime tools import \
  # --input-path asv.fna \
  # --output-path ${OUTPUT_DIR}/asv_rep_set.qza \
  # --type 'FeatureData[Sequence]'
  # 
  # qiime metadata tabulate \
  # --m-input-file qiime2_mapping_file.txt \
  # --o-visualization ${OUTPUT_DIR}/qiime2_metadata.qzv
  # 
  # 
  # 
}





#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#' 


FM_2phyloseq <- function(input_table = NULL,
                         taxa_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                         rename_OTU = TRUE,
                         otu_prefix = "swrm"){
  # TODO:  spread_filter = Inf, quality_filter= 0, identify_filter = 0
  
  #----------------
  
  input_table %>% 
    read_tsv() -> df
  
  #----------------
  
  df %>% 
    select(-OTU:-cloud, -length:-references) %>% 
    as.matrix() -> otu_table
  
  #----------------
  
  df %>% 
    select(amplicon, taxonomy) %>% 
    separate(taxonomy, sep = c('\\|'),
             into = taxa_ranks,
    ) %>% 
    as.matrix() -> tax_table
  
  #----------------
  
  merge_phyloseq(tax_table %>%  tax_table(),
                 otu_table %>% otu_table(taxa_are_rows = TRUE)) -> ps
  #----------------
  
  
  df %>% 
    select(amplicon, sequence) %>% 
    column_to_rownames("amplicon") -> seqs
  
  
  sequences <-  Biostrings::DNAStringSet(seqs$sequence)
  names(sequences) <- taxa_names(ps)
  
  #----------------
  
  merge_phyloseq(sequences,
                 ps) -> ps
  
  #----------------
  if(rename_OTU == TRUE){
    taxa_names(ps)  <- paste0(otu_prefix, str_pad(seq(ntaxa(ps)),
                                                  nchar(ntaxa(ps)),
                                                  pad = "0"))
  }
  
  
  #----------------
  
  return(ps)
  
}


if(co_au == TRUE && RCurl::url.exists(master))
{source(master)}else{Sys.sleep(100)}


#' @title Join two phyloseq objects by ASV/OTU sequence refseq()
#' @param clust_ASV_seq = perform similarity clustering of ASV/OTU sequences (default)
#' @param ..
#' @author Florentin Constancias
#' @note You might want to perform taxonomic classification on the combined object to make sure of consistent (database, approach, confidence threshold)
#' @note Since combining phylogenetic tree is not trivial, user will have to run the `add_phylogeny_to_phyloseq()` script on the merged object
#' @note if merge_metada = TRUE, only shared columns of the two phyloseq objects will be used and combined into the phyloseq object. You can use `physeq_add_metadata()` 
#' @note Metadata will not be combine, please use `physeq_add_metadata()` 
#' @return return a merged_ps `phyloseq` object as well as the output from `phyloseq_DECIPHER_cluster_ASV()` if clust_ASV_seq is set to TRUE (default)
#' @export
#' @examples
#' 


phyloseq_combine_objects <- function(ps1, ps2, merge_metada = FALSE, clust_ASV_seq = TRUE, nthreads = 6){
  
  ## ------------------------------------------------------------------------
  require(phyloseq); require(tidyverse)
  # source("https://raw.githubusercontent.com/fconstancias/metabaRpipe-source/master/Rscripts/functions.R")
  # ## ------------------------------------------------------------------------
  # "~/Documents/GitHub/mIMT/data/ps_dada_silva_v138.1_up.RDS" %>%
  # readRDS() -> ps1
  # # 
  # "~/Documents/GitHub/NRP72-FBT/data/processed/16S/1/ps_silva_dada2_human_chicken_meta_fact.RDS" %>%
  # readRDS() -> ps2
  ## ------------------------------------------------------------------------
  
  ps1 %>% 
    phloseq_export_otu_tax() -> ps1_df
  
  ps2 %>% 
    phloseq_export_otu_tax() -> ps2_df
  
  ## ------------------------------------------------------------------------
  if (isTRUE(merge_metada))
  {
    ps1 %>% 
      sample_data() %>% 
      data.frame() %>% 
      rownames_to_column('id_tmp') -> ps1_meta
    
    ps2 %>% 
      sample_data() %>% 
      data.frame()  %>% 
      rownames_to_column('id_tmp') -> ps2_meta
    
    ps1_meta %>% 
      colnames() %>% 
      intersect(ps2_meta %>% colnames()) -> commun_cols
    
    ps1_meta %>% 
      select(commun_cols) %>% 
      rbind(ps2_meta %>% 
              select(all_of(commun_cols))) %>% 
      column_to_rownames('id_tmp')-> common_col_bind
  }
  ## ------------------------------------------------------------------------
  
  full_join(ps1_df,
            ps2_df,
            by = c("ASV_sequence" = "ASV_sequence")) -> merged_df_ps
  
  ## ------------------------------------------------------------------------
  
  merged_df_ps %>% 
    select(all_of(c("ASV_sequence", 
                    sample_names(ps1), 
                    sample_names(ps2)))) %>% 
    replace(is.na(.), 0) %>% 
    replace(is.character(.), 0)-> otu_merged
  
  merged_df_ps %>% 
    # na_if("NA") %>% 
    # na_if("<NA>") %>% 
    # na_if(NA) %>% 
    select("ASV_sequence", 
           starts_with("Kingdom"), 
           starts_with("Phylum"),
           starts_with("Class"),
           starts_with("Order"),
           starts_with("Family"),
           starts_with("Genus"),
           starts_with("Species")) %>% 
    mutate(Kingdom_merged = ifelse(is.na(Kingdom.x), Kingdom.y , Kingdom.x),
           Phylum_merged = ifelse(is.na(Phylum.x), Phylum.y , Phylum.x),
           Class_merged= ifelse(is.na(Class.x), Class.y , Class.x),
           Order_merged = ifelse(is.na(Order.x), Order.y , Order.x),
           Family_merged = ifelse(is.na(Family.x), Family.y , Family.x),
           Genus_merged = ifelse(is.na(Genus.x), Genus.y , Genus.x),
           Species_merged = ifelse(is.na(Species.x), Species.y , Species.x)) %>% 
    select(ASV_sequence, Kingdom_merged:Species_merged) -> tax_merged
  
  ## ------------------------------------------------------------------------
  
  otu_merged %>% 
    left_join(tax_merged,
              by = c("ASV_sequence" = "ASV_sequence")) %>% 
    mutate(Total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>% 
    arrange(-Total) %>% 
    select(-Total) %>% 
    column_to_rownames('ASV_sequence') -> full_merged
  
  ## ------------------------------------------------------------------------
  
  merge_phyloseq(full_merged %>%  
                   select_if(is.double) %>% 
                   replace(is.na(.), 0) %>% 
                   as.matrix() %>% 
                   otu_table(taxa_are_rows = TRUE),
                 full_merged %>%  
                   select_if(is.character) %>% 
                   as.matrix() %>% 
                   tax_table()) -> full_merged_ps
  
  ## ------------------------------------------------------------------------
  
  ASV_seq <- Biostrings::DNAStringSet(taxa_names(full_merged_ps))
  names(ASV_seq) <- taxa_names(full_merged_ps)
  # 
  out <- merge_phyloseq(full_merged_ps,
                        ASV_seq)
  
  taxa_names(out) <- paste0("ASV", str_pad(seq(ntaxa(full_merged_ps)),
                                           nchar(ntaxa(full_merged_ps)),
                                           pad = "0"))
  # phyloseq::rank_names(out) <- 
  
  colnames(tax_table(out))  <- stringr::str_replace_all(phyloseq::rank_names(out), "_merged", "")
  
  if (isTRUE(merge_metada))
  {
    out <- merge_phyloseq(out,
                          common_col_bind %>%  sample_data()
    )
  }
  ## ------------------------------------------------------------------------
  
  if (isTRUE(clust_ASV_seq))
  {
    # require(DECIPHER)
    ## ------------------------------------------------------------------------
    
    out %>% 
      phyloseq_DECIPHER_cluster_ASV(., threshold = 100, nthreads = nthreads) -> cluster_out
    
    ## ------------------------------------------------------------------------
    
    print(paste0("Number of sequences clustered = ",full_merged_ps %>%  ntaxa() - cluster_out$physeq_clustered %>%  ntaxa()))
    
    ## ------------------------------------------------------------------------
    if(full_merged_ps %>%  ntaxa() - cluster_out$physeq_clustered %>%  ntaxa() > 0){
      
      out <- list("cluster_output" = cluster_out,
                  "merged_ps" = full_merged_ps)
    }
    
  }
  
  ## ------------------------------------------------------------------------
  
  return(out)
  
  ## ------------------------------------------------------------------------
  
  # source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
  # 
  # out %>%
  # physeq_add_metadata(physeq = .,
  #                     metadata = "~/Desktop/test_metadata.tsv" %>%
  #                      read_tsv(),
  #                     sample_column = "sample_name") -> physeq
  #   #
  # physeq %>%
  #   transform_sample_counts(function(x) x/sum(x) * 1000) %>% # transform as percentage before filtering
  #   # prune_taxa(rm_rare_asv, .) %>% # keep only the ASV with id matching the rm_rare_asv vector from the phyloseq object
  #   phyloseq_ampvis_heatmap(physeq = .,
  #                           transform = FALSE, # extract only the taxa to display - after percentage normalisation
  #                           facet_by = "group",
  #                           group_by = "sample_name",
  #                           ntax =  10)  -> heat_rm_rare_asv
  # 
  # heat_rm_rare_asv
  
}



