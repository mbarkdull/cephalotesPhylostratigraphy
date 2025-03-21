library(dplyr)
library(readr)
library(magrittr)
library(knitr)
library(reshape2)
library(phylostratr)
library(ggplot2)
library(ggtree)
library(ggstatsplot)

#### Read in all the results and get a list of contrasts ####
allResults <- read_csv("allPhylostratigraphyAndDE.csv")
listOfContrasts <- unique(allResults$contrast) %>%
  na.omit() %>%
  as.character()

#### Look at all differentially expressed genes ####
hyperGeometricAllDE <- function(specificContrast) {
  # Filter to one contrast 
  contrastGenes <- allResults %>%
    filter(contrast == specificContrast)
  
  # Summarize the number of differentially expressed genes at each stratum
  summarizedResults <- contrastGenes %>% 
    count(ps, 
          differentiallyExpressed) %>%
    tidyr::pivot_wider(names_from = differentiallyExpressed,
                       values_from = n)
  
  summarizedResults <- summarizedResults %>% 
    mutate(across(where(is.numeric), 
                  tidyr::replace_na, 0))
  
  summarizedResults$ps <- as.character(summarizedResults$ps)
  
  summarizedResults <- summarizedResults %>%
    bind_rows(summarise(.,
                        across(where(is.numeric), sum),
                        across(where(is.character), ~"Total")))
  
  summarizedResults$ps <- gsub(pattern = "300",
                               replacement = "Total",
                               summarizedResults$ps)
  #summarizedResults$yes <- summarizedResults$yesShared + summarizedResults$yesUnique
  #summarizedResults$proportion <- (summarizedResults$yesUnique + summarizedResults$yesShared)/(summarizedResults$yesUnique + summarizedResults$yesShared + summarizedResults$no + summarizedResults$`NA`)
  summarizedResults$proportion <- (summarizedResults$yes)/(summarizedResults$yes + summarizedResults$no + summarizedResults$`NA`)
  proportions <- summarizedResults
  
  # Make the proportion numeric:
  proportions$proportion <- as.numeric(proportions$proportion)
  
  # Pad the phylostratum so it will sort correctly
  proportions$phylostratum <- stringr::str_pad(proportions$ps, 
                                               width=2, 
                                               side = "left", 
                                               pad = "0")
  
  # Get a table of MRCA names, instead of just strata numbers:
  names <- allResults %>%
    dplyr::select(ps,
                  mrca_name) %>%
    distinct() %>%
    arrange(as.numeric(ps))
  names$ps <- as.character(names$ps)
  names$mrca_name <- factor(names$mrca_name,
                            levels = names$mrca_name)
  
  # Add the names to the proportions:
  proportions <- left_join(proportions,
                           names, 
                           by = c("ps" = "ps"))
  proportions$mrca_name <- as.character(proportions$mrca_name) %>%
    tidyr::replace_na("Total")
  proportions$mrca_name <- factor(proportions$mrca_name,
                                  levels = unique(proportions$mrca_name))
  
  #### Do hypergeometric tests: ####
  pupa <- proportions %>%
    filter(ps != "Total") %>%
    dplyr::select(c(ps,
                    mrca_name,
                    yes,
                    no,
                    "NA"))
  
  pupa$yes <- as.numeric(pupa$yes)
  pupa$no <- as.numeric(pupa$no)
  pupa$`NA` <- as.numeric(pupa$`NA`)
  
  pupa$total <- pupa$yes + pupa$no + pupa$`NA`
  
  pupa <- pupa %>%
    dplyr::select(c(ps,
                    mrca_name,
                    yes,
                    total))
  
  pupa$totalDE <- sum(pupa$yes)
  pupa$totalOverall <- sum(pupa$total)
  pupa
  
  # hypergeometric calculation
  # q = (overlap between DE and stratum-1)
  # m = number of genes in experiment 1, differential expression
  # n = (total number of genes - m)
  # k = number of genes in experiment 2, that stratum. 
  # dhyper: Calculates the probability of obtaining exactly x successes in a sample of size k 
  # drawn from a population with m successes and n failures.
  # Probability of getting at least as many genes as I see, but no more:
  CDFHyper <- phyper(q = pupa$yes - 1, 
                     m = pupa$totalDE, 
                     n = pupa$totalOverall - pupa$totalDE, 
                     k = pupa$total)
  
  # Probability of getting exactly as many genes as I see:
  PDFHyper <- dhyper(x = pupa$yes - 1, 
                     m = pupa$totalDE, 
                     n = pupa$totalOverall - pupa$totalDE, 
                     k = pupa$total)
  
  # Probability of getting as many or more genes as I see:
  CDFHyperOver <- (1 - CDFHyper) + PDFHyper
  
  # Get the lowest of the two tails:
  raw_p_value <- pmin(CDFHyper, CDFHyperOver)*2
  
  # Get all of the data:
  pupaHyper <- cbind.data.frame(pupa, CDFHyper, CDFHyperOver, raw_p_value)
  
  pupaHyper <- pupaHyper %>%
    mutate(signficant = case_when(raw_p_value <= 0.05 & CDFHyperOver < CDFHyper ~ "MoreThanExpected",
                                  raw_p_value <= 0.05 & CDFHyperOver > CDFHyper ~ "FewerThanExpected",
                                  raw_p_value > 0.05 ~ "noDifference"))
  pupaHyper$category <- "all"
  pupaHyper$contrast <- specificContrast
  
  #### Calculate log-odds at each phylostratum: ####
  # Get just the key columns from the table
  phylostrata <- as.numeric(pupa$ps)
  quant <- pupa$yes
  sample <- pupa$total
  hit <- pupa$totalDE
  total <- pupa$totalOverall
  
  # And combine columns from table into dataset
  dataset <- cbind.data.frame(phylostrata, quant, sample, hit, total)
  
  odds_sample <- quant / (sample - quant)
  odds_rest <- (hit - quant) / (total - hit - sample + quant)
  real_log_odds <- log(odds_sample/odds_rest)
  dataset <- cbind.data.frame(dataset, odds_sample, odds_rest, real_log_odds)
  logOdds <- dataset
  rm(dataset)
  logOdds$phylostrata <- as.character(logOdds$phylostrata)
  
  #### Combine hypergeomtric results with log-odds: ####
  allStats <- full_join(pupaHyper,
                        logOdds, 
                        by = c("ps" = "phylostrata")) %>%
    dplyr::select(-c("quant",
                     "sample",
                     "hit",
                     "total.y"))
  allStats$contrast <- specificContrast
  
  return(allStats)
  
}

allHyper <- purrr::map(listOfContrasts,
                       hyperGeometricAllDE)
allHyper <- as.data.frame(do.call(rbind, allHyper)) 

# Make a nice table for export, following column names from https://academic.oup.com/mbe/article/32/2/299/1058654#supplementary-data:
library(gt)

allHyperTable <- allHyper %>%
  dplyr::select(-c("totalDE",
                   "totalOverall",
                   "category",
                   "odds_sample",
                   "odds_rest"))
colnames(allHyperTable) <- c("Phylostratum",
                             "MRCA name",
                             "Differentially expresed genes",
                             "Total genes",
                             "CDFHyper",
                             "CDFHyperOver",
                             "p-value",
                             "Direction of significance",
                             "contrast",
                             "log-odds")

hypergeometricTests <- gt(allHyperTable,
                          groupname_col = "contrast")

hypergeometricTests

gtsave(hypergeometricTests, 
       filename = "./hypergeometricTest.docx")

#### Look at uniquely differentially expressed genes ####
hyperGeometricUniqueDE <- function(specificContrast) {
  # Filter to one contrast 
  contrastGenes <- allResults %>%
    filter(contrast == specificContrast)
  
  # Summarize the number of differentially expressed genes at each stratum
  summarizedResults <- contrastGenes %>% 
    count(ps, 
          uniquelyDifferentiallyExpressed) %>%
    tidyr::pivot_wider(names_from = uniquelyDifferentiallyExpressed,
                       values_from = n)
  
  summarizedResults <- summarizedResults %>% 
    mutate(across(where(is.numeric), tidyr::replace_na, 0))
  
  summarizedResults$ps <- as.character(summarizedResults$ps)
  
  summarizedResults <- summarizedResults %>%
    bind_rows(summarise(.,
                        across(where(is.numeric), sum),
                        across(where(is.character), ~"Total")))
  
  summarizedResults$ps <- gsub(pattern = "300",
                               replacement = "Total",
                               summarizedResults$ps)
  summarizedResults$yes <- summarizedResults$yesShared + summarizedResults$yesUnique
  summarizedResults$proportion <- (summarizedResults$yesUnique + summarizedResults$yesShared)/(summarizedResults$yesUnique + summarizedResults$yesShared + summarizedResults$no + summarizedResults$`NA`)
  proportions <- summarizedResults
  
  # Make the proportion numeric:
  proportions$proportion <- as.numeric(proportions$proportion)
  
  # Pad the phylostratum so it will sort correctly
  proportions$phylostratum <- stringr::str_pad(proportions$ps, 
                                               width=2, 
                                               side = "left", 
                                               pad = "0")
  
  # Get a table of MRCA names, instead of just strata numbers:
  names <- allResults %>%
    dplyr::select(ps,
                  mrca_name) %>%
    distinct() %>%
    arrange(as.numeric(ps))
  names$ps <- as.character(names$ps)
  names$mrca_name <- factor(names$mrca_name,
                            levels = names$mrca_name)
  
  # Add the names to the proportions:
  proportions <- left_join(proportions,
                           names, 
                           by = c("ps" = "ps"))
  proportions$mrca_name <- as.character(proportions$mrca_name) %>%
    tidyr::replace_na("Total")
  proportions$mrca_name <- factor(proportions$mrca_name,
                                  levels = unique(proportions$mrca_name))
  
  # Do hypergeometric tests:
  pupa <- proportions %>%
    filter(ps != "Total") %>%
    dplyr::select(c(ps,
                    mrca_name,
                    yes,
                    yesUnique,
                    no,
                    "NA"))
  
  pupa$yesUnique <- as.numeric(pupa$yesUnique)
  pupa$yes <- as.numeric(pupa$yes)
  pupa$no <- as.numeric(pupa$no)
  pupa$`NA` <- as.numeric(pupa$`NA`)
  
  pupa$total <- pupa$yes + pupa$no + pupa$`NA`
  
  pupa <- pupa %>%
    dplyr::select(c(ps,
                    mrca_name,
                    yesUnique,
                    total))
  
  pupa$totalDE <- sum(pupa$yesUnique)
  pupa$totalOverall <- sum(pupa$total)
  pupa
  
  # hypergeometric calculation
  # q = (overlap between DE and stratum-1)
  # m = number of genes in experiment 1, differential expression
  # n = (total number of genes - m)
  # k = number of genes in experiment 2, that stratum. 
  # dhyper: Calculates the probability of obtaining exactly x successes in a sample of size k 
  # drawn from a population with m successes and n failures.
  # Probability of getting at least as many genes as I see, but no more:
  CDFHyper <- phyper(q = pupa$yesUnique - 1, 
                     m = pupa$totalDE, 
                     n = pupa$totalOverall - pupa$totalDE, 
                     k = pupa$total)
  
  # Probability of getting exactly as many genes as I see:
  PDFHyper <- dhyper(x = pupa$yesUnique - 1, 
                     m = pupa$totalDE, 
                     n = pupa$totalOverall - pupa$totalDE, 
                     k = pupa$total)
  
  # Probability of getting as many or more genes as I see:
  CDFHyperOver <- (1 - CDFHyper) + PDFHyper
  
  # Get the lowest of the two tails:
  raw_p_value <- pmin(CDFHyper, CDFHyperOver)*2
  
  # Get all of the data:
  pupaHyper <- cbind.data.frame(pupa, CDFHyper, CDFHyperOver, raw_p_value)
  
  pupaHyper <- pupaHyper %>%
    mutate(signficant = case_when(raw_p_value <= 0.05 & CDFHyperOver < CDFHyper ~ "MoreThanExpected",
                                  raw_p_value <= 0.05 & CDFHyperOver > CDFHyper ~ "FewerThanExpected",
                                  raw_p_value > 0.05 ~ "noDifference"))
  pupaHyper$category <- "unique"
  pupaHyper <- rename(pupaHyper, yes = yesUnique)
  pupaHyper$contrast <- specificContrast
  
  return(pupaHyper)
}

uniqueHyper <- purrr::map(listOfContrasts,
                          hyperGeometricUniqueDE)
uniqueHyper <- as.data.frame(do.call(rbind, uniqueHyper)) 
uniqueHyper$odds_sample <- ""
uniqueHyper$odds_rest <- ""
uniqueHyper$real_log_odds <- ""
colnames(uniqueHyper) <- c("ps",
                           "mrca_name",
                           "yes",
                           "total.x",
                           "totalDE",
                           "totalOverall",
                           "CDFHyper",
                           "CDFHyperOver",
                           "raw_p_value",
                           "signficant",
                           "category",
                           "contrast",
                           "odds_sample",
                           "odds_rest",
                           "real_log_odds")

#### Combine all statistical test results ####
hyper <- rbind(allHyper,
               uniqueHyper)

#### Plot raw gene counts for all differentially expressed genes only: ####
plotProportions <- function(specificContrast) {
  contrastData <- filter(hyper,
                         contrast == specificContrast &
                           category == "all")
  
  contrastData$yes <- as.numeric(contrastData$yes)
  contrastData$total <- as.numeric(contrastData$total.x)
  
  # Convert the signifcance info to a factor:
  contrastData$signficant <- factor(contrastData$signficant,
                                    levels = c("MoreThanExpected", 
                                               "FewerThanExpected",
                                               "noDifference"),
                                    labels = c("More genes than expected", 
                                               "Fewer genes than expected",
                                               "No"))
  contrastDataTest <- contrastData %>%
    tidyr::pivot_longer(cols = c(yes, total.x),
                        names_to = "Type",
                        values_to = "Counts")
  
  contrastDataTest$resultCategory <- paste(contrastDataTest$signficant,
                                           " ",
                                           contrastDataTest$Type,
                                           " ",
                                           contrastDataTest$category)
  
  contrastDataTest$resultCategory <- factor(contrastDataTest$resultCategory,
                                            levels = c("More genes than expected   total   all",
                                                       "No   total   all",
                                                       "Fewer genes than expected   total   all",
                                                       "No   total   unique",
                                                       "More genes than expected   total   unique",
                                                       "Fewer genes than expected   total   unique",
                                                       
                                                       "More genes than expected   yes   all",
                                                       "No   yes   all",
                                                       "Fewer genes than expected   yes   all",
                                                       
                                                       "No   yes   unique",
                                                       "More genes than expected   yes   unique",
                                                       "Fewer genes than expected   yes   unique"),
                                            
                                            labels = c("All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       
                                                       "Differentially expressed genes are over-enriched",
                                                       "Differentially expressed genes are not enriched",
                                                       "Differentially expressed genes are under-enriched",
                                                       
                                                       "Uniquely differentially expressed genes are not enriched",
                                                       "Uniquely differentially expressed genes are over-enriched",
                                                       "Uniquely differentially expressed genes are under-enriched"))
  
  colors <- c("All genes at this phylostratum" = "#e6ded3",
              "Differentially expressed genes are over-enriched" = "#783064",
              "Differentially expressed genes are not enriched" = "#c2b19b",
              "Differentially expressed genes are under-enriched" = "#cfaec6")
  
  ageAndSignificance <- ggplot(data = contrastDataTest %>%
                                 filter(mrca_name != "Atta") %>%
                                 arrange(-Counts)) +
    geom_col(mapping = aes(x = mrca_name,
                           y = Counts,
                           fill = resultCategory),
             position = "identity") +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(axis.text = ggtext::element_markdown(angle = 45,
                                               hjust = 1),
          legend.title = element_blank()) + 
    ggtitle(paste("Differential expression in",
                  tolower(specificContrast))) + 
    labs(x = "Gene age category",
         y = "Number of differentially expressed genes",
         fill = "Representation of\ndifferentially expressed genes") +
    scale_y_log10(expand = c(0, 0))
  
  ageAndSignificance
}

plots <- purrr::map(listOfContrasts,
                    plotProportions)

allPlotsPatchwork <- patchwork::wrap_plots(plots, 
                                           ncol = 2) + 
  patchwork::plot_annotation(title = 'Differential expression across phylostrata',
                             theme = theme(plot.title = element_text(size = 19))) +
  patchwork::plot_annotation(tag_levels = 'A',
                             tag_suffix = '.') + 
  patchwork::plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        plot.tag.position = c(0.05, 0.985)) 

plot(allPlotsPatchwork)

barPlots <- (allPlotsPatchwork[[1]] + allPlotsPatchwork[[2]]) + 
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

barPlots

ggsave(filename = "enrichmentOfPhylostrata.png",
       width = 12,
       height = 6,
       units = "in",
       dpi = 900)

#### Plot log-odds for all differentially expressed genes ####
plotLogOdds <- function(specificContrast) {
  contrastData <- filter(allHyper,
                         contrast == specificContrast &
                           category == "all")
  
  contrastData$yes <- as.numeric(contrastData$yes)
  contrastData$total.x <- as.numeric(contrastData$total.x)
  
  # Convert the signifcance info to a factor:
  contrastData$resultCategory <- factor(contrastData$signficant,
                                        levels = c("MoreThanExpected", 
                                                   "FewerThanExpected",
                                                   "noDifference"),
                                        labels = c("Differentially expressed genes are over-enriched", 
                                                   "Differentially expressed genes are under-enriched",
                                                   "Differentially expressed genes are neither\nover- nor under-enriched"))
  
  colors <- c("Differentially expressed genes are over-enriched" = "#7e2f87",
              "Differentially expressed genes are neither\nover- nor under-enriched" = "#c2b19b",
              "Differentially expressed genes are under-enriched" = "#cfaec6")
  
  ageAndSignificance <- ggplot(data = contrastData %>%
                                 filter(mrca_name != "Atta")) +
    geom_col(mapping = aes(x = mrca_name,
                           y = real_log_odds,
                           fill = resultCategory),
             position = "identity") +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(axis.text = ggtext::element_markdown(angle = 45,
                                               hjust = 1),
          legend.title = element_blank()) + 
    ggtitle(paste("Ages of genes differentially expressed\nin",
                  tolower(specificContrast))) + 
    labs(x = "Gene age category",
         y = "Log-odds of differential expression",
         fill = "Representation of\ndifferentially expressed genes") 
  
  ageAndSignificance
}

plots <- purrr::map(listOfContrasts,
                    plotLogOdds)

allPlotsPatchwork <- patchwork::wrap_plots(plots, 
                                           ncol = 2) + 
  patchwork::plot_annotation(title = 'Differential expression across phylostrata',
                             theme = theme(plot.title = element_text(size = 19))) +
  patchwork::plot_annotation(tag_levels = 'A',
                             tag_suffix = '.') + 
  patchwork::plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        plot.tag.position = c(0.05, 0.985)) 

plot(allPlotsPatchwork)



barPlots <- (allPlotsPatchwork[[1]] / allPlotsPatchwork[[2]]) + 
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Write a function to make scales and axes consistent across plots:
makePlotsConsistent <- function(plot){
  # Get the total number of plots that are combined together:
  num_plots <- length(plot)
  
  # Fix x limits and bubble sizes:
  # Get the minimum and maximum values of geneRatioDecimal for each plot:
  yLimits <- lapply(1:num_plots, function(x) ggplot_build(plot[[x]])$layout$panel_scales_y[[1]]$range$range)
  # Get the minimum and maximum x values for each plot:
  minY <- min(unlist(yLimits))
  maxY <- max(unlist(yLimits))
  
  
  plot & 
    ylim(floor(minY), 
         ceiling(maxY)) 
}

# Plot all results with consistent scales and a single legend:
logOddsPlotsList <- makePlotsConsistent(barPlots)

logOddsPlotsList

ggsave(filename = "enrichmentOfPhylostrataLogOdds.png",
       width = 6,
       height = 12,
       units = "in",
       dpi = 900)

ggsave(filename = "enrichmentOfPhylostrataLogOdds.pdf",
       width = 6,
       height = 12,
       units = "in")

#### Plot for all and uniquely DE genes: ####
plotStackedProportions <- function(specificContrast) {
  contrastData <- filter(hyper,
                         contrast == specificContrast)
  
  contrastData$yes <- as.numeric(contrastData$yes)
  contrastData$total <- as.numeric(contrastData$total)
  
  # Convert the signifcance info to a factor:
  contrastData$signficant <- factor(contrastData$signficant,
                                    levels = c("MoreThanExpected", 
                                               "FewerThanExpected",
                                               "noDifference"),
                                    labels = c("More genes than expected", 
                                               "Fewer genes than expected",
                                               "No"))
  contrastDataTest <- contrastData %>%
    tidyr::pivot_longer(cols = c(yes, total),
                        names_to = "Type",
                        values_to = "Counts")
  
  contrastDataTest$resultCategory <- paste(contrastDataTest$signficant,
                                           " ",
                                           contrastDataTest$Type,
                                           " ",
                                           contrastDataTest$category)
  
  contrastDataTest$resultCategory <- factor(contrastDataTest$resultCategory,
                                            levels = c("More genes than expected   total   all",
                                                       "No   total   all",
                                                       "Fewer genes than expected   total   all",
                                                       "No   total   unique",
                                                       "More genes than expected   total   unique",
                                                       "Fewer genes than expected   total   unique",
                                                       
                                                       "More genes than expected   yes   all",
                                                       "No   yes   all",
                                                       "Fewer genes than expected   yes   all",
                                                       
                                                       "No   yes   unique",
                                                       "More genes than expected   yes   unique",
                                                       "Fewer genes than expected   yes   unique"),
                                            
                                            labels = c("All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       "All genes at this phylostratum",
                                                       
                                                       "Differentially expressed genes are over-enriched",
                                                       "Differentially expressed genes are not enriched",
                                                       "Differentially expressed genes are under-enriched",
                                                       
                                                       "Uniquely differentially expressed genes are not enriched",
                                                       "Uniquely differentially expressed genes are over-enriched",
                                                       "Uniquely differentially expressed genes are under-enriched"))
  
  colors <- c("All genes at this phylostratum" = "#e8e6cc",
              "Differentially expressed genes are over-enriched" = "#f5b83d",
              "Differentially expressed genes are not enriched" = "#c7c181",
              "Differentially expressed genes are under-enriched" = "#74cc9d",
              "Uniquely differentially expressed genes are not enriched" = "#968e36",
              "Uniquely differentially expressed genes are over-enriched" = "#d18f08",
              "Uniquely differentially expressed genes are under-enriched" = "#49a876")
  
  ageAndSignificance <- ggplot(data = contrastDataTest %>%
                                 filter(mrca_name != "Atta") %>%
                                 arrange(-Counts)) +
    geom_col(mapping = aes(x = mrca_name,
                           y = Counts,
                           fill = resultCategory),
             position = "identity") +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(axis.text = ggtext::element_markdown(angle = 45,
                                               hjust = 1),
          legend.title = element_blank()) + 
    ggtitle(paste("Differential expression in",
                  tolower(specificContrast))) + 
    labs(x = "Gene age category",
         y = "Number of differentially expressed genes",
         fill = "Representation of\ndifferentially expressed genes") +
    scale_y_log10(expand = c(0, 0))
  
  ageAndSignificance
}

plots <- purrr::map(listOfContrasts,
                    plotStackedProportions)

allPlotsPatchwork <- patchwork::wrap_plots(plots, 
                                           ncol = 2) + 
  patchwork::plot_annotation(title = 'Differential expression across phylostrata',
                             theme = theme(plot.title = element_text(size = 19))) +
  patchwork::plot_annotation(tag_levels = 'A',
                             tag_suffix = '.') + 
  patchwork::plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom',
        plot.tag.position = c(0.05, 0.985)) 

plot(allPlotsPatchwork)

#### Do phylostrata differ in magnitude of differential expression or in their overall expression level? ####
# Read in the phylostratigraphy results and get rid of isoforms:
phylostratigraphyResults <- read_csv("allPhylostratigraphyAndDE.csv") %>%
  dplyr::select(gene_name,
                ps,
                mrca_name)

# Read in the differential expression results:
differentialExpressionResults <- read_csv("allDifferentialExpressionResultsAndFunctions.csv")

# Combine the two:
allResults <- right_join(differentialExpressionResults,
                         phylostratigraphyResults,
                         by = c("gene_name" = "gene_name")) %>%
  dplyr::select(gene_name,
                baseMean,
                log2FoldChange,
                padj,
                ps,
                mrca_name, 
                contrast) %>%
  distinct()
allResults$mrca_name <- gsub("Cyphomyrmex costatus", 
                             "\\*Cephalotes varians\\*", 
                             allResults$mrca_name)

# Quickly generate volcano plots for each phylostratum:
ggplot(data = allResults) +
  geom_point(mapping = aes(x = log2FoldChange, 
                           y = padj),
             size = 0.5) +
  scale_y_continuous(trans = scales::compose_trans("log10", 
                                           "reverse"),
                     labels = scales::label_log()) +
  facet_wrap(~ mrca_name)

# Do a partial Spearman correlation to control for expression levels (in the function):
library(ppcor)

# Write a function to get correlation coefficients and plot expression levels and differential expression:
expressionPlotsContrast <- function(specificContrast) {
  allResults <- filter(allResults, 
                       contrast == specificContrast) %>%
    distinct()
  
  spearmanMagnitudeDE <- cor.test(allResults$ps,
                                  abs(allResults$log2FoldChange),
                                  method = "spearman",
                                  exact = TRUE)
  spearmanMagnitudeDE
  
  spearmanTotalExpression <- cor.test(allResults$ps,
                                      allResults$baseMean,
                                      method = "spearman",
                                      exact = TRUE)
  spearmanTotalExpression
  
  test <- allResults %>% 
    na.omit()
  
  # Test the correlation between x and y, while controlling for z:
  partialSpearman <- pcor.test(x = test$ps,
                               y = abs(test$log2FoldChange),
                               z = test$baseMean,
                               method = "spearman")
  
  partialSpearman
  
  # Group by phylostratum and get the mean expression:
  phylostratumExpression <- allResults %>% 
    group_by(ps) %>% 
    filter(padj < 0.05) %>%
    summarise(meanMagnitude = mean(abs(log2FoldChange),
                                   na.rm = TRUE),
              meanExpression = mean(baseMean,
                                    na.rm = TRUE)) %>%
    as.data.frame() %>%
    left_join(allResults) %>%
    dplyr::select(-c(baseMean,
              gene_name,
              log2FoldChange,
              padj)) %>%
    distinct()
  
  phylostratumExpression$mrca_name <- factor(phylostratumExpression$mrca_name,
                                             levels = phylostratumExpression$mrca_name)
  
  # Plot:
  magnitudeDEColor <- "#8b9448"
  totalExpressionColor <- "#c97718"
  
  contrastForPlotTitle <- gsub(pattern = "vs.",
                               replacement = "and",
                               specificContrast) %>%
    tolower()
  
  
  ggplot(data = phylostratumExpression,
         mapping = aes(x = mrca_name)) +
    geom_point(mapping = aes(y = meanMagnitude),
               color = magnitudeDEColor) +
    geom_point(mapping = aes(y = meanExpression/(max(meanExpression)/max(meanMagnitude))),
               color = totalExpressionColor) +
    theme_bw() +
    annotate(geom = "label",
             x = 1, 
             y = max(phylostratumExpression$meanMagnitude) * 1.4, 
             label = paste("Correlation between phylostratum and magnitude of differential expression,\n controlling for expression level (partial Spearman's rank correlation):\n",
                           "ρ: ",
                           format(round(partialSpearman$estimate, 
                                        3), 
                                  nsmall = 2),
                           "; p-value: ",
                           formatC(partialSpearman$p.value, 
                                   format = "e", 
                                   digits = 2),
                           
                           sep = ""), 
             fill = "white",
             hjust = 0) + 
    theme(axis.text = ggtext::element_markdown(angle = 45,
                                               hjust = 1),
          axis.title.y = element_text(color = magnitudeDEColor),
          axis.title.y.right = element_text(color = totalExpressionColor)) +
    labs(title = paste("Relationship between gene age and expression patterns across\n",
                       contrastForPlotTitle),
         x = "Phylostratum") + 
    scale_y_continuous(name = "Mean magnitude of differential expression",
                       limits = c(0, max(phylostratumExpression$meanMagnitude) * 1.5),
                       sec.axis = sec_axis(~.*(max(phylostratumExpression$meanExpression)/max(phylostratumExpression$meanMagnitude)), 
                                           name = "Mean total expression level")) 
}

expressionPlotsContrast("Worker pupae vs. soldier pupae")

ggsave(filename = "./magnitudeOfDEAndTotalExpressionPupae.png",
       height = 6,
       width = 8,
       units = "in",
       dpi = 900)

spearmanPlots <- purrr::map(na.omit(unique(allResults$contrast)),
                            expressionPlotsContrast)

patchwork::wrap_plots(spearmanPlots)

ggsave(filename = "./partialSpearmanExpressionLevelsAndPhylostratum.png",
       height = 10,
       width = 14,
       units = "in",
       dpi = 900)

#### Combine plots ####
allPlotsPatchwork[[1]]
allPlotsPatchwork[[2]]
spearmanPlots[[1]]
spearmanPlots[[2]]


library(patchwork)
barPlots <- (allPlotsPatchwork[[1]] + allPlotsPatchwork[[2]]) + 
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "none")

dotPlots <- (spearmanPlots[[1]] + spearmanPlots[[2]]) + 
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "none")

barLegend <- ggpubr::get_legend(allPlotsPatchwork[[1]]) %>%
  ggpubr::as_ggplot()
dotLegend <- ggpubr::get_legend(spearmanPlots[[1]]) %>%
  ggpubr::as_ggplot()

(barPlots / barLegend / dotPlots) + 
  patchwork::plot_layout(heights = c(1, 0.35, 1)) 

ggsave("./combinedPlots.png", 
       width = 35, 
       height = 25, 
       units = "cm", 
       dpi = 1800)



# Group by phylostratum and get the mean expression:
phylostratumExpression <- allResults %>% 
  group_by(ps) %>% 
  summarise(meanExpression = mean(baseMean,
                                  na.rm = TRUE)) %>%
  as.data.frame() %>%
  left_join(allResults) %>%
  select(-c(baseMean,
            gene_name)) %>%
  distinct()

phylostratumExpression$mrca_name <- factor(phylostratumExpression$mrca_name,
                                           levels = phylostratumExpression$mrca_name)

# Check significance of the correlation between continuous + ordinal variables:
meanSpearmanCoefficient <- cor.test(phylostratumExpression$ps,
                                phylostratumExpression$meanExpression,
                                method = "spearman",
                                exact = TRUE)
meanSpearmanCoefficient

# Plot:
ggplot(data = phylostratumExpression,
       mapping = aes(x = mrca_name, 
                     y = meanExpression)) +
  geom_point() +
  theme_bw() +
  annotate(geom = "label",
           x = 16, 
           y = 3250, 
           label = paste("Spearman's rank correlation:\np-value:",
                         formatC(spearmanCoefficient$p.value, 
                                 format = "e", 
                                 digits = 2),
                         "\nρ:",
                         format(round(spearmanCoefficient$estimate, 
                                      2), 
                                nsmall = 2)), 
           fill = "white",
           hjust = 0) + 
  theme(axis.text = ggtext::element_markdown(angle = 45,
                                             hjust = 1)) +
  labs(x = "Phylostratum",
       y= "Mean expression level")

ggsave(filename = "./expressionLevelsAverage.png",
       height = 6,
       width = 8,
       units = "in",
       dpi = 900)

# Plot as a boxplot:
allResults$mrca_name <- factor(allResults$mrca_name,
                               levels = phylostratumExpression$mrca_name)

ggplot(data = allResults,
       mapping = aes(x = mrca_name, 
                     y = baseMean)) +
  geom_boxplot(size = 0.5,
               outlier.size = 0.3,
               fill = "#e8e6cc") +
 # geom_jitter(size = 0.5,
  #           alpha = 0.05) +
  theme_bw() +
  annotate(geom = "label",
           x = 3, 
           y = 900000, 
           label = paste("Spearman's rank correlation:\np-value:",
                         formatC(spearmanCoefficient$p.value, 
                                 format = "e", 
                                 digits = 2),
                         "\nρ:",
                         format(round(spearmanCoefficient$estimate, 
                                      2), 
                                nsmall = 2)), 
           fill = "white",
           hjust = 0) + 
  theme(axis.text = ggtext::element_markdown(angle = 45,
                                             hjust = 1)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits = c(0, 3015713),
                     expand = c(0, 0)) +
  labs(x = "Phylostratum",
       y= "Expression level")

ggsave(filename = "./expressionLevelsBoxplot.png",
       height = 6,
       width = 8,
       units = "in",
       dpi = 900)

# Or try plotting with ggscatterstats:
# *r* gives effect size
ggscatterstats(phylostratumExpression,
               x = ps,
               y = meanExpression,
               label.var = mrca_name,
               type = "nonparametric")



#### Look at genes that are up or downregulated ####
hyperGeometricDirectional <- function(specificContrast, direction ) {
  if (direction == "up") {
    contrastGenes <- allResults %>%
      filter(contrast == specificContrast) %>%
      mutate(differentiallyExpressed = case_when(differentiallyExpressed == "yes" &
                                                   log2FoldChange < 0 ~ "no",
                                                 differentiallyExpressed == "yes" &
                                                   log2FoldChange > 0 ~ "yes",
                                                 TRUE ~ "no"))
  } else {
    contrastGenes <- allResults %>%
      filter(contrast == specificContrast) %>%
      mutate(differentiallyExpressed = case_when(differentiallyExpressed == "yes" &
                                                   log2FoldChange < 0 ~ "yes",
                                                 differentiallyExpressed == "yes" &
                                                   log2FoldChange > 0 ~ "no",
                                                 TRUE ~ "no"))
  }
  
  
  
  
  # Summarize the number of differentially expressed genes at each stratum
  summarizedResults <- contrastGenes %>% 
    count(ps, 
          differentiallyExpressed) %>%
    tidyr::pivot_wider(names_from = differentiallyExpressed,
                       values_from = n)
  
  summarizedResults <- summarizedResults %>% 
    mutate(across(where(is.numeric), tidyr::replace_na, 0))
  
  summarizedResults$ps <- as.character(summarizedResults$ps)
  
  summarizedResults <- summarizedResults %>%
    bind_rows(summarise(.,
                        across(where(is.numeric), sum),
                        across(where(is.character), ~"Total")))
  
  summarizedResults$ps <- gsub(pattern = "300",
                               replacement = "Total",
                               summarizedResults$ps)
  
  summarizedResults$proportion <- (summarizedResults$yes)/(summarizedResults$yes + summarizedResults$no)
  proportions <- summarizedResults
  
  # Make the proportion numeric:
  proportions$proportion <- as.numeric(proportions$proportion)
  
  # Pad the phylostratum so it will sort correctly
  proportions$phylostratum <- stringr::str_pad(proportions$ps, 
                                               width=2, 
                                               side = "left", 
                                               pad = "0")
  
  # Get a table of MRCA names, instead of just strata numbers:
  names <- allResults %>%
    dplyr::select(ps,
                  mrca_name) %>%
    distinct() %>%
    arrange(as.numeric(ps))
  names$ps <- as.character(names$ps)
  names$mrca_name <- factor(names$mrca_name,
                            levels = names$mrca_name)
  
  # Add the names to the proportions:
  proportions <- left_join(proportions,
                           names, 
                           by = c("ps" = "ps"))
  proportions$mrca_name <- as.character(proportions$mrca_name) %>%
    tidyr::replace_na("Total")
  proportions$mrca_name <- factor(proportions$mrca_name,
                                  levels = unique(proportions$mrca_name))
  
  # Do hypergeometric tests:
  pupa <- proportions %>%
    filter(ps != "Total") %>%
    dplyr::select(c(ps,
                    mrca_name,
                    yes,
                    no))
  
  pupa$yes <- as.numeric(pupa$yes)
  pupa$no <- as.numeric(pupa$no)
  
  pupa$total <- pupa$yes + pupa$no
  
  pupa <- pupa %>%
    dplyr::select(c(ps,
                    mrca_name,
                    yes,
                    total))
  
  pupa$totalDE <- sum(pupa$yes)
  pupa$totalOverall <- sum(pupa$total)
  pupa
  
  # hypergeometric calculation
  # q = (overlap between DE and stratum-1)
  # m = number of genes in experiment 1, differential expression
  # n = (total number of genes - m)
  # k = number of genes in experiment 2, that stratum. 
  # dhyper: Calculates the probability of obtaining exactly x successes in a sample of size k 
  # drawn from a population with m successes and n failures.
  # Probability of getting at least as many genes as I see, but no more:
  CDFHyper <- phyper(q = pupa$yes - 1, 
                     m = pupa$totalDE, 
                     n = pupa$totalOverall - pupa$totalDE, 
                     k = pupa$total)
  
  # Probability of getting exactly as many genes as I see:
  PDFHyper <- dhyper(x = pupa$yes - 1, 
                     m = pupa$totalDE, 
                     n = pupa$totalOverall - pupa$totalDE, 
                     k = pupa$total)
  
  # Probability of getting as many or more genes as I see:
  CDFHyperOver <- (1 - CDFHyper) + PDFHyper
  
  # Get the lowest of the two tails:
  raw_p_value <- pmin(CDFHyper, CDFHyperOver)*2
  
  # Get all of the data:
  pupaHyper <- cbind.data.frame(pupa, CDFHyper, CDFHyperOver, raw_p_value)
  
  pupaHyper <- pupaHyper %>%
    mutate(signficant = case_when(raw_p_value <= 0.05 & CDFHyperOver < CDFHyper ~ "MoreThanExpected",
                                  raw_p_value <= 0.05 & CDFHyperOver > CDFHyper ~ "FewerThanExpected",
                                  raw_p_value > 0.05 ~ "noDifference"))
  pupaHyper$category <- "all"
  pupaHyper$contrast <- specificContrast
  return(pupaHyper)
}