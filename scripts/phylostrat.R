library(dplyr)
library(readr)
library(magrittr)
library(knitr)
library(reshape2)
#library(taxizedb)
#library(phylostratr)
library(ggplot2)
#library(myTAI)
library(ggtree)

#### Run the phylostratum age analysis ####
# Create a vector of weights indicating how much to prefer a reference proteome when one is available. 
weights <- uniprot_weight_by_ref()

# Set the focal taxon to 456900, the Cyphomyrmex costatus taxon id:
focal_taxid <- '456900'

# Create a directory into which proteome sequences will be downloaded:
dir.create("./proteomeSequences2")

# Get stratified relatives of the focal taxon that are represented in UniProt,
# then select a diverse subset of 5 or fewer representatives from each stratum:
strata <- uniprot_strata(focal_taxid, 
                         from = 2) %>%
  strata_apply(f = diverse_subtree, 
               n = 5, 
               weights = weights) 

# Add in a prebuilt set of prokaryotic species, and add yeast and human:
strata <- strata %>%
  use_recommended_prokaryotes() %>%
  add_taxa(c('4932', 
             '9606'))

# Download proteomes, storing the filenames in the strata object:
strata <- strata %>%
  uniprot_fill_strata(dir = "./proteomeSequences2")

# Some proteomes produce errors, so remove those species:
warnings <- warnings() 
warnings <- names(warnings) %>%
  as.vector()
warnings <- stringr::str_extract(warnings, "\\d+")

for (i in warnings) {
  strata <- prune(strata, 
                  i, 
                  type = 'name')
}

## Replace the UniProt Cyphomyrmex costatus with my Cephalotes varians proteome:
strata@data$faa[['456900']] <- '/home/mb2337/Megan/Chapter2/CVAR/CVAR_OGS_v1.0_pep.fasta'

# Plot all selected species:
strata %>% 
  strata_convert(target = 'all', 
                 to = 'name') %>% 
  sort_strata() %>% 
  plot(cex = 0.5)

# Run orthology inference
# Use orthofinder instead of BLAST via phylostratr:
#system("./runOrthofinder/")

# BLAST against each target genome (this will take a few hours)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
                        "/programs/ncbi-blast-2.13.0/bin/", 
                        sep= .Platform$path.sep))

strataFilled <- strata_blast(strata, 
                             blast_args = list(nthreads = 64)) %>% 
  strata_besthits()

# Merge results into a single hittable
results <- merge_besthits(strataFilled)

# Add some additional data:
strataFilled <- strataFilled %>%
  add_proteome_stats()

# Plot the proteome diagnostics:
prot <- proteome_stats_table(strataFilled)
strata2 <- strata_convert(strataFilled, 
                          target = 'all', 
                          to = 'name')
g1 <- plot_proteome_stats(strata2)
g2 <- plot_proteome_lengths(strata2)
# # make interactive in browser
plotly::ggplotly(g1)
plotly::ggplotly(g2)

# Get the number of genes in each phylostratum:
ph <- stratify(results, 
               classify_by_adjusted_pvalue(0.001))
write_csv(x = ph,
          file = "./phylostratumResults.csv")
ph$locus <- sub('\\.[0-9]+',
                '',
                ph$qseqid,
                perl = TRUE)
summaryTable <- ph %>%
  dplyr::select(-qseqid) %>%
  dplyr::distinct() %>%
  dplyr::group_by(mrca_name, 
                  ps) %>%
  dplyr::summarize(n = length(ps))
summaryTable

table(ph$mrca_name)

ggplot(data = ph) +
  geom_bar(mapping = aes(x = ps))

#### Plot the tree from the analysis ####
namedTree <- strataFilled %>% 
  strata_convert(target = 'all', 
                 to = 'name') %>% 
  sort_strata() 
namedTree %>% 
  plot(cex = 0.5)

tree <- namedTree@tree
tree$tip.label <- gsub(pattern = "Cyphomyrmex costatus",
                       replacement = "Cephalotes varians",
                       tree$tip.label)

ggtree(tree) + 
  geom_tiplab(size = 1.5, 
              offset = 0.1, 
              hjust = -0.15, 
              fontface = 'italic') + 
  geom_text(aes(label = node),
            size = 1.5) +
  xlim(-0.1, 35) 




#### Combine phylostratigraphy and differential expression results ####
# Read in the phylostratigraphy results and get rid of isoforms:
phylostratigraphyResults <- read_csv("./phylostratumResults.csv") %>%
  filter(stringr::str_detect(qseqid, "-RA")) %>%
  arrange(qseqid)
phylostratigraphyResults$qseqid <- gsub(pattern = "-RA",
                                        replacement = "",
                                        x = phylostratigraphyResults$qseqid)

# Read in the differential expression results:
differentialExpressionResults <- read_csv("allDifferentialExpressionResultsAndFunctions.csv")

# Combine the two:
allResults <- right_join(differentialExpressionResults,
                         phylostratigraphyResults,
                         by = c("gene_name" = "qseqid")) %>%
  select(gene_name,
         log2FoldChange,
         padj,
         contrast,
         ps,
         mrca_name) %>%
  distinct()
allResults$mrca_name <- gsub("Cyphomyrmex costatus", 
                             "\\*Cephalotes varians\\*", 
                             allResults$mrca_name)

# Create a column classifying genes as differentially expressed or not:
allResults <- allResults %>%
  mutate(differentiallyExpressed = case_when(padj <= 0.05 ~ "yes",
                                             padj > 0.05 ~ "no"))
#### Do some plotting of expression and gene age ####
# Write a function to get proportions of genes differential expressed across strata:
getProportionDE <- function(specificContrast) {
  contrastGenes <- allResults %>%
    filter(contrast == specificContrast)
  summarizedResults <- contrastGenes %>% 
    count(ps, differentiallyExpressed) %>%
    tidyr::pivot_wider(names_from = differentiallyExpressed,
                       values_from = n)
  
  summarizedResults$ps <- as.character(summarizedResults$ps)
  
  summarizedResults <- summarizedResults %>%
    bind_rows(summarise(.,
                        across(where(is.numeric), sum),
                        across(where(is.character), ~"Total")))
  
  summarizedResults$ps <- gsub(pattern = "300",
                               replacement = "Total",
                               summarizedResults$ps)
  
  summarizedResults$proportion <- summarizedResults$yes/(summarizedResults$yes + summarizedResults$no + summarizedResults$`NA`)
  summarizedResults$contrast <- specificContrast
  return(summarizedResults)
}

# Make it safer:
possiblyGetProportionDE <- purrr::possibly(getProportionDE,
                                           otherwise = "Error")

# Run it across all contrasts:
proportions <- purrr::map(unique(allResults$contrast),
                          possiblyGetProportionDE)

proportions <- as.data.frame(do.call(rbind, proportions)) %>%
  filter(proportion != "Error")

# Do some processing for plotting:
proportions$proportion <- as.numeric(proportions$proportion)
proportions$phylostratum <- stringr::str_pad(proportions$ps, 
                                             width=2, 
                                             side = "left", 
                                             pad = "0")

# Get a table of MRCA names, instead of just strata numbers:
names <- allResults %>%
  select(ps,
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

# Plot differential expression and age:
ggplot(data = proportions) +
  geom_col(mapping = aes(x = mrca_name, 
                         y = proportion)) +
  facet_wrap(~contrast) +
  theme_bw() +
  theme(axis.text = element_text(angle = 45,
                                 hjust = 1))

#### Statistically assess the relationship ####

# A two-tailed hypergeometric test with Bonferroni correction was employed 
# to test whether gene evolution of each tissue-specific gene set differed significantly from background gene origin 
# in each of the phylostrata. (Hilgers et al. 2018)

doHypergeometricTests <- function(specificContrast) {
  pupa <- proportions %>%
    filter(contrast == specificContrast,
           ps != "Total") %>%
    select(c(ps,
             mrca_name,
             yes,
             no,
             "NA"))
  
  pupa$yes <- as.numeric(pupa$yes)
  pupa$no <- as.numeric(pupa$no)
  pupa$`NA` <- as.numeric(pupa$`NA`)
  
  pupa$total <- pupa$yes + pupa$no + pupa$`NA`
  
  pupa <- pupa %>%
    select(c(ps,
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
  pupaHyper$contrast <- specificContrast
  
  return(pupaHyper)
}

possiblyDoHypergeometricTest <- purrr::possibly(doHypergeometricTests, 
                                                otherwise = "Error")

allHypergeometricTests <- purrr::map(unique(allResults$contrast),
                                     possiblyDoHypergeometricTest)

allHypergeometricTests <- as.data.frame(do.call(rbind, 
                                                allHypergeometricTests)) %>%
  filter(ps != "Error")


plotHypergeometric <- function(specificContrast) {
  contrastData <- filter(allHypergeometricTests,
                         contrast == specificContrast)
  
  contrastData$yes <- as.numeric(contrastData$yes)
  contrastData$total <- as.numeric(contrastData$total)
  contrastData$signficant <- factor(contrastData$signficant,
                                    levels = c("MoreThanExpected", 
                                               "FewerThanExpected",
                                               "noDifference"),
                                    labels = c("More genes than expected", 
                                               "Fewer genes than expected",
                                               "No"))
  
  colors <- c("More genes than expected" = "#BAB700", 
              "Fewer genes than expected" = "#14342B",
              "No" = "#BBDFC5")
  
  ageAndSignificance <- ggplot(data = filter(contrastData,
                                             mrca_name != "Atta")) +
    geom_col(mapping = aes(x = mrca_name,
                           y = yes,
                           fill = signficant)) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(axis.text = element_text(angle = 45,
                                   hjust = 1)) + 
    ggtitle(specificContrast) + 
    labs(x = "Gene age category",
         y = "Number of differentially expressed genes",
         fill = "Significant overlap")
  
  ageAndSignificance
  return(ageAndSignificance)
}

possiblyPlotHypergeometric <- purrr::possibly(plotHypergeometric,
                                              otherwise = "Error")

allHypergeometricPlots <- purrr::map(na.omit(unique(allResults$contrast)),
                                     possiblyPlotHypergeometric)


allHypergeometricPlotsPatchwork <- patchwork::wrap_plots(allHypergeometricPlots, 
                                                         ncol = 2) + 
  patchwork::plot_annotation(title = 'Differential expression across phylostrata',
                             theme = theme(plot.title = element_text(size = 19))) + 
  patchwork::plot_layout(guides = "collect")
plot(allHypergeometricPlotsPatchwork)







# Play around with stacked bar plots 
proportionHypergeometric <- function(specificContrast) {
  contrastData <- filter(allHypergeometricTests,
                         contrast == specificContrast)
  
  contrastData$yes <- as.numeric(contrastData$yes)
  contrastData$total <- as.numeric(contrastData$total)
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
  contrastDataTest$category <- paste(contrastDataTest$signficant,
                                     " ",
                                     contrastDataTest$Type)
  contrastDataTest$category <- factor(contrastDataTest$category,
                                      levels = c("More genes than expected   total",
                                                 "No   total",
                                                 "Fewer genes than expected   total",
                                                 "More genes than expected   yes",
                                                 "No   yes",
                                                 "Fewer genes than expected   yes"),
                                      labels = c("All phylostratum genes",
                                                 "All phylostratum genes",
                                                 "All phylostratum genes",
                                                 "More genes than expected",
                                                 "No enrichment",
                                                 "Fewer genes than expected"))
  
  colors <- c("More genes than expected" = "#F98948",
              "All phylostratum genes" = "#F9EA9A",
              "No enrichment" = "#684E32",
              "Fewer genes than expected" = "#9D9B05")
  
  ageAndSignificance <- ggplot(data = contrastDataTest %>%
                                 filter(mrca_name != "Atta") %>%
                                 arrange(-Counts)) +
    geom_col(mapping = aes(x = mrca_name,
                           y = Counts,
                           fill = category),
             position = "identity") +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(axis.text = ggtext::element_markdown(angle = 45,
                                               hjust = 1)) + 
    ggtitle(specificContrast) + 
    labs(x = "Gene age category",
         y = "Number of differentially expressed genes",
         fill = "Representation of\ndifferentially expressed genes") +
    scale_y_log10()
  
  ageAndSignificance
  return(ageAndSignificance)
}

allHypergeometricProportionPlots <- purrr::map(na.omit(unique(allResults$contrast)),
                                               proportionHypergeometric)


patchwork::wrap_plots(allHypergeometricProportionPlots, 
                      ncol = 2) + 
  patchwork::plot_annotation(title = 'Differential expression across phylostrata',
                             theme = theme(plot.title = element_text(size = 19))) + 
  patchwork::plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')


#### Get general numbers for results section ####
length(unique(allResults$gene_name))

singleContrastGenes <- allResults %>%
  filter(contrast == specificContrast)

summarizedResults <- singleContrastGenes %>% 
  count(ps, differentiallyExpressed) %>%
  tidyr::pivot_wider(names_from = differentiallyExpressed,
                     values_from = n)


summarizedResults <- allResults %>%
  tidyr::pivot_wider(names_from = contrast,
                     values_from = c(ps,
                                     log2FoldChange,
                                     padj,
                                     mrca_name,
                                     differentiallyExpressed)) %>%
  count(`ps_Adult workers vs. adult soldiers`,
        `mrca_name_Adult workers vs. adult soldiers`)


# Plot all genes across phylostrata:
phylostratigraphyResultsPlotting <- phylostratigraphyResults
phylostratigraphyResultsPlotting$mrca_name <- gsub("Cyphomyrmex costatus", 
                                                   "\\*Cephalotes varians\\*", 
                                                   phylostratigraphyResultsPlotting$mrca_name)
# Do some processing for plotting:
phylostratigraphyResultsPlotting$ps <- stringr::str_pad(phylostratigraphyResultsPlotting$ps, 
                                                        width = 2, 
                                                        side = "left", 
                                                        pad = "0")

# Get a table of MRCA names, instead of just strata numbers:
names <- phylostratigraphyResultsPlotting %>%
  select(ps,
         mrca_name) %>%
  distinct() %>%
  arrange(as.numeric(ps))
names$ps <- as.character(names$ps)
names$mrca_name <- factor(names$mrca_name,
                          levels = names$mrca_name)

# Add the names to the proportions:
phylostratigraphyResultsPlotting <- left_join(phylostratigraphyResultsPlotting,
                                              names, 
                                              by = c("ps" = "ps"))

phylostratigraphyResultsPlotting$mrca_name <- factor(phylostratigraphyResultsPlotting$mrca_name.y,
                                                     levels = levels(phylostratigraphyResultsPlotting$mrca_name.y))

ggplot(data = phylostratigraphyResultsPlotting) +
  geom_bar(mapping = aes(x = mrca_name),
           fill = "#BBDFC5") +
  theme_bw() +
  theme(axis.text = ggtext::element_markdown(angle = 45,
                                             hjust = 1),
        title = ggtext::element_markdown()) + 
  labs(title = "Distribution of gene ages in the *Cephalotes varians* genome",
       x = "Gene age category",
       y = "Number of genes") 


