colors <- c("All phylostratum genes" = "#d6c6b8",
            "Differentially expressed genes are over-enriched" = "#b3c480",
            "Differentially expressed genes are not enriched" = "#8f6d4f",
            "Differentially expressed genes are under-enriched" = "#4c84e6",
            "Uniquely differentially expressed genes are not enriched" = "#6e4827",
            "Uniquely differentially expressed genes are over-enriched" = "#779714",
            "Uniquely differentially expressed genes are under-enriched" = "#043487")

ageAndSignificance <- ggplot(data = contrastDataTest %>%
                               filter(mrca_name != "Atta",
                                      resultCategory == "All phylostratum genes") %>%
                               arrange(-Counts)) +
  geom_col(mapping = aes(x = mrca_name,
                         y = Counts,
                         fill = resultCategory),
           position = "identity") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.text = ggtext::element_markdown(angle = 45,
                                             hjust = 1),
        legend.title = element_blank(),
        legend.position = "none") + 
  ggtitle(paste("Differential expression in",
                tolower(specificContrast))) + 
  labs(x = "Gene age category",
       y = "Number of differentially expressed genes",
       fill = "Representation of\ndifferentially expressed genes") +
  scale_y_continuous(expand = c(0, 0))

ageAndSignificance

ggsave(filename = "nonLogGenesPresentation.png",
       plot = ageAndSignificance,
       width = 12,
       height = 8, 
       units = "in", 
       dpi = 600)
