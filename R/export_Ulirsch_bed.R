library(dplyr)

grs <- readRDS("./data/DHSmerged_ROADMAPK27ac_CAK27ac.rds")
chromosomes <- c(as.character(seq(22)), "X")

purrr::map2(grs, list("DHSmerged", "Roadmap_H3K27ac", "CA_H3K27ac"), function(gr, name) {
  df <- GenomicRanges::as.data.frame(gr) %>%
    dplyr::select(seqnames, start, end) %>%
    dplyr::mutate(seqnames = stringr::str_remove(seqnames, "^chr")) %>%
    dplyr::filter(seqnames %in% chromosomes) %>%
    dplyr::mutate(seqnames = factor(seqnames, levels = chromosomes)) %>%
    dplyr::arrange(seqnames, start)
  write.table(df, sprintf("./data/%s_Ulirsch.bed", name), row.names = F, col.names = F, quote = F, sep = "\t")
})