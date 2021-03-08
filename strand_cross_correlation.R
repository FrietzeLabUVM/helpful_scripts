library(seqsetvis)
library(peakrefine)
library(ggplot2)

#measure strand cross correlation for a bam file within peaks called
bam_file = "/slipstream/galaxy/uploads/working/qc_framework/output/MCF10A_H3K4ME3_R1/MCF10A_H3K4ME3_R1.bam"
np_file = "/slipstream/galaxy/uploads/working/qc_framework/output/MCF10A_H3K4ME3_R1/MCF10A_H3K4ME3_R1_peaks.narrowPeak"

qgr = seqsetvis::easyLoad_narrowPeak(bam_file)[[1]]
options(mc.cores = 20)
scc_res = peakrefine::calcSCCMetrics(np_file, qgr, frag_sizes = seq(50, 300, 10), n_splits = 20)

scc_res$read_length
p1 = ggplot(scc_res$average_correlation, aes(x = shift, y = correlation)) + geom_path() +
    annotate("line", x = rep(scc_res$read_length, 2), y = range(scc_res$average_correlation$correlation), linetype = 2, color = "red") +
    annotate("text", x = scc_res$read_length + 3, y = mean(scc_res$average_correlation$correlation), label = paste0("read length\n", scc_res$read_length), hjust = 0, vjust = .5, color = "red") +
    annotate("line", x = rep(scc_res$fragment_length, 2), y = range(scc_res$average_correlation$correlation), linetype = 2, color = "blue") +
    annotate("text", x = scc_res$fragment_length + 3, y = mean(scc_res$average_correlation$correlation), label = paste0("fragment length estimate\n", scc_res$fragment_length), hjust = 0, vjust = .5, color = "blue") +
    labs(title = "Strand Cross Correlation (SCC)", subtitle = paste(basename(bam_file), "@", basename(np_file)))


# plot(scc_res$read_correlation$correlation,
#      scc_res$stable_fragment_correlation$correlation, xlab = "read length correlation", ylab = "fragment length correlation") 
# rect(xleft = .9, xright = 1.03, ybottom = min(scc_res$read_correlation$correlation), ytop = max(scc_res$read_correlation$correlation), col = "#FF000055")
# title("Peak in the red zone have quite high correlation at read length and are likely artifacts")

p2 = ggplot() +
    annotate("rect", xmin = .9, xmax = 1.03, 
             ymin = min(scc_res$read_correlation$correlation), 
             ymax = max(scc_res$read_correlation$correlation), fill = "#FF000055") +
    annotate("point", x = scc_res$read_correlation$correlation, y =  scc_res$stable_fragment_correlation$correlation, size = .7, alpha = .5) +
    labs(title = "Peaks in the red zone have quite high\ncorrelation at read length and are likely artifacts",
         x = "read length SCC",
         y = "fragment length estimate SCC")
    
cowplot::plot_grid(p1, p2)
