Off-target screening — artifact filtering note
================================================
Automated screening flagged low-support candidate loci in most samples. However,
the same loci recur across UNRELATED sample lines, which is the signature of
alignment artifacts (mapping to repetitive / low-complexity regions), not true
random integration — a genuine off-target event is private to a single sample.

Recurrent loci removed as artifacts (appearing in >=2 unrelated samples), e.g.:
  chr7:90091184 (5 samples), chr6:128806035, chr14:46050272, chr11:22923177,
  chr10:121599139 (4 samples each), and additional 2-sample loci.

After removing recurrent-across-sample loci, NO credible off-target integration
remains in any of the six samples.

Per-sample raw candidate files are provided for transparency:
  <sample>.offtarget_screen.tsv, <sample>.candidate_integration_sites.tsv
