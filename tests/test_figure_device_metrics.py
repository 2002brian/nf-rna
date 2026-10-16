"""Stale text metrics across graphics devices (Phase 4B): the guard that keeps PDF/TIFF layouts consistent."""

from __future__ import annotations

from test_figures import run_r


def test_every_device_the_helpers_open_is_primed_before_anything_is_drawn_or_measured(tmp_path):
    result = run_r('''
calls <- 0L; real <- nf_rna_prime_device_metrics
nf_rna_prime_device_metrics <- function(family="") { calls <<- calls + 1L; real(family) }
ctx <- nf_rna_figure_context(list(figure_style=list(preset="legacy", profile="general")))
p <- ggplot(data.frame(x=1:3, y=1:3), aes(x, y)) + geom_point()
dir.create("out")
before <- calls
nf_rna_save_figure(ctx, p, file.path("out", "f"), 60, 50, "out")
saved <- calls - before
before <- calls
nf_rna_with_probe_device(ctx, function() NULL)
cat(toJSON(list(per_save = saved, per_probe = calls - before, formats = length(ctx$records[[1]]$files)), auto_unbox = TRUE))''', tmp_path)
    assert result["formats"] == 3 and result["per_save"] == 3 and result["per_probe"] == 1  # PNG, TIFF, PDF and the probe device


def test_priming_removes_the_stale_metric_that_follows_a_device_of_another_resolution(tmp_path):
    # Without priming, the first 7 pt text height on a device opened right after another device (other resolution) that
    # ended with the same font is intermittently returned in the previous device's units.  A garbage collection or a short
    # busy period between the devices makes that frequent, so 20 rounds expose it; with priming every round must agree.
    result = run_r('''
fam <- "sans"
open <- function(d) switch(d, png=png(tempfile(fileext=".png"), width=630, height=665, res=150, type="cairo", family=fam),
  tiff=tiff(tempfile(fileext=".tiff"), width=1260, height=1330, res=300, type="cairo", compression="lzw", family=fam),
  pdf=cairo_pdf(tempfile(fileext=".pdf"), width=4.2, height=4.45, family=fam, onefile=FALSE))
m7 <- function() grid::convertHeight(grid::grobHeight(grid::textGrob("M", gp=grid::gpar(fontsize=7, fontfamily=fam))), "pt", valueOnly=TRUE)
run <- function(prime) {
  seen <- list()
  for (round in 1:20) for (d in c("pdf", "png", "tiff", "pdf")) {
    open(d); if (prime) nf_rna_prime_device_metrics(fam)
    seen[[d]] <- c(seen[[d]], round(m7(), 3))
    invisible(gc()); t0 <- proc.time()[[3]]; while (proc.time()[[3]] - t0 < 0.02) NULL
    dev.off()
  }
  vapply(seen, function(v) length(unique(v)), 1L)
}
cat(toJSON(list(primed = run(TRUE)), auto_unbox = FALSE))''', tmp_path)
    assert all(count == 1 for count in result["primed"]), result  # one value per device type in every round
