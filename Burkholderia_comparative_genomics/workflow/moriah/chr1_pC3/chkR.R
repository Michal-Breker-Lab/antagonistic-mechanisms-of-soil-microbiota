cat(R.version.string, "\n")
for (p in c("ape","phytools","caper","geiger","phylolm","nlme")) {
  cat(sprintf("%-10s %s\n", p, ifelse(requireNamespace(p, quietly=TRUE), "present", "MISSING")))
}
