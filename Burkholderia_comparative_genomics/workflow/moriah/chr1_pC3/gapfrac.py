import statistics as st
p="/sci/backup/ofinkel/moshea/burkholderia_c3/results/chr1_core_763.aln"
name=None; g=0; n=0; out=[]
for ln in open(p):
    if ln.startswith(">"):
        if name: out.append((g/n, name))
        name=ln[1:].strip(); g=0; n=0
    else:
        s=ln.strip(); n+=len(s); g+=s.count("-")+s.upper().count("N")
if name: out.append((g/n, name))
out.sort(reverse=True)
print("n_seq", len(out))
print("median gap%%: %.2f   >50%%: %d   >90%%: %d" % (st.median(f for f,_ in out)*100,
      sum(1 for f,_ in out if f>.5), sum(1 for f,_ in out if f>.9)))
for f,nm in out[:20]: print("%6.2f  %s" % (f*100, nm))
