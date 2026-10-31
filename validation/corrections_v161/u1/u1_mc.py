# U-1: size of the differences against the Monte Carlo error that fgsea reports (log2err), and
# against the seed-1 vs seed-42 variation of the reference.
import csv,glob,math
from pathlib import Path
V=Path("/home/brian/nf-rna-validation")
l2=Path(glob.glob(str(V/"v161/nfrna/D2/runs/FIX-D2/*/downstream/l2"))[0])
rd=lambda p,k:{r[k]:r for r in csv.DictReader(open(p),delimiter="\t")}
for coll,sub in (("GO_BP","BP"),("GO_MF","MF")):
    f=rd(l2/"enrichment/gsea_go/LP_vs_ML"/sub/"all_terms.tsv","ID")
    r=rd(V/"v1.6.0/reference/outputs/D2/contrasts/LP_vs_ML/gsea"/coll/"seed_1.tsv","set_id")
    r42=rd(V/"v1.6.0/reference/outputs/D2/contrasts/LP_vs_ML/gsea"/coll/"seed_42.tsv","set_id")
    worst=0;na=0;n=0;nes=[];seedratio=[];near=[]
    for k,b in r.items():
        a=f[k];p1,p0=float(a["pvalue"]),float(b["pval"])
        if abs(p1/p0-1)<1e-12 and abs(float(a["NES"])/float(b["NES"])-1)<1e-12: continue
        n+=1
        lr=abs(math.log2(p1/p0))
        if b["log2err"] in("","NA"): na+=1
        else: worst=max(worst,lr/float(b["log2err"]))
        nes.append(abs(float(a["NES"])/float(b["NES"])-1))
        s=abs(math.log2(float(r42[k]["pval"])/p0)); seedratio.append(s)
        pa,pr=float(a["p.adjust"]),float(b["padj"])
        if (pa<=0.05)!=(pr<=0.05): near.append(k)
    seedratio.sort()
    print(coll,"sets differing",n,"| max |log2 p ratio|/log2err",round(worst,4),"| log2err NA",na,"| max NES rel",max(nes,default=0),
          "| seed1 vs seed42 median |log2 p ratio| for same sets",round(seedratio[len(seedratio)//2],4),"| significance flips",len(near))
