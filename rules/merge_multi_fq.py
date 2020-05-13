#############################################
### merge fq.gz
### tianliu@genomics.cn
### 2019-11-27
#############################################

# python merge_multi_fq.py sample_dup.txt sample_dup_merge.txt

import sys
import os

infile = sys.argv[1]
outfile = sys.argv[2]

if not os.path.exists("1.assay/00.tmp"):
    os.system("mkdir -p 1.assay/00.tmp")

f_o = open(outfile, "w")
f_o.write("id\tfq1\tfq2\n")

path_set = set()
for i in open(infile).readlines()[1:]:
    sample, r1, r2 = i.strip().split("\t")
    path_set.add("%s\t1.assay/00.tmp/%s.1.fq.gz\t1.assay/00.tmp/%s.2.fq.gz\n" %(sample, sample, sample))
    r1_rmhost = "1.assay/02.rmhost/%s.rmhost.1.fq.gz" %(sample)
    if os.path.exists(r1_rmhost):
        continue
    else:
        os.system("cat %s >> 1.assay/00.tmp/%s.1.fq.gz" %(r1, sample))
        os.system("cat %s >> 1.assay/00.tmp/%s.2.fq.gz" %(r2, sample))

for path in path_set:
    f_o.write(path)

f_o.close()
