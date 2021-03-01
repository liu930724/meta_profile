### 1. meta_profile简介及使用
meta_profile v1.0，集合trimming，remove host，metaphlan3的宏基因组快速profile流程。

更全面的宏基因组分析流程见朱杰的metapi。

#### 1.2 安装
+ 使用conda安装相关依赖。conda相关教程见 https://biogit.cn/TianLiu/meta_course/wikis/c1.p3.Conda.
```shell
conda env create -n meta_profile -f rules/env.yaml
conda activate meta_profile
```

+ 如果你在BGI的集群上，可直接加载工作环境。
```shell
# solution1 : use full env path
conda activate /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_profile

# solution2 : use ~/.conda/environments.txt
echo "/ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_profile" >> ~/.conda/environments.txt
conda activate meta_profile
```

#### 1.3 数据库
+ 宿主的bowtie2索引文件
BGI集群上默认为人的GRCh38.p13参考基因组。该基因组已经包括不同地区人的序列亚型。
```
#https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/

bowtie2-build --threads 4 GCF_000001405.39_GRCh38.p13_genomic.fna.gz hg38
```
索引建好后，在config.yaml里更新bowtie2_index的路径。

+ metaphlan 数据库
BGI集群上默认为mpa_v30_CHOCOPhlAn_201901版本数据库。
```
#mpa_v30_CHOCOPhlAn_201901
mkdir metaphlan_database && cd metaphlan_database
wget http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901.tar .
wget http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901.md5 .
wget http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases/mpa_latest .

metaphlan install --bowtie2db .
```
安装完毕后，在config.yaml里更新metaphlan3的bowtie2db和index信息。

#### 1.3 配置文件

**sample.txt** : 待处理的样本路径，包括了两个测试样本。

**rules/profile.smk** : 主程序，记录有流程的执行规则 。

**config.yaml** : 配置文件，可在此设定各步骤使用的CPU等相关参数，host index以及结果保存的路径。

**默认宿主为人(hg38.p13)，若是其他的宿主，请在config.yaml中更换index**

**cluster.yaml** : 投递任务的配置文件，可在此设置项目编号，任务队列，任务所需资源等参数。

**work.sh** : 投递任务的脚本。


#### 1.4 运行

软件默认使用hg38对人体共生微生物进行过滤，若宿主为其他类型，请在config.yaml文件中替换基因组的bowtie2索引路径。

将需要运行样本的ID和路径追加到sample.txt中，注意分隔符为\t，可从excle里直接把数据粘贴过来。

```shell
### dry run, 检测下机文件路径是否存在：
snakemake --snakefile rules/profile.smk -n

### 小样本测试节点运行，确保流程无误
snakemake --snakefile rules/profile.smk --core 24 2> smk.log &

### 大规模投任务，记得更改cluster.yaml文件中的配置信息
nohup sh work.sh &
```

**1.5 结果路径**

过滤掉宿主的样本fq结果目录 ：**1.assay/02.rmhost/**

样本metaphlan3的结果目录 : **1.assay/03.profile/metaphlan3/**

\*vir.profile为包括病毒的profile结果；\*unknown.profile在前者基础上加了UNKNONW。

过滤统计结果 : **2.result/filter_summary.txt**

合并后的profile结果: **2.result/metaphlan3.profile.merge.txt**



### 2. 工作流程详解

#### 2.1 Trimming

使用[fastp](https://github.com/OpenGene/fastp)对下机数据进行修剪。

早期BGISEQ测序平台下机reads中间部分会出现零散的质量较差的碱基，cOMG流程中的OAs1修剪方案专门针对此情况进行了优化。OAs1会从头识别整段read，若中间出现低质量碱基，则会修剪掉后半段碱基。该策略会尽可能保留多的reads，但会损失大量的base(>10%)。Reads更多有利于profiling，但损失大量的base不利于组装。

现在BGISEQ测序平台的测序质量已经大大提高，总体Q30%已经到了90%。并且SPAdes组装器会对reads先纠错后组装。因此本流程过滤采用fastp的默认模式，仅对reads两端的低质量序列进行修剪。BGISEQ测序平台下机默认去除了接头序列，无需再对接头进行处理。

#### 2.2 Remove host component

cOMG中使用SOAP比对到参考基因组用来去宿主，soap使用的seed为30，且需要保存中间结果.soap。计算速度慢，空间占用大，不推荐使用。

在实际使用中发现，seed值会大大影响宿主率的结果。在高宿主率样本中，seed 30会遗漏很多人源的reads。当然，seed过低也会将微生物的reads误认为是人源的，在部分粪便样本的stLFR数据中发现，BWA默认seed 19会导致10%以上的宿主率，大量reads的match只有19bp，而将seed改到23宿主率则会回归正常的1-3%水平。关于过滤参数的讨论，可以进一步阅读:
[Aligner optimization increases accuracy and decreases compute times in multi-species sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5643015/)

本流程使用bowtie2的--very-sensitive模式比对到参考基因组上，然后通过管道传给samtools -f 12提取PE reads都没有比对上参考基因组的结果。--very-sensitive比对模式在准确度和性能上有较好的平衡。

顺带一提，若组装样本中含有大量真菌，需要考虑对宿主基因组的高度同源区进行屏蔽，避免与宿主同源的真菌序列(eg. 核糖体)在过滤宿主中丢失，影响真菌的组装。该策略详见[BBmap](http://seqanswers.com/forums/showthread.php?t=42552)，该策略不包含在本流程中。

#### 2.3 MetaPhlAn3

MetaPhlAn基于marker基因，可以快速生成样本的profile，并且计算资源消耗低，现已推出第三版，新版本的bowtie2中间结果文件与MetaPhlAn2并不兼容。

MetaPhlAn relies on unique clade-specific marker genes identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, ~3,500 viral, and ~110 eukaryotic).

**What's new in version 3**

+ New MetaPhlAn marker genes extracted with a newer version of ChocoPhlAn based on UniRef

+ Estimation of metagenome composed by unknown microbes with parameter `--unknown_estimation`

+ Automatic retrieval and installation of the latest MetaPhlAn database with parameter `--index latest`

+ Virus profiling with `--add_viruses`

- Calculation of metagenome size for improved estimation of reads mapped to a given clade

- Inclusion of NCBI taxonomy ID in the ouput file

- CAMI (Taxonomic) Profiling Output Format included

- Removal of reads with low MAPQ values

由于数据库中病毒的clade基因较少，metaphlan得到的病毒的profile并不准确，仅供参考。病毒的profile现在暂无公认流程。

### 3. Troubleshooting
#### 3.1 fastp缺少zlib文件。
``` shell
./fastp: /lib64/libz.so.1: version `ZLIB_1.2.3.5' not found (required by ./fastp)
```
解决方法: 自行下载zlib，或在～/.bashrc添加我的zlib。
```shell
export LD_LIBRARY_PATH=/share/app/gcc-5.2.0/lib:/share/app/gcc-5.2.0/lib64:$LD_LIBRARY_PATH
ZLIB_HOME=/ldfssz1/ST_META/P18Z10200N0127_VC/tianliu/bin/lib/zlib-1.2.11
export LD_LIBRARY_PATH=$ZLIB_HOME:$LD_LIBRARY_PATH
```

#### 3.2 MetaPhlAn3无访问database的权限
```shell
ERROR: The directory is not writeable: /ldfssz1/ST_META/SHR/opt/MetaPhlAn/metaphlan/metaphlan_databases. Please modify the permissions.
```
解决方法:
database需要写入权限, 将拷贝metaphlan3到工作目录并在config.yaml里替换路径。
```cp -r /ldfssz1/ST_META/SHR/opt/MetaPhlAn/metaphlan /your/tools/path```


#### 3.3 下机数据路径不存在
snakemake在运行前会检查所有数据路径是否存在。若不存在则会报缺失数据的错误。
```shell
MissingInputException in line 4 of rules/profile.smk:
Missing input files for rule filter:
0.data/test3.fq.gz
```
请依次检查sample.txt是否是按制表符分隔，样本的路径是否完整，路径是否更新过(长周期项目往往会转移数据)。

可以通过check_PE_reads_exist.py检查哪些下机数据路径缺失:
```shell
python rules/check_PE_reads_exist.py sample.txt
```
sample.txt.exist为数据路径正常样本；sample.txt.noexist为数据路径缺失样本。

#### 3.4 同一个样本有多个文库的下机数据
将同一样本的多个下机数据先合并，再跑流程。提供merge_multi_fq.py脚本合并多个下机数据的样本。

将所有多个下机数据的样本按sample.txt的格式整理为sample_dup.txt文件，该文件中不包含单样本单文库的样本。merge_multi_fq.py脚本将按id合并fq文件，合并后的下机数据保存在1.assay/00.tmp中。

如test样本测了两个文库，希望将这两个文库合并后再跑程序，则sample_dup.txt的内容应为:
id	fq1	fq2
test	0.data/test1_1.fq.gz	0.data/test1_2.fq.gz
test	0.data/test2_1.fq.gz	0.data/test2_2.fq.gz
```shell
python rules/merge_multi_fq.py sample_dup.txt sample_dup_merge.txt
```
将rules/profile.smk的21行中的sample.txt替换成sample_dup_merge.txt后，再运行流程即可。

#### 3.5 追加样本后显示Nothing to be done

流程每次检查2.result/filter_summary.txt，和metaphlan3.profile.merge.txt来判断是否执行完毕。因此若在上一批样本跑完了全流程后，再追加样本时由于上两个文件已经生成了，就会显示Nothing to be done。

若已经生成了结果文件，则删去结果文件夹即可。``` rm -r 2.result```。

### 4. 更新计划
1. 增加BWA去除宿主的选项。
BWA在宏基因组比对过程中，会出现大量仅有局部比对上的soft clipping比对结果，该比对结果与reference的identity往往较低。此种策略在粪便，舌苔等低宿主率的样本中会将大量微生物的Reads丢掉。因此BWA适用于单基因组比对而不适合宏基因组比对。
而生殖道等样本的宿主率高达98%以上，可近似看成单基因组的样本。为了尽可能避免人源reads对后续分析造成干扰，生殖道，血液等样本即可使用BWA比对，同时也可以与GATK等变异检测流程对接。

2. 增加功能profile。
HUMAnN，存储的中间结果过多，大量消耗盘阵，还需优化。

### 5. 进一步阅读
1. 宏基因组快速处理流程
Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3
2. profile工具评估
Benchmarking Metagenomics Tools for Taxonomic Classification
3. 种子序列长度与灵敏度之间的关系
Aligner optimization increases accuracy and decreases compute times in multi-species sequence data
4. 人和微生物序列的混杂对后续分析的干扰
Contaminating DNA in human saliva alters the detection of variants from whole genome sequencing
