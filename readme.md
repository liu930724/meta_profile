### 1. meta_profile简介及使用

meta_profile v1.0，集合trimming，remove host，metaphlan3的宏基因组快速profile流程。

#### 1.1 依赖软件:

snakemake (v5.14.0)

python3 (v3.6.10),  pandas

bowtie2 (v2.3.5.1)

samtools (v1.9)

metaphlan3 (v3.0)

seqkit (v0.12.1)

#### 1.2 安装

如果你在BGI的集群上:

```shell
#1.加载朱杰的工作环境到～/.bashrc
export PATH="/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin:$PATH"

#2.拷贝流程到工作目录
cp -r /ldfssz1/ST_META/share/User/tianliu/pipline/meta_profile/* /your/path

#3.mp3暂时有bug, database需要写入权限, 将拷贝metaphlan3到工作目录并在config.yaml里替换路径。
cp -r /ldfssz1/ST_META/SHR/opt/MetaPhlAn/metaphlan /your/tools/path

```

如果不是的话，推荐使用conda安装管理上述所有的软件。

#### 1.3 配置文件

**sample.txt** : 待处理的样本路径，包括了两个测试样本。

**rules/profile.smk** : 主程序，记录有流程的执行规则 。

**config.yaml** : 配置文件，可在此设定各步骤使用的CPU等相关参数，host index以及结果保存的路径。

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

现在BGISEQ测序平台的测序质量已经大大提高，并且SPAdes组装器会对reads先纠错后组装。因此本流程过滤采用fastp的默认模式，仅对reads两端的低质量序列进行修剪。BGISEQ测序平台下机默认去除了街头序列。

#### 2.2 Remove host component

cOMG中使用SOAP比对到参考基因组用来去宿主，soap使用的seed为30，且需要保存中间结果.soap。计算速度慢，空间占用大，不推荐使用。

在实际使用中发现，seed值会大大影响宿主率的结果。在高宿主率样本中，seed 30会遗漏很多人源的reads。当然，seed过低也会将微生物的reads误认为是人源的，在部分粪便样本的stLFR数据中发现，BWA默认seed 19会导致10%以上的宿主率，大量reads的match只有19bp，而将seed改到23宿主率则会回归正常的1-3%水平。关于过滤参数的讨论，可以进一步阅读:
[Aligner optimization increases accuracy and decreases compute times in multi-species sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5643015/)

本流程使用bowtie2(BWA/bowtie2按个人喜好即可)的--very-sensitive模式比对到参考基因组上，然后通过管道传给samtools -f 12提取PE reads都没有比对上参考基因组的结果。--very-sensitive比对模式在准确度和性能上有较好的平衡。

顺带一提，若组装样本中含有大量真菌，需要考虑对宿主基因组的高度同源区进行屏蔽，避免与宿主同源的真菌序列(eg. 核糖体)在过滤宿主中丢失。该策略详见[BBmap](http://seqanswers.com/forums/showthread.php?t=42552)，不包含在本流程中。

#### 2.3 MetaPhlAn3

MetaPhlAn基于marker基因，可以快速生成样本的profile，并且计算资源消耗低，现已推出第三版，新版本的bowtie2中间结果文件与MetaPhlAn2并不兼容。

**What's new in version 3**

+ New MetaPhlAn marker genes extracted with a newer version of ChocoPhlAn based on UniRef

+ Estimation of metagenome composed by unknown microbes with parameter `--unknown_estimation`

+ Automatic retrieval and installation of the latest MetaPhlAn database with parameter `--index latest`

+ Virus profiling with `--add_viruses`

- Calculation of metagenome size for improved estimation of reads mapped to a given clade

- Inclusion of NCBI taxonomy ID in the ouput file

- CAMI (Taxonomic) Profiling Output Format included

- Removal of reads with low MAPQ values

  

使用时需注意将整个文件夹拷贝到自己的工作目录下，否则会出现database无法写入的错误。

MetaPhlAn3默认加入了UNKNOWN的结果，在1.assay/03.profile/metaphlan3/下同时保存着有无UNKNOWN的结果，merge_profile结果是按有UNKNOWN结果合并的，大家按需取用。

更多profile的方法见综述[Benchmarking Metagenomics Tools for Taxonomic Classification](https://www.cell.com/cell/fulltext/S0092-8674(19)30775-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419307755%3Fshowall%3Dtrue)。

### 3. Dedug
1. fastp缺少zlib文件。
``` shell
./fastp: /lib64/libz.so.1: version `ZLIB_1.2.3.5' not found (required by ./fastp)
```
解决方法: 自行下载zlib，或在～/.bashrc添加我的zlib。
```shell
export LD_LIBRARY_PATH=/share/app/gcc-5.2.0/lib:/share/app/gcc-5.2.0/lib64:$LD_LIBRARY_PATH
ZLIB_HOME=/ldfssz1/ST_META/P18Z10200N0127_VC/tianliu/bin/lib/zlib-1.2.11
export LD_LIBRARY_PATH=$ZLIB_HOME:$LD_LIBRARY_PATH
```

2. MetaPhlAn3无访问database的权限
```shell
ERROR: The directory is not writeable: /ldfssz1/ST_META/SHR/opt/MetaPhlAn/metaphlan/metaphlan_databases. Please modify the permissions.
```
解决方法:
将MetaPhlAn3拷贝到自己的工作目录下。


