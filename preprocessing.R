library(minfi)
library(sva)
library(limma)
#1、造一个SampleSheet.csv,里面有三列，分别是Sample_Name(样本名)Slide(片编号)Array(阵列编号),
#eg:GSM2405657(样本名)_2002?7210073(片编号)_R05C01(阵列编号)_Red
#2、把SampleSheet.csv和数据放在一个目录下
#3、导入数据
#basedir是数据的目录
basedir='/home/chenyx/dnameth/gse90496/raw_data'
library(readr)
batch_file <- read_csv("/home/chenyx/dnameth/gse90496/batch_file.csv")
#targets会把SampleSheet.csv读进来并添加一列，写有每个样本的路径
targets <- read.metharray.sheet(basedir)
#根据targets表把数据读进来，现在的数据只是探针水平，没有比对到基因组
single_sample<-function(i){
RGset<-read.metharray.exp(targets = targets[i,])
#背景校正
RGset<-bgcorrect.illumina(RGset)
intensity<-10000
Green<-getGreen(RGset)
Red<-getRed(RGset)    
Green_control_probe_adress<-getControlAddress(RGset,controlType = c("NORM_C","NORM_G"))
Red_control_probe_adress<-getControlAddress(RGset,controlType = c("NORM_A","NORM_T"))
Green_avg<-colMeans(Green[Green_control_probe_adress, , drop = FALSE])
Red_avg<-colMeans(Red[Red_control_probe_adress, , drop = FALSE])
Green_unit<-intensity/Green_avg
Red_unit<-intensity/Red_avg
Green<-sweep(Green,2,FUN = "*",Green_unit)
Red<-sweep(Red,2,FUN = "*",Red_unit)
RGset@assays@data@listData[["Green"]]<-Green
RGset@assays@data@listData[["Red"]]<-Red
#5、把光信号转变成探针信号
MsetEx<-preprocessRaw(RGset)
return(MsetEx)}
MsetEx<-single_sample(1)
for(i in 2:2801){
  MsetEx<-cbind(MsetEx,single_sample(i))
}
#6、文献提到limma包的removeBatchEffect函数去除冷冻组织和石蜡包埋之间的批次效应，
#但是我们的样本可能没有这种情况，但是由于被测样本在不同的芯片，我们可以对这种情况做类似的处理，
#把片的差异类比为冷冻和石蜡包埋的差异

#先把矩阵及列名提出来
col_names<-colnames(MsetEx@assays@data@listData[["Meth"]])
Meth<-MsetEx@assays@data@listData[["Meth"]]
Unmeth<-MsetEx@assays@data@listData[["Unmeth"]]
#再把片编号提出来当成批次编号，由于示例数据的样本文件名eg:"5723646052_R02C02"这种形式，
#与文献的数据名格式不一致，具体做的时候需要纠正一下
batch<-batch_file$`!Sample_characteristics_ch1_1`
#log2变换、去批次
removed_batch_Meth<- 2^removeBatchEffect(log2(Meth+1),batch)
removed_batch_Unmeth<-2^removeBatchEffect(log2(Unmeth+1),batch)
MsetEx@assays@data@listData[["Meth"]]<-removed_batch_Meth
rm(Meth)
gc()
MsetEx@assays@data@listData[["Unmeth"]]<-removed_batch_Unmeth
rm(Unmeth)
gc()
print('55th out')
#5、过滤探针
#获取注释信息
annotation<-getAnnotation(MsetEx, what = "everything")
#去掉性染色体探针
sex_probe<-rownames(annotation)[annotation$chr %in% c("chrX", "chrY")]
keep<-!(featureNames(MsetEx) %in% sex_probe)
MsetEx<-MsetEx[keep, ]
gc()
#将探针比对到基因组
GMsetEx<-mapToGenome(MsetEx, mergeManifest = FALSE)
rm(MsetEx)
gc()
print('67th out')
#去掉snp位点
GMsetEx<-dropLociWithSnps(GMsetEx, snps = c("CpG"), maf = 0, snpAnno = NULL)
gc()
library(readr)
GPL450k<- read_delim("/home/chenyx/dnameth/450k.txt","\t", escape_double = FALSE, trim_ws = TRUE)
#去掉不能唯一比对到基因组的探针
Ununique_probe_index<-c()
for(i in 1:dim(GPL450k)[1])
{temp_vector<-unique(strsplit(GPL450k$UCSC_RefGene_Name[i],";")[[1]])
if(length(temp_vector)>1)
{Ununique_probe_index<-c(Ununique_probe_index,i)}
print(i)}
Ununique_probe<-GPL450k$ID[Ununique_probe_index]
temp_vector<-na.omit(match(Ununique_probe,rownames(GMsetEx)))
GMsetEx<-GMsetEx[-temp_vector,]
gc()
rm(GPL450k)
print('75th out')
#去掉EPIC中没有的探针
GPL850k <- read_csv("/home/chenyx/dnameth/850k.csv")
Not_in_EPIC_probe<-setdiff(rownames(GMsetEx),GPL850k$IlmnID)
temp_vector<-na.omit(match(Not_in_EPIC_probe,rownames(GMsetEx)))
GMsetEx<-GMsetEx[-temp_vector,]
gc()
#6、获取beta值
beta<-getBeta(GMsetEx,offset=100)
#这一段为sva，判断信息是否冗余，我认为是多余的
# mod = model.matrix(~as.factor(batch))
# mod0 = model.matrix(~1,data=as.data.frame(batch))
# n.sv = num.sv(beta,mod,method="leek")
# svobj = sva(beta,mod,mod0,n.sv=n.sv)
write.csv(beta,'/home/chenyx/dnameth/remake/beta/gse90496.csv')
#代码重复跑之前需要把beta.csv移除，否则读取sample_sheet的函数无法判断哪个是样本信息表