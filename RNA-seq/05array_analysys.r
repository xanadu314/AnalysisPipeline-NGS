#仅针对芯片数据
#下游分析同RNA seq
#部分jimmy包需要fq下载
library(GEOmirror)
## Setting options('download.file.method.GEOquery'='auto')
## Setting options('GEOquery.inmemory.gpl'=FALSE)
## GEOmirror v 0.1.0  welcome to use GEOmirror!
## If any suggestion please feel free to email to jmzeng1314@163.com!
## If you use GEOmirror in published research, please acknowledgements:
## We thank Dr.Jianming Zeng(University of Macau), and all the members of his bioinformatics team, biotrainee, for generously sharing their experience and codes.
library(GEOquery)
## Loading required package: Biobase
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which, which.max, which.min
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
eSet=geoChina("GSE33532")
## you can also use getGEO from GEOquery, by 
## getGEO('GSE33532', destdir=".", AnnotGPL = F, getGPL = F)
head(eSet)
## $GSE33532_series_matrix.txt.gz
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 25906 features, 100 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: GSM835268 GSM835269 ... GSM835367 (100 total)
##   varLabels: title geo_accession ... tumor_stage:ch1 (41 total)
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation: GPL570
eSet=eSet[[1]]
head(eSet)
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 1 features, 100 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: GSM835268 GSM835269 ... GSM835367 (100 total)
##   varLabels: title geo_accession ... tumor_stage:ch1 (41 total)
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation: GPL570
probes_expr <- exprs(eSet) #获取探针表达矩阵
head(probes_expr)
##           GSM835268 GSM835269 GSM835270 GSM835271 GSM835272 GSM835273 GSM835274
## 1007_s_at  9.872989 10.434562 10.093931  9.830508  9.571424  9.904023  9.736608
## 1053_at    6.648595  6.773601  6.643051  6.854348  5.558082  7.364951  7.313815
## 121_at     7.185752  6.980180  7.153527  7.199128  6.664586  7.811820  7.762098
## 1255_g_at  4.117641  4.318899  3.796970  4.030910  3.716171  4.319518  4.400075
## 1294_at    7.529764  7.387681  7.500212  7.671593  8.550646  7.885696  7.801964
## 1316_at    5.255940  4.946495  5.527133  5.217120  6.163936  5.177423  5.282608
##           GSM835275 GSM835276 GSM835277 GSM835278 GSM835279 GSM835280 GSM835281
## 1007_s_at  9.787895  9.719861  9.829164 10.763830 10.640284 10.742810 10.659695
## 1053_at    7.310051  6.958884  6.142737  7.243693  7.246412  7.315477  7.120801
## 121_at     7.655468  7.744584  7.563815  6.918900  6.928518  6.820209  6.912232
## 1255_g_at  4.493545  5.010725  3.933282  3.905165  3.796434  3.840558  3.966955
## 1294_at    7.956659  7.966873  8.351416  7.855689  7.805189  7.763224  7.878375
## 1316_at    5.310621  5.517881  6.152662  5.116082  5.740664  5.645054  6.025233
##           GSM835282 GSM835283 GSM835284 GSM835285 GSM835286 GSM835287 GSM835288
## 1007_s_at  9.828964  9.918608  9.888298 10.268530  9.801217  9.404580 10.227295
## 1053_at    5.815041  5.834042  5.899908  5.932526  5.695452  6.072709  7.299339
## 121_at     6.557610  7.468079  7.418716  7.387166  7.420231  7.080029  7.398838
## 1255_g_at  3.665432  4.749513  4.407903  5.050961  4.207951  3.660254  4.122986
## 1294_at    8.686691  7.335405  7.228783  7.197245  7.495444  8.661079  6.179481
## 1316_at    6.184049  5.610292  5.689444  5.678307  5.835542  6.114158  5.121633
##           GSM835289 GSM835290 GSM835291 GSM835292 GSM835293 GSM835294 GSM835295
## 1007_s_at 10.254358 10.294003 10.087927  9.583512  9.979222 10.348294 10.047674
## 1053_at    7.202675  7.052240  6.857298  5.929094  6.965075  7.273838  6.800071
## 121_at     7.435463  7.463110  7.411235  7.163564  7.081120  7.023661  7.129218
## 1255_g_at  4.048515  4.174773  3.773285  3.881541  4.621946  3.958328  4.240645
## 1294_at    6.461623  6.338560  7.078877  8.044031  7.083999  7.051791  7.207025
## 1316_at    5.106532  4.979522  4.908048  5.449323  5.654336  5.700402  5.538671
##           GSM835296 GSM835297 GSM835298 GSM835299 GSM835300 GSM835301 GSM835302
## 1007_s_at 10.283157  9.561983 10.002821  9.958619  9.989986  9.962364  9.382145
## 1053_at    6.957979  5.927921  6.541091  6.225832  6.724558  6.035666  5.997237
## 121_at     7.086807  7.103749  7.013140  6.984259  6.920129  6.871263  6.875759
## 1255_g_at  4.173323  3.867787  4.111473  3.950358  4.221346  3.910653  3.702536
## 1294_at    7.392710  8.200965  7.011403  7.586410  7.390797  7.526705  8.417145
## 1316_at    5.580544  5.835348  5.505628  5.212945  5.711949  5.379999  6.159012
##           GSM835303 GSM835304 GSM835305 GSM835306 GSM835307 GSM835308 GSM835309
## 1007_s_at 10.148137 10.379708 10.398468 10.529892  9.601174  9.843413  9.913944
## 1053_at    6.420194  6.437205  6.236970  6.090034  5.893326  6.399296  6.479475
## 121_at     7.502800  7.380074  7.485525  7.640895  7.959765  7.175029  7.288277
## 1255_g_at  5.735638  4.715692  5.491776  4.101910  3.890514  3.809013  3.785137
## 1294_at    7.143606  6.782291  7.151813  6.634158  8.411357  6.785528  6.873914
## 1316_at    5.477775  5.401552  5.467527  5.365946  5.339052  5.517225  5.831077
##           GSM835310 GSM835311 GSM835312 GSM835313 GSM835314 GSM835315 GSM835316
## 1007_s_at  9.893826 10.009412  9.534127 10.677097  9.958462 10.768931 10.704244
## 1053_at    6.321759  6.333644  6.008040  7.711346  6.680224  7.582875  7.264642
## 121_at     7.288135  7.339073  6.954613  7.274443  7.281953  7.282770  7.332462
## 1255_g_at  3.745945  3.830927  3.682324  4.779105  3.982447  4.728782  4.819307
## 1294_at    7.194554  7.136078  8.553064  7.563668  8.094164  7.554775  8.003657
## 1316_at    5.462127  5.395751  6.273538  5.777604  5.789652  5.698221  5.756313
##           GSM835317 GSM835318 GSM835319 GSM835320 GSM835321 GSM835322 GSM835323
## 1007_s_at  9.407688  9.085575  9.059848  9.019359  9.388711  9.700179 10.720057
## 1053_at    6.023851  6.240819  6.186271  6.312622  6.257651  5.907675  6.344311
## 121_at     7.104356  7.455106  7.610141  7.429427  7.321183  7.275903  7.474187
## 1255_g_at  3.591072  4.780106  5.130357  5.510398  4.062747  3.781078  4.659060
## 1294_at    8.536748  8.061009  7.877922  7.761523  7.774723  8.480853  7.700691
## 1316_at    5.709324  5.569938  5.820410  5.955582  5.828894  6.157210  5.791349
##           GSM835324 GSM835325 GSM835326 GSM835327 GSM835328 GSM835329 GSM835330
## 1007_s_at 10.480827 11.033909 10.832540  9.320314 10.764057 10.505017 10.333511
## 1053_at    6.165422  6.372226  6.600278  5.845705  6.403352  6.601547  6.214928
## 121_at     7.354236  7.576010  7.527176  6.753347  7.674834  7.592396  7.759395
## 1255_g_at  4.167277  4.495923  4.224582  3.551185  5.936536  6.750587  5.057516
## 1294_at    7.871611  7.231740  7.683526  8.583053  7.393855  7.401539  7.555430
## 1316_at    5.735119  5.793075  5.807199  5.864810  5.417522  5.465534  5.458177
##           GSM835331 GSM835332 GSM835333 GSM835334 GSM835335 GSM835336 GSM835337
## 1007_s_at 10.313240  9.051339 10.368356 10.337808 10.516741 10.316717  9.614298
## 1053_at    6.166719  6.061229  7.105526  7.035240  7.329425  7.128353  5.970752
## 121_at     7.719251  7.397093  7.001350  6.897826  6.919596  6.991412  6.638600
## 1255_g_at  6.146795  3.953080  5.693654  5.237110  4.452377  4.770525  3.746660
## 1294_at    7.701038  8.145234  7.046271  6.640916  6.688493  6.957180  8.495759
## 1316_at    5.373341  5.992198  5.527030  5.772854  5.614471  5.577908  6.071189
##           GSM835338 GSM835339 GSM835340 GSM835341 GSM835342 GSM835343 GSM835344
## 1007_s_at 10.059178  9.884711  9.923491  9.872638  9.675588  9.584881  9.401099
## 1053_at    5.474820  5.522895  5.449596  5.464653  5.960509  5.871726  6.257143
## 121_at     7.209450  7.231131  7.177149  7.060270  6.952749  7.109454  7.091216
## 1255_g_at  3.728187  3.907912  3.877566  3.754272  3.640967  3.989741  4.141127
## 1294_at    8.107388  8.024072  8.187213  8.056627  8.793582  7.612829  7.564685
## 1316_at    5.645931  5.764033  5.921983  5.781348  6.420803  5.598685  5.942978
##           GSM835345 GSM835346 GSM835347 GSM835348 GSM835349 GSM835350 GSM835351
## 1007_s_at  9.698175  9.844314  9.678102 10.335628 10.718721 10.816063 10.858924
## 1053_at    5.564758  6.164904  5.960457  5.351954  5.370681  5.225038  5.118605
## 121_at     7.374186  7.302561  7.132505  7.093121  6.981173  7.148704  7.132393
## 1255_g_at  3.909426  3.843012  3.766365  3.816334  3.767770  3.625902  3.832739
## 1294_at    7.387527  7.310905  8.535733  7.642486  7.215705  7.003715  6.927074
## 1316_at    5.540392  5.557051  5.794943  5.905031  6.177188  5.965282  6.021482
##           GSM835352 GSM835353 GSM835354 GSM835355 GSM835356 GSM835357 GSM835358
## 1007_s_at  9.852283 11.334220 10.895149 10.482554 11.012359  9.295507 10.959623
## 1053_at    5.670191  5.955507  6.248399  6.255390  5.616761  5.815319  7.295961
## 121_at     6.913990  7.380207  7.351131  7.521451  7.427777  6.694294  6.736662
## 1255_g_at  3.723612  3.902949  3.800862  3.925481  4.136847  3.512346  3.959987
## 1294_at    8.586559  7.248416  7.384903  7.358016  7.409237  8.812530  6.633219
## 1316_at    6.028830  5.531689  5.723834  5.409114  5.649910  6.064351  5.119413
##           GSM835359 GSM835360 GSM835361 GSM835362 GSM835363 GSM835364 GSM835365
## 1007_s_at 10.936236 10.848241 11.136892  9.313465  9.107971  9.307671  9.220561
## 1053_at    7.260342  7.652184  7.067730  5.802689  6.495575  6.670689  6.359590
## 121_at     6.920597  6.762258  6.857436  6.615234  7.040572  7.004058  7.095032
## 1255_g_at  3.896977  3.899963  4.590393  3.778235  4.518251  4.852265  4.327886
## 1294_at    6.849581  6.672233  6.425828  8.352538  7.416545  7.709097  7.549512
## 1316_at    5.254927  5.180164  5.325834  5.882690  5.725747  5.614676  5.642507
##           GSM835366 GSM835367
## 1007_s_at  9.386351  9.575864
## 1053_at    6.307005  5.768602
## 121_at     7.208753  7.102090
## 1255_g_at  4.497089  3.682845
## 1294_at    7.618242  8.591617
## 1316_at    5.993895  6.116266
dim(probes_expr)
## [1] 25906   100
head(probes_expr[,1:4])
##           GSM835268 GSM835269 GSM835270 GSM835271
## 1007_s_at  9.872989 10.434562 10.093931  9.830508
## 1053_at    6.648595  6.773601  6.643051  6.854348
## 121_at     7.185752  6.980180  7.153527  7.199128
## 1255_g_at  4.117641  4.318899  3.796970  4.030910
## 1294_at    7.529764  7.387681  7.500212  7.671593
## 1316_at    5.255940  4.946495  5.527133  5.217120
boxplot(probes_expr,las=2) #查看有无异常值
phenoDat <- pData(eSet) #获取分组信息
head(phenoDat[,1:4])
##                                     title geo_accession                status
## GSM835268  patient 02, tumor sub-sample A     GSM835268 Public on Sep 23 2014
## GSM835269  patient 02, tumor sub-sample B     GSM835269 Public on Sep 23 2014
## GSM835270  patient 02, tumor sub-sample C     GSM835270 Public on Sep 23 2014
## GSM835271  patient 02, tumor sub-sample D     GSM835271 Public on Sep 23 2014
## GSM835272 patient 02, matched normal lung     GSM835272 Public on Sep 23 2014
## GSM835273  patient 07, tumor sub-sample A     GSM835273 Public on Sep 23 2014
##           submission_date
## GSM835268     Nov 17 2011
## GSM835269     Nov 17 2011
## GSM835270     Nov 17 2011
## GSM835271     Nov 17 2011
## GSM835272     Nov 17 2011
## GSM835273     Nov 17 2011
tumor=rownames(phenoDat[grepl('tumor',as.character(phenoDat$title)),])
normal=rownames(phenoDat[grepl('normal',as.character(phenoDat$title)),])
dat=probes_expr
dat=dat[,c(tumor,normal)] #排除一些不要的分组
group_list=c(rep('tumor',length(tumor)),rep('normal',length(normal))) #这个具体情况具体写
table(group_list)
## group_list
## normal  tumor 
##     20     80
dim(dat)
## [1] 25906   100
dat[1:4,1:4] 
##           GSM835268 GSM835269 GSM835270 GSM835271
## 1007_s_at  9.872989 10.434562 10.093931  9.830508
## 1053_at    6.648595  6.773601  6.643051  6.854348
## 121_at     7.185752  6.980180  7.153527  7.199128
## 1255_g_at  4.117641  4.318899  3.796970  4.030910
dat <- as.data.frame(dat)
save(dat,group_list,file = 'step1.1_dat_group.Rdata')

#下载注释文件
gpl=eSet@annotation
library(AnnoProbe)#用了jimmy的包
## AnnoProbe v 0.1.0  welcome to use AnnoProbe!
## If you use AnnoProbe in published research, please acknowledgements:
## We thank Dr.Jianming Zeng(University of Macau), and all the members of his bioinformatics team, biotrainee, for generously sharing their experience and codes.
## 
## Attaching package: 'AnnoProbe'
## The following objects are masked from 'package:GEOmirror':
## 
##     filterEM, geoChina
checkGPL(gpl)
## [1] TRUE
printGPLInfo(gpl)
##    132                                                           
## V1 "GPL570"                                                      
## V2 "[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array"
## V3 "Homo sapiens"
probe2gene=idmap(gpl)
genes_expr <- filterEM(probes_expr,probe2gene) #Jimmy的函数，用于转换芯片表达矩阵到基因表达矩阵
## input expression matrix is 25906 rows(genes or probes) and 100 columns(samples).
## input probe2gene is 41937 rows(genes or probes)
## after remove NA or useless probes for probe2gene, 23739 rows(genes or probes) left
## There are 23739 of 25906 probes can be annotated.
## output expression matrix is 13836 rows(genes or probes) and 100 columns(samples).
head(genes_expr)
##        GSM835268 GSM835269 GSM835270 GSM835271 GSM835272 GSM835273 GSM835274
## ZZZ3    7.673059  7.818656  7.791604  7.559840  7.591636  8.245246  8.190175
## ZZEF1   7.125457  7.361978  7.072040  7.279635  8.396008  7.880050  7.835793
## ZYX     8.356221  8.705160  8.447408  8.535484  9.577757  9.023629  9.019582
## ZYG11B  8.957735  8.456601  8.917016  8.973982  8.916972  8.645889  8.500038
## ZYG11A  4.399679  4.887158  4.652484  4.812378  3.659176  6.444072  6.064733
## ZXDC    6.811726  6.995278  6.755707  6.727731  7.694904  7.145818  6.973559
##        GSM835275 GSM835276 GSM835277 GSM835278 GSM835279 GSM835280 GSM835281
## ZZZ3    8.199950  7.975647  7.281434  8.224006  8.050748  8.081576  7.962779
## ZZEF1   7.895610  7.733515  7.741254  7.569365  7.440815  7.560796  7.651183
## ZYX     9.202452  9.164714  9.715941  8.323588  8.045569  8.342762  8.138669
## ZYG11B  8.494820  8.731307  8.749599  9.183484  9.177427  9.079588  9.108718
## ZYG11A  6.024293  5.798648  4.097012  4.066768  3.986605  4.160754  4.177033
## ZXDC    6.891685  6.842296  7.262184  7.606981  7.594644  7.799865  7.489626
##        GSM835282 GSM835283 GSM835284 GSM835285 GSM835286 GSM835287 GSM835288
## ZZZ3    7.469818  8.297214  8.208226  8.494908  7.763506  7.485167  8.515434
## ZZEF1   8.260506  7.067298  6.972346  7.017068  7.298076  8.086405  7.283029
## ZYX     9.467965  8.591807  8.434475  8.739671  8.576007  9.567128  7.677324
## ZYG11B  8.990652  8.992174  8.973087  8.907123  8.936169  9.117315  8.760104
## ZYG11A  3.858654  3.816472  3.716875  4.274437  3.875055  3.688318  6.808303
## ZXDC    7.792710  6.789350  6.836264  6.952628  6.758253  7.426808  7.252696
##        GSM835289 GSM835290 GSM835291 GSM835292 GSM835293 GSM835294 GSM835295
## ZZZ3    8.360835  8.559062  8.113384  6.786749  8.809625  8.903971  8.934150
## ZZEF1   7.296481  6.985219  7.294649  8.122786  7.862684  7.944346  7.324586
## ZYX     8.082029  7.669410  8.808698  9.740701  9.014335  9.104131  8.350796
## ZYG11B  8.538537  8.784737  8.450530  8.175981  8.473360  8.540240  8.874637
## ZYG11A  6.385760  6.526270  6.214472  3.787070  6.654523  6.746409  5.692094
## ZXDC    7.321210  7.192039  7.131870  7.305082  7.078034  7.276836  7.108951
##        GSM835296 GSM835297 GSM835298 GSM835299 GSM835300 GSM835301 GSM835302
## ZZZ3    8.882163  7.233024  8.029047  7.994040  8.074936  8.019510  7.578642
## ZZEF1   7.962352  8.174577  6.450465  6.759040  7.050566  7.029562  7.952365
## ZYX     8.828801  9.545185  8.504548  8.240996  8.252823  8.274565 10.745388
## ZYG11B  8.528809  8.665179  9.252465  9.286891  9.358776  9.198875  8.975644
## ZYG11A  6.758886  4.094473  5.292898  4.757278  5.239955  4.875352  3.803032
## ZXDC    7.273769  7.329695  6.863476  7.028728  7.202059  7.181948  7.279548
##        GSM835303 GSM835304 GSM835305 GSM835306 GSM835307 GSM835308 GSM835309
## ZZZ3    8.333811  8.585953  8.338957  8.467705  7.596215  8.102888  8.127631
## ZZEF1   6.949255  7.093205  6.900908  7.131913  8.286312  7.070391  7.067336
## ZYX     9.399831  9.069051  9.823166  9.710215 10.118693  7.610372  7.522256
## ZYG11B  8.685566  8.897112  8.677644  8.786226  8.916736  8.248983  8.035623
## ZYG11A  7.426830  7.630317  6.876357  6.924667  3.667200  3.616949  3.759101
## ZXDC    6.855783  7.117934  7.102443  7.379343  7.335582  7.183526  7.335349
##        GSM835310 GSM835311 GSM835312 GSM835313 GSM835314 GSM835315 GSM835316
## ZZZ3    7.738268  7.883992  7.604689  7.082252  7.778259  7.252756  7.420010
## ZZEF1   7.261158  7.163885  8.213308  7.084362  7.525691  7.235120  7.444703
## ZYX     8.325413  8.002064  9.764640  9.009980  9.759429  8.803754  9.495095
## ZYG11B  8.156311  7.867879  8.930634  8.408472  8.713087  8.465734  8.467983
## ZYG11A  3.763976  3.947396  3.690078  3.784395  3.897233  4.306080  3.959607
## ZXDC    7.200090  7.298159  7.428722  6.778446  6.695263  6.606462  6.572803
##        GSM835317 GSM835318 GSM835319 GSM835320 GSM835321 GSM835322 GSM835323
## ZZZ3    7.560903  7.399466  7.548903  7.561303  7.612310  7.439854  7.760924
## ZZEF1   7.602119  7.558293  7.201299  7.159182  7.383451  8.209269  7.260081
## ZYX     9.960238  9.274590  8.668641  8.563979  8.410804  9.843377  9.970337
## ZYG11B  8.851189  8.286721  8.728788  8.827857  8.812867  9.192642  8.399064
## ZYG11A  3.745164  3.783248  4.074701  3.817972  3.949947  3.782239  4.029883
## ZXDC    7.362240  6.848145  6.744611  6.747545  6.901562  7.434379  7.051481
##        GSM835324 GSM835325 GSM835326 GSM835327 GSM835328 GSM835329 GSM835330
## ZZZ3    7.676407  7.957212  7.887984  7.629436  8.269976  8.385118  8.092783
## ZZEF1   7.461855  6.915614  7.400144  8.050700  6.965791  6.998528  6.950403
## ZYX     9.831176  9.786776  9.829971  9.904446  9.348009  9.485039  9.913077
## ZYG11B  8.519453  8.461715  8.380016  9.059824  8.220660  8.363566  8.276595
## ZYG11A  3.927088  4.000322  4.173895  3.556417  4.736709  4.700651  4.766753
## ZXDC    7.145330  7.403411  7.410607  7.467465  6.592886  6.464449  6.382408
##        GSM835331 GSM835332 GSM835333 GSM835334 GSM835335 GSM835336 GSM835337
## ZZZ3    8.010752  8.074552  7.645786  7.820727  7.744477  7.594398  7.649380
## ZZEF1   6.802666  7.854054  7.200594  7.187599  7.149817  7.030431  8.258221
## ZYX    10.230395  9.429178  9.702519  9.097689  9.498158  9.798381  9.795085
## ZYG11B  8.386879  8.659098  8.104119  8.219194  8.058674  7.987895  8.937572
## ZYG11A  4.463988  3.601377  5.830541  6.153377  6.216556  5.943647  3.589968
## ZXDC    6.267071  6.946675  6.408753  6.215946  6.501323  6.268121  7.499155
##        GSM835338 GSM835339 GSM835340 GSM835341 GSM835342 GSM835343 GSM835344
## ZZZ3    7.390777  7.188347  7.363591  7.508110  7.525855  7.463871  7.417055
## ZZEF1   7.491115  7.536219  7.620021  7.654789  8.160394  7.656506  7.531271
## ZYX     8.909242  9.447489  9.237781  9.183886  9.817378  9.037027  8.996077
## ZYG11B  8.538681  8.431099  8.289017  8.548352  8.936119  8.477815  8.308486
## ZYG11A  3.872615  3.882731  3.689939  3.683286  3.594567  4.690363  4.696860
## ZXDC    6.659791  6.540723  6.595842  6.596836  7.482971  7.016846  6.816992
##        GSM835345 GSM835346 GSM835347 GSM835348 GSM835349 GSM835350 GSM835351
## ZZZ3    7.325485  7.138167  7.539135  7.430069  7.304708  7.188530  7.322041
## ZZEF1   7.152219  7.519776  8.082285  7.665253  7.345474  7.560945  7.583971
## ZYX     9.071607  8.672927 10.447565  8.488436  8.615942  8.265500  8.096064
## ZYG11B  8.021958  7.968042  8.850053  8.470150  8.403535  8.494958  8.268716
## ZYG11A  5.215761  5.154171  3.742386  3.618956  3.697355  3.983545  3.868999
## ZXDC    6.848619  7.115784  7.301046  6.665137  6.757307  6.705439  6.852236
##        GSM835352 GSM835353 GSM835354 GSM835355 GSM835356 GSM835357 GSM835358
## ZZZ3    7.435328  7.787018  7.829318  8.074228  7.671735  7.732395  8.253871
## ZZEF1   8.124649  6.965811  7.210269  6.977914  6.510555  8.283491  7.150780
## ZYX    10.172886  8.506672  8.511145  8.466718  9.071480 10.123169  8.499685
## ZYG11B  8.797698  8.133985  8.368218  8.315942  8.086288  8.853830  8.676368
## ZYG11A  3.709823  4.099518  3.872731  4.250887  4.266970  3.737552  5.324720
## ZXDC    7.639972  8.022907  7.573398  7.405813  7.846055  7.676875  7.889629
##        GSM835359 GSM835360 GSM835361 GSM835362 GSM835363 GSM835364 GSM835365
## ZZZ3    8.234161  8.217027  8.416368  7.430859  7.439546  7.479553  7.390038
## ZZEF1   7.292403  7.293753  7.143220  7.979738  6.867336  7.133145  6.968093
## ZYX     8.762709  8.622939  8.884885  9.631262  9.596766  9.418197  9.510312
## ZYG11B  8.550713  8.558573  8.549370  8.703057  8.728532  8.838374  8.675086
## ZYG11A  5.088814  5.284311  4.893824  3.713303  4.458707  4.310558  4.288368
## ZXDC    7.573219  8.034135  7.592597  7.572877  6.908549  6.941536  6.931567
##        GSM835366 GSM835367
## ZZZ3    7.354228  7.316619
## ZZEF1   6.698191  8.104184
## ZYX     9.527225  9.954073
## ZYG11B  8.661859  8.917929
## ZYG11A  4.587503  3.892232
## ZXDC    6.799619  7.449807
dim(probe2gene)
## [1] 41937     2
dim(probes_expr)
## [1] 25906   100
dim(genes_expr)
## [1] 13836   100

# 绘图
library("FactoMineR")
library("factoextra")
## Loading required package: ggplot2
## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa
dat.pca <- PCA(t(genes_expr) , graph = FALSE)
dat.pca
## **Results for the Principal Component Analysis (PCA)**
## The analysis was performed on 100 individuals, described by 13836 variables
## *The results are available in the following objects:
## 
##    name               description                          
## 1  "$eig"             "eigenvalues"                        
## 2  "$var"             "results for the variables"          
## 3  "$var$coord"       "coord. for the variables"           
## 4  "$var$cor"         "correlations variables - dimensions"
## 5  "$var$cos2"        "cos2 for the variables"             
## 6  "$var$contrib"     "contributions of the variables"     
## 7  "$ind"             "results for the individuals"        
## 8  "$ind$coord"       "coord. for the individuals"         
## 9  "$ind$cos2"        "cos2 for the individuals"           
## 10 "$ind$contrib"     "contributions of the individuals"   
## 11 "$call"            "summary statistics"                 
## 12 "$call$centre"     "mean of the variables"              
## 13 "$call$ecart.type" "standard error of the variables"    
## 14 "$call$row.w"      "weights for the individuals"        
## 15 "$call$col.w"      "weights for the variables"
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             legend.title = "Groups"
)

library(limma)
## 
## Attaching package: 'limma'
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
design=model.matrix(~factor(group_list))
design
##     (Intercept) factor(group_list)tumor
## 1             1                       1
## 2             1                       1
## 3             1                       1
## 4             1                       1
## 5             1                       1
## 6             1                       1
## 7             1                       1
## 8             1                       1
## 9             1                       1
## 10            1                       1
## 11            1                       1
## 12            1                       1
## 13            1                       1
## 14            1                       1
## 15            1                       1
## 16            1                       1
## 17            1                       1
## 18            1                       1
## 19            1                       1
## 20            1                       1
## 21            1                       1
## 22            1                       1
## 23            1                       1
## 24            1                       1
## 25            1                       1
## 26            1                       1
## 27            1                       1
## 28            1                       1
## 29            1                       1
## 30            1                       1
## 31            1                       1
## 32            1                       1
## 33            1                       1
## 34            1                       1
## 35            1                       1
## 36            1                       1
## 37            1                       1
## 38            1                       1
## 39            1                       1
## 40            1                       1
## 41            1                       1
## 42            1                       1
## 43            1                       1
## 44            1                       1
## 45            1                       1
## 46            1                       1
## 47            1                       1
## 48            1                       1
## 49            1                       1
## 50            1                       1
## 51            1                       1
## 52            1                       1
## 53            1                       1
## 54            1                       1
## 55            1                       1
## 56            1                       1
## 57            1                       1
## 58            1                       1
## 59            1                       1
## 60            1                       1
## 61            1                       1
## 62            1                       1
## 63            1                       1
## 64            1                       1
## 65            1                       1
## 66            1                       1
## 67            1                       1
## 68            1                       1
## 69            1                       1
## 70            1                       1
## 71            1                       1
## 72            1                       1
## 73            1                       1
## 74            1                       1
## 75            1                       1
## 76            1                       1
## 77            1                       1
## 78            1                       1
## 79            1                       1
## 80            1                       1
## 81            1                       0
## 82            1                       0
## 83            1                       0
## 84            1                       0
## 85            1                       0
## 86            1                       0
## 87            1                       0
## 88            1                       0
## 89            1                       0
## 90            1                       0
## 91            1                       0
## 92            1                       0
## 93            1                       0
## 94            1                       0
## 95            1                       0
## 96            1                       0
## 97            1                       0
## 98            1                       0
## 99            1                       0
## 100           1                       0
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$`factor(group_list)`
## [1] "contr.treatment"
fit=lmFit(genes_expr,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)
##             logFC   AveExpr         t      P.Value    adj.P.Val        B
## ZNF71   0.5732053  6.335189  8.159428 9.886311e-13 1.367870e-08 18.33355
## LPAR2  -0.4645595  6.882403 -7.347358 5.417315e-11 3.747699e-07 14.58418
## NOS2   -1.0231316  4.687544 -7.093366 1.854077e-10 7.072127e-07 13.43122
## IFT52   0.5484985  7.743018  7.015059 2.702441e-10 7.072127e-07 13.07814
## SUSD4  -1.4528754  6.202299 -6.996490 2.954455e-10 7.072127e-07 12.99459
## CHMP4B  0.5183351 10.296150  6.988710 3.066837e-10 7.072127e-07 12.95960
df=DEG
attach(df)
df$v= -log10(P.Value)
df$g=ifelse(df$P.Value>0.05,'stable',
            ifelse( df$logFC >1,'up',
                    ifelse( df$logFC < -1,'down','stable') )
)
table(df$g)
## 
##   down stable     up 
##     72  13727     37
df$name=rownames(df)
head(df)
##             logFC   AveExpr         t      P.Value    adj.P.Val        B
## ZNF71   0.5732053  6.335189  8.159428 9.886311e-13 1.367870e-08 18.33355
## LPAR2  -0.4645595  6.882403 -7.347358 5.417315e-11 3.747699e-07 14.58418
## NOS2   -1.0231316  4.687544 -7.093366 1.854077e-10 7.072127e-07 13.43122
## IFT52   0.5484985  7.743018  7.015059 2.702441e-10 7.072127e-07 13.07814
## SUSD4  -1.4528754  6.202299 -6.996490 2.954455e-10 7.072127e-07 12.99459
## CHMP4B  0.5183351 10.296150  6.988710 3.066837e-10 7.072127e-07 12.95960
##                v      g   name
## ZNF71  12.004966 stable  ZNF71
## LPAR2  10.266216 stable  LPAR2
## NOS2    9.731872   down   NOS2
## IFT52   9.568244 stable  IFT52
## SUSD4   9.529523   down  SUSD4
## CHMP4B  9.513309 stable CHMP4B
library(ggpubr)
ggpubr::ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
                  label = "name", repel = T,
                  label.select =head(rownames(df)),
                  palette = c("#00AFBB", "#E7B800", "#FC4E07") )

detach(df)
x=DEG$logFC
names(x)=rownames(DEG)
cg=c(names(head(sort(x),100)),
     names(tail(sort(x),100)))
cg
##   [1] "KRT13"     "SPRR2C"    "NTRK2"     "KRT15"     "NTS"       "DAPL1"    
##   [7] "SPRR1B"    "DSC3"      "CALML3"    "KRT5"      "UPK1B"     "GBP6"     
##  [13] "TP63"      "SERPINB3"  "DSG3"      "WIF1"      "PCP4L1"    "PITX2"    
##  [19] "SERPINB13" "DDX43"     "CYP4X1"    "PTPRZ1"    "SUSD4"     "SPRR3"    
##  [25] "IGF2BP3"   "DEFB1"     "PLA2G2A"   "SPINK5"    "KLHL13"    "MAGEA4"   
##  [31] "SERPINB2"  "HMGA2"     "RASSF9"    "SPRR1A"    "TMPRSS4"   "DDX3Y"    
##  [37] "ZNF750"    "LYPD6B"    "SCGB2A1"   "CSTA"      "HEY1"      "CDH26"    
##  [43] "P2RY1"     "PCSK2"     "SCGB1A1"   "PKP1"      "MUC5B"     "PCDHB2"   
##  [49] "EIF1AY"    "SIX2"      "MRAP2"     "CCDC190"   "CHST9"     "SERPINB4" 
##  [55] "EPHA3"     "POF1B"     "FOXE1"     "NRN1"      "ALDH1A1"   "FRAS1"    
##  [61] "GDA"       "RPS4Y1"    "SERPINI1"  "TMEM246"   "GSTT1"     "FST"      
##  [67] "RGMA"      "FAT2"      "C12orf54"  "NOS2"      "AOC1"      "RELN"     
##  [73] "AADACP1"   "CYP2C18"   "TMPRSS11D" "CLCA2"     "ALOX12"    "ARTN"     
##  [79] "CEL"       "PI3"       "FGFR2"     "CERS3"     "AADAC"     "GJB5"     
##  [85] "ADH7"      "CFAP53"    "HAS3"      "STXBP6"    "SOCS2"     "FAM110C"  
##  [91] "DQX1"      "ABCC5"     "FXYD3"     "KDM5D"     "CYP2S1"    "MYT1"     
##  [97] "JAG1"      "SLC16A1"   "SHISA2"    "CPA3"      "KCNMB4"    "EDIL3"    
## [103] "MAP3K9"    "SFRP1"     "MEGF9"     "GSE1"      "DKK1"      "COL8A1"   
## [109] "BMP2"      "PARD6B"    "CRIP1"     "ZSCAN18"   "DIRAS3"    "TDRKH"    
## [115] "GPAT2"     "ATP11A"    "PPP1R9A"   "COL10A1"   "COL11A1"   "AREG"     
## [121] "ZNF300"    "PXDN"      "COL4A4"    "NPW"       "COL12A1"   "FHOD3"    
## [127] "KBTBD11"   "DCBLD2"    "CLDN4"     "LMO3"      "MYO5C"     "SCG5"     
## [133] "KCNN4"     "GLIS3"     "NDNF"      "WWC1"      "PON3"      "HLA-DQA1" 
## [139] "VGLL3"     "DPT"       "AZGP1"     "FHL2"      "CRABP2"    "HORMAD1"  
## [145] "VWDE"      "ACSM3"     "C4BPA"     "FGA"       "ADGRA3"    "FMO5"     
## [151] "CXCL5"     "PTPN21"    "GPR39"     "CDA"       "UCHL1"     "PCOLCE2"  
## [157] "ANPEP"     "TOX"       "RASGEF1A"  "PDK4"      "MACROD2"   "CPE"      
## [163] "FLRT3"     "MET"       "STEAP4"    "DSEL"      "CHGB"      "MACC1"    
## [169] "TUBB2B"    "HLA-DRB4"  "OPN3"      "CGN"       "PRKAR2B"   "S100A7"   
## [175] "SGMS2"     "DDAH1"     "HSPA2"     "MUM1L1"    "ABCC3"     "PCSK1"    
## [181] "PNMA2"     "FGG"       "AMIGO2"    "SFTA2"     "CYP27C1"   "NT5E"     
## [187] "CST6"      "LONRF2"    "EREG"      "COBL"      "FGB"       "KRT7"     
## [193] "SCEL"      "KRT80"     "CES1"      "GRB14"     "TOX3"      "TFPI2"    
## [199] "CPS1"      "MAP2"
library(pheatmap)
n=t(scale(t(genes_expr[cg,])))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
##         GSM835268  GSM835269  GSM835270  GSM835271
## KRT13  -0.3931053 -0.4092995  0.1243226 -0.3518128
## SPRR2C -0.5212176 -0.4635166 -0.4436278 -0.6081929
## NTRK2  -0.8121804 -0.8454570 -0.8304742 -0.9134054
## KRT15   0.4560598  0.7134532  0.6573773  0.3768809
ac=data.frame(groupList=group_list)
rownames(ac)=colnames(n)  
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac)
