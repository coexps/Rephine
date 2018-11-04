
##########simulate data
## require packages
library(lmtest)
library(psych)
library(ROCR)

##########simulate data
library(lmtest)
library(psych)
library(ROCR)
####define functions(1) correlated variated generation
correlatedValue = function(r,vari ){
  r2 = r**2
  ve = 1-r2
  SD = sqrt(ve)
  e  = rnorm(length(vari), mean=0, sd=SD)
  y  = r* vari + e
  return(y)
}

####define functions(2) enrichment scores
fun3<-function(mm, Par_score, RP_score)
{

RP_score_adj<-RP_score/max(RP_score)


sum(abs(Par_score)* RP_score_adj)->NR_sum
order(Par_score,decreasing=T)->order_Par_score
sort_Par_score<-Par_score[order_Par_score]
RP_score_adj[order_Par_score]-> sort_RP_score_adj
P_hit<-sort_RP_score_adj*(abs(sort_Par_score)^mm)/NR_sum
Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
P_miss<-(1-sort_RP_score_adj)/Non_Sum
P_hit-P_miss

}

fun3_GSEA_1<-function(mm, Par_Score, RP_score)
{
    sign(RP_score)->RP_score
    RP_score[RP_score<0]=0
    RP_score_adj<-RP_score/max(RP_score)
    
    
    sum(abs(Par_Score)* RP_score_adj)->NR_sum
    order(Par_Score,decreasing=T)->order_Par_Score
    sort_Par_Score<-Par_Score[order_Par_Score]
    RP_score_adj[order_Par_Score]-> sort_RP_score_adj
    P_hit<-sort_RP_score_adj*(abs(sort_Par_Score)^mm)/NR_sum
    Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
    P_miss<-(1-sort_RP_score_adj)/Non_Sum
    P_hit-P_miss

    
}



fun3_GSEA_2<-function(mm, Par_score, RP_score)
{
    
    as.numeric(RP_score>sort(RP_score,decreasing=T)[1001])->RP_score
    RP_score_adj<-RP_score/max(RP_score)
    
    sum(abs(Par_score)* RP_score_adj)->NR_sum
    order(Par_score,decreasing=T)->order_Par_score
    sort_Par_score<-Par_score[order_Par_score]
    RP_score_adj[order_Par_score]-> sort_RP_score_adj
    P_hit<-sort_RP_score_adj*(abs(sort_Par_score)^mm)/NR_sum
    Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
    P_miss<-(1-sort_RP_score_adj)/Non_Sum
    P_hit-P_miss

    
}


####define functions(3) permutation input RP_matrix

fun_random_1000 <-function(mm, Par_score, RP_score)
{

RP_score_adj<-RP_score/max(RP_score)
sum(abs(Par_score)* RP_score_adj)->NR_sum
order(Par_score,decreasing=T)->order_Par_score
sort_Par_score<-Par_score[order_Par_score]
RP_score_adj[order_Par_score]-> sort_RP_score_adj
P_hit<-sort_RP_score_adj*(abs(sort_Par_score)^mm)/NR_sum
Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
P_miss<-(1-sort_RP_score_adj)/Non_Sum
cumsum(P_hit-P_miss)->ress
return(c(max(ress),min(ress)))
}



fun_random_1000_GSEA_1 <-function(mm, Par_Score, RP_score)
{
    sign(RP_score)->RP_score
    RP_score[RP_score<0]=0

    sort_RP_score_adj<-RP_score/max(RP_score)
    
    
    sum(abs(Par_score)* RP_score_adj)->NR_sum
    order(Par_score,decreasing=T)->order_Par_Score
    sort_Par_Score<-Par_Score[order_Par_score]
    RP_score_adj[order_Par_Score]-> sort_RP_score_adj
    P_hit<-sort_RP_score_adj*(abs(sort_Par_Score)^mm)/NR_sum
    Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
    P_miss<-(1-sort_RP_score_adj)/Non_Sum
    cumsum(P_hit-P_miss)->ress
    return(c(max(ress),min(ress)))
}

fun_random_1000_GSEA_2 <-function(mm, Par_score, RP_score)
{



    as.numeric(RP_score>sort(RP_score,decreasing=T)[1001])->RP_core
    RP_score_adj<-RP_score/max(RP_score)

    
    sum(abs(Par_score)* RP_score_adj)->NR_sum
    order(Par_score,decreasing=T)->order_Par_score
    sort_Par_score<-Par_score[order_Par_score]
    RP_score_adj[order_Par_score]-> sort_RP_score_adj
    P_hit<-sort_RP_score_adj*(abs(sort_Par_score)^mm)/NR_sum
    Non_Sum<-length(RP_score_adj)-sum(sort_RP_score_adj)
    P_miss<-(1-sort_RP_score_adj)/Non_Sum
    cumsum(P_hit-P_miss)->ress
    return(c(max(ress),min(ress)))
}


######start steps simulate Rephine performance with different levels of noise
#generation data
drug = rnorm(300,mean=5)  ###generate drug value

#difine a function in for one simulation run 1) simulate expression data for TR, target and non-target expression, random noise and RP scores.2) calculte the Rephine significance univariate P and permutation P as well as the correlation P.
perm<-function(drug , r_effic, sd_class, chip_sd, target_freq) {
TF = correlatedValue(vari=drug, r= r_effic) ###generate TF exp
non_TF = correlatedValue(vari=drug, r=0)
c(runif(target_freq*1e4,0.000001,10),rep(0,1e4-target_freq*1e4))+rnorm(10000,0, chip_sd)->RP_scores ## ChIP_rp generation
cor_noise<-as.vector(c(fisherz2r(rnorm(target_freq*1e4,fisherz(r_effic),1)), fisherz2r (rnorm(1e4-target_freq*1e4,fisherz(0),1))))   ###generate target_correlation with noise
nonTF_cor_noise<-c(fisherz2r(rnorm(10000,fisherz(0),1)))  ###generate nonTFtarget_correlation with noise
noise,function(x) correlatedValue(x, drug)))-> target_matrix #generate TF target expression
t(sapply(nonTF_cor_noise,function(x) correlatedValue(x, drug)))-> nonTF_all_matrix #generate nonTF target expression
#adding noise
target_matrix2<-t(apply(target_matrix[,],1,function(x) as.numeric(x+rnorm(300,0, sd_class)))) ##add noise
noneTF_target_matrix2<-t(apply(nonTF_all_matrix,1,function(x) as.numeric(x+rnorm(300,0, sd_class))))##add noise


#calculation of Rephine univariate P-value
cor.test(drug,TF+rnorm(300,0, sd_class))->TF_P
cor.test(drug, non_TF +rnorm(300,0, sd_class))->nonTF_P
as.numeric(cor(drug,t(target_matrix2)))->cor_res #####partial correlation
as.numeric(cor(drug,t(noneTF_target_matrix2)))->noncor_res #####partial correlation

summary(glm(cor_res ~RP_scores))$coefficients[2,4]->target_P #univariate P of the Rephine of positve TR
summary(glm(noncor_res ~RP_scores))$coefficients[2,4]->non_target_P  #univariate P of the Rephine of non-related TR

#calculation of Rephine permutation P-value
t(sapply(1:1000,function(x) sample(RP_scores)))->randim_RP # permutate the RP scores
as.numeric(cor(drug,t(target_matrix2)))->cor_res
apply(randim_RP,1,function(x) fun_random_1000(1, cor_res,x))->TF_maxmin #permutate for 1000 times
apply(randim_RP,1,function(x) fun_random_1000(1, noncor_res,x))->nonTF_maxmin #permutate for 1000 times


cumsum(fun3(1, cor_res, RP_scores))->trueT #calculate the enrichment scores
max(trueT)->max_ture
min(trueT)->min_ture
cumsum(fun3(1, noncor_res, RP_scores))->falseF
max(falseF)->max_false
min(falseF)->min_false

min(1-sum(max_ture>TF_maxmin[1,])/1000,1-sum(min_ture<TF_maxmin[2,])/1000)->target_P_en #calculate the permutation significance.
min(1-sum(max_false> nonTF_maxmin[1,])/1000,1-sum(min_false <nonTF_maxmin[2,])/1000)->nontarget_P_en

return(c(TF_P$p.val, nonTF_P$p.val, target_P, non_target_P, target_P_en, nontarget_P_en))
}




#########results1 compare with noise

#sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm

res_perm_bak1-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#AB82FF",lwd=2, lty=2)
auc <- performance(pred, "auc") #0.89 #res_perm_bak1

apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#749EC8",lwd=2)
auc2 <- performance(pred2, "auc") #0.962

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=2,target_freq=0.05))->res_perm
#res_perm_bak2-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFDAB9",lwd=2,add=T,lty=2)
auc <- performance(pred, "auc") #0.6144  #res_perm_bak2

apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#00B2EE",lwd=2)
auc2 <- performance(pred2, "auc") #0.9415

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=1,chip_sd=2,target_freq=0.05))->res_perm
#res_perm_bak3-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=T,lty=2)
auc <- performance(pred, "auc") #0.9975

#res_perm-> res_perm_bak3
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#FFC1C1",lwd=2)
auc2 <- performance(pred2, "auc") #0.9497

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=10,chip_sd=2,target_freq=0.05))->res_perm
#res_perm->res_perm_bak4
#res_perm-> res_perm_bak4
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=T,lty=2)
auc <- performance(pred, "auc") #0.5733

apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#FF8247",lwd=2)
auc2 <- performance(pred2, "auc") #0.8904

#sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=20,chip_sd=2,target_freq=0.05))->res_perm

#res_perm-> res_perm_bak5
res_perm_bak5-> res_perm
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#90EE90",lwd=2)
auc2 <- performance(pred2, "auc")  #0.71275





sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=15,chip_sd=2,target_freq=0.05))->res_perm
#res_perm->res_perm_bak6
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=T,lty=2)
auc <- performance(pred, "auc") #0.5111

res_perm-> res_perm_bak3
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#FFC1C1",lwd=2)
auc2 <- performance(pred2, "auc") #0.8254

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=30,chip_sd=2,target_freq=0.05))->res_perm
#res_perm->res_perm_bak7
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=T,lty=2)
auc <- performance(pred, "auc") #0.5175

#res_perm-> res_perm_bak3
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#FFC1C1",lwd=2)
auc2 <- performance(pred2, "auc") #0.6585

#draw the roc curve
TF_cor_AUC<-c(0.8922, 0.6144, 0.9975, 0.5733,0.5038,0.5111, 0.5175)
target_cor_AUC<-c(0.962, 0.9415, 0.9497, 0.8904, 0.71275, 0.8254, 0.6585)
core<-c(0.3,0.3,0.3,0.3,0.3,0.3,0.3)
sd_noise_plot<-c(2,5,1,10,20,15,30)
cbind(TF_cor_AUC, target_cor_AUC, core, sd_noise_plot)->compare_res
compare_res[order(compare_res[,4]),]-> compare_res2
bp<-barplot(rbind(compare_res2[,1],compare_res2[,2]),beside=T,ylim=c(0,1.3),col=c("#B8FDAF","#FEFDC5"),border="grey")
axis(4,at=seq(0,1.2,0.3),label=seq(0,1.2,0.3)*30)
lines(bp[1,]+0.5,compare_res2[,4]/30,col="#CEB5FB",lwd=4,lty=3)






#######compare 2 performance diff cor


sapply(1:100,function(x) perm(drug,r_effic=0.4,sd_class=5,chip_sd=2,target_freq=0.05))->res_perm
#res_perm->res_perm_roc_1
res_perm_roc_1 -> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#A9A9A9",add=F,lty=3,lwd=3) #0.8352
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#FF0000",lwd=2)
auc2 <- performance(pred2, "auc") #0.9976

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=2,target_freq=0.05))->res_perm
#res_perm-> res_perm_roc_2
res_perm_roc_2-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#FF9301",lwd=2)
auc2 <- performance(pred2, "auc") #0.93635

sapply(1:100,function(x) perm(drug,r_effic=0.2,sd_class=5,chip_sd=2,target_freq=0.05))->res_perm
#res_perm-> res_perm_roc_3
res_perm_roc_3-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=4,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#6560FF",lwd=2)
auc2 <- performance(pred2, "auc") #0.8114



sapply(1:100,function(x) perm(drug,r_effic=0.25,sd_class=5,chip_sd=2,target_freq=0.05))->res_perm
#res_perm-> res_perm_roc_4
res_perm_roc_4-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#0C9508",lwd=2)
auc2 <- performance(pred2, "auc") #0.889


########## add chip noise

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=1,target_freq=0.05))->res_perm
#res_perm-> res_perm_chip_1
res_perm_chip_1-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=F,col="#FF0000",lwd=2)
auc2 <- performance(pred2, "auc") #0.9896



sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=5,target_freq=0.05))->res_perm
#res_perm-> res_perm_chip_2
res_perm_chip_2-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#6560FF",lwd=2)
auc2 <- performance(pred2, "auc") #0.6642

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=3,target_freq=0.05))->res_perm
#res_perm-> res_perm_chip_3
res_perm_chip_3-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.8372

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=4,target_freq=0.05))->res_perm
#res_perm-> res_perm_chip_4
res_perm_chip_4-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#0C9508",lwd=2)
auc2 <- performance(pred2, "auc") #0.75145

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=2,target_freq=0.05))->res_perm
#res_perm-> res_perm_chip_5
res_perm_chip_5-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#FF9301",lwd=2)
auc2 <- performance(pred2, "auc") #0.9378


sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=5,target_freq=0.1))->res_perm
#res_perm-> res_perm_chip_6
res_perm_chip_6-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#A9A9A9",lwd=3,lty=3)
auc2 <- performance(pred2, "auc") #0.8694




######## chip target diff


sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=2,target_freq=0.1))->res_perm
#res_perm-> res_perm_target_1
res_perm_target_1-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=F,col="#FF0000",lwd=2)
auc2 <- performance(pred2, "auc") #0.9962


#sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=2,target_freq=0.03))->res_perm
#res_perm-> res_perm_target_2
res_perm_target_2-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.81485


sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=2,target_freq=0.05))->res_perm
#res_perm-> res_perm_target_3
res_perm_target_3-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#FF9301",lwd=2)
auc2 <- performance(pred2, "auc") #0.9393

sapply(1:100,function(x) perm(drug,r_effic=0.3,sd_class=5,chip_sd=2,target_freq=0.02))->res_perm
#res_perm-> res_perm_target_4
res_perm_target_4-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
#plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5175
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#6560FF",lwd=2)
auc2 <- performance(pred2, "auc") #0.68285

######simulation code for performance comparison between RePhine and traditional correlation method with different confounding existing


#generation data
library(ppcor)
library(lmtest)
library(psych)
library(ROCR)

drug = rnorm(300,mean=5)  ###generate drug value


perm_confound<-function(drug,drug_confound, confounding_dummy, r_effic, sd_class, chip_sd, target_freq) {
TF = correlatedValue(vari=drug, r= r_effic) ###generate TF exp
non_TF = correlatedValue(vari=drug, r=0)
c(runif(target_freq*1e4,0.000001,10),rep(0,1e4-target_freq*1e4))+rnorm(10000,0, chip_sd)->RP_scores ## ChIP_rp generation
cor_noise<-as.vector(c(fisherz2r(rnorm(target_freq*1e4,fisherz(r_effic),1)), fisherz2r (rnorm(1e4-target_freq*1e4,fisherz(0),1))))   ###generate target_correlation with noise
nonTF_cor_noise<-c(fisherz2r(rnorm(10000,fisherz(0),1)))  ###generate nonTFtarget_correlation with noise
t(sapply(cor_noise,function(x) correlatedValue(x, drug)))-> target_matrix #generate TF target expression
t(sapply(nonTF_cor_noise,function(x) correlatedValue(x, drug)))-> nonTF_all_matrix #generate nonTF target expression
#adding noise
target_matrix2<-t(apply(target_matrix[,],1,function(x) as.numeric(x+rnorm(300,0, sd_class)))) ##add noise

noneTF_target_matrix2<-t(apply(nonTF_all_matrix,1,function(x) as.numeric(x+rnorm(300,0, sd_class))))##add noise

#calculation
TF+rnorm(300,0, sd_class)->TF_exp
cor.test(drug_confound,TF_exp)->TF_P
non_TF +rnorm(300,0, sd_class)->non_TF_exp
cor.test(drug_confound, non_TF_exp)->nonTF_P

#sapply(1:dim(target_matrix2)[1],function(x) pcor.test(drug_confound_last,target_matrix2[x,], confounding_dummy)[1]$estimate)->cor_res 
#####partial correlation
apply(target_matrix2,1,function(x) pcor.test(drug_confound,x, confounding_dummy)[1]$estimate)->cor_res
#apply(target_matrix2,1,function(x) cor.test(drug_confound,x)$estimate)->normal_cor_res
as.numeric(cor(drug_confound,t(target_matrix2)))->normal_cor_res

sapply(1:dim(target_matrix2)[1],function(x) pcor.test(drug_confound, noneTF_target_matrix2[x,], confounding_dummy)[1]$estimate)-> noncor_res

#sapply(1:dim(target_matrix2)[1],function(x) cor.test(drug_confound, noneTF_target_matrix2[x,])$estimate)-> normal_noncor_res
as.numeric(cor(drug_confound,t(noneTF_target_matrix2[,])))->normal_noncor_res


summary(glm(cor_res ~RP_scores))$coefficients[2,4]->target_P
summary(glm(noncor_res ~RP_scores))$coefficients[2,4]->non_target_P

t(sapply(1:1000,function(x) sample(RP_scores)))->randim_RP


#as.numeric(cor(drug,t(target_matrix2)))->cor_res

apply(randim_RP,1,function(x) fun_random_1000(1, cor_res,x))->TF_maxmin
apply(randim_RP,1,function(x) fun_random_1000(1, noncor_res,x))->nonTF_maxmin

cumsum(fun3(1, cor_res, RP_scores))->trueT
max(trueT)->max_ture
min(trueT)->min_ture
cumsum(fun3(1, noncor_res, RP_scores))->falseF
max(falseF)->max_false
min(falseF)->min_false

min(1-sum(max_ture>TF_maxmin[1,])/1000,1-sum(min_ture<TF_maxmin[2,])/1000)->target_P_en
min(1-sum(max_false> nonTF_maxmin[1,])/1000,1-sum(min_false <nonTF_maxmin[2,])/1000)->nontarget_P_en

#apply(randim_RP,1,function(x) fun_random_1000_GSEA_1(1, cor_res,x))->TF_maxmin_GSEA_1
#apply(randim_RP,1,function(x) fun_random_1000_GSEA_1(1, noncor_res,x))->nonTF_maxmin_GSEA_1

#cumsum(fun3(1, normal_cor_res, RP_scores))->trueT_GSEA_1
#max(trueT_GSEA_1)->max_ture_GSEA_1
#min(trueT_GSEA_1)->min_ture_GSEA_1
#cumsum(fun3(1, normal_noncor_res, RP_scores))->falseF_GSEA_1
#max(falseF_GSEA_1)->max_false_GSEA_1
#min(falseF_GSEA_1)->min_false_GSEA_1

#min(1-sum(max_ture_GSEA_1>TF_maxmin_GSEA_1[1,])/1000,1-sum(min_ture_GSEA_1<TF_maxmin_GSEA_1[2,])/1000)->target_P_en_GSEA_1
#min(1-sum(max_false_GSEA_1> nonTF_maxmin_GSEA_1[1,])/1000,1-sum(min_false_GSEA_1 <nonTF_maxmin_GSEA_1[2,])/1000)->nontarget_P_en_GSEA_1


apply(randim_RP,1,function(x) fun_random_1000_GSEA_2(1, normal_cor_res,x))->TF_maxmin_GSEA_2
apply(randim_RP,1,function(x) fun_random_1000_GSEA_2(1, normal_noncor_res,x))->nonTF_maxmin_GSEA_2

cumsum(fun3_GSEA_2(1, normal_cor_res, RP_scores))->trueT_GSEA_2
max(trueT_GSEA_2)->max_ture_GSEA_2
min(trueT_GSEA_2)->min_ture_GSEA_2
cumsum(fun3_GSEA_2(1, normal_noncor_res, RP_scores))->falseF_GSEA_2
max(falseF_GSEA_2)->max_false_GSEA_2
min(falseF_GSEA_2)->min_false_GSEA_2
min(1-sum(max_ture_GSEA_2>TF_maxmin_GSEA_2[1,])/1000,1-sum(min_ture_GSEA_2<TF_maxmin_GSEA_2[2,])/1000)->target_P_en_GSEA_2
min(1-sum(max_false_GSEA_2> nonTF_maxmin_GSEA_2[1,])/1000,1-sum(min_false_GSEA_2 <nonTF_maxmin_GSEA_2[2,])/1000)->nontarget_P_en_GSEA_2
##logistic_ID
summary(glm(as.numeric(drug_confound>median(drug_confound))~TF_exp,family="binomial"))$coefficients[2,4]->logistic_TF_P
summary(glm(as.numeric(drug_confound>median(drug_confound))~non_TF_exp,family="binomial"))$coefficients[2,4]->logistic_non_TF_P
return(c(TF_P$p.val, nonTF_P$p.val, target_P, non_target_P, target_P_en, nontarget_P_en,target_P_en_GSEA_2,nontarget_P_en_GSEA_2,logistic_TF_P,logistic_non_TF_P))
}

 ###one confounders
(which(runif(300,0,1)>0.9))->confund_list1
rep(0,300)-> confound1
confound1[confund_list1]=1

drug+ confound1*5->drug_con1

drug_con1-> drug_confound_last
confound1-> confound_last

#drug_confound<-drug_confound_last
#confounding_dummy<-confound_last


sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm

pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.7056
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.95905
#GSEA
pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.94605
##logsitc
pred4<-prediction(log(c(res_perm[9,],res_perm[10,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf4 <- performance(pred4, "tpr", "fpr")
auc4 <- performance(pred4, "auc") #0.747




 ###two confounders
(which(runif(300,0,1)>0.9))->confund_list1
rep(0,300)-> confound1
confound1[confund_list1]=1
(which(runif(300,0,1)>0.9))->confund_list2
rep(0,300)-> confound2
confound2[confund_list2]=1

drug+ confound1*5+confound2*5->drug_con

drug_con-> drug_confound_last
cbind(confound1, confound2)-> confound_last

drug_confound<-drug_confound_last
confounding_dummy<-confound_last


#sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm

pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.8588
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.9679
#GSEA
pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.95365
##logsitc
pred4<-prediction(log(c(res_perm[9,],res_perm[10,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf4 <- performance(pred4, "tpr", "fpr")
auc4 <- performance(pred4, "auc") #0.7194




####three confounders
(which(runif(300,0,1)>0.9))->confund_list1
rep(0,300)-> confound1
confound1[confund_list1]=1
(which(runif(300,0,1)>0.9))->confund_list2
rep(0,300)-> confound2
confound2[confund_list2]=1
rep(0,300)-> confound3
(which(runif(300,0,1)>0.9))->confund_list3
confound3[confund_list3]=1

drug+ confound1*5+confound2*5+confound3*5->drug_con

drug_con-> drug_confound_last
cbind(confound1, confound2,confound3)-> confound_last

drug_confound<-drug_confound_last
confounding_dummy<-confound_last


sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm
res_perm->res_perm_found3



#res_perm_found3-> res_perm

pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.6401
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.9514
#GSEA
pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.89
##logsitc
pred4<-prediction(log(c(res_perm[9,],res_perm[10,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf4 <- performance(pred4, "tpr", "fpr")
auc4 <- performance(pred4, "auc") #0.7375

#auc2 <- performance(pred2, "auc")

sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=3,target_freq=0.05))->res_perm


#confounder 5
(which(runif(300,0,1)>0.9))->confund_list1
rep(0,300)-> confound1
confound1[confund_list1]=1
(which(runif(300,0,1)>0.9))->confund_list2
rep(0,300)-> confound2
confound2[confund_list2]=1
rep(0,300)-> confound3
(which(runif(300,0,1)>0.9))->confund_list3
confound3[confund_list3]=1
rep(0,300)-> confound4
(which(runif(300,0,1)>0.9))->confund_list4
confound4[confund_list4]=1
rep(0,300)-> confound5
(which(runif(300,0,1)>0.9))->confund_list5
confound5[confund_list5]=1

drug+ confound1*5+confound2*5+confound3*5+confound4*5+confound5*5->drug_con

drug_con-> drug_confound_last
cbind(confound1, confound2,confound3,confound4,confound5)-> confound_last

drug_confound<-drug_confound_last
confounding_dummy<-confound_last

sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm


pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.6165 0.5191
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.9686 0.9819
#GSEA
pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.89835 0.84505
##logsitc
pred4<-prediction(log(c(res_perm[9,],res_perm[10,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf4 <- performance(pred4, "tpr", "fpr")
auc4 <- performance(pred4, "auc") #0.5957 0.6562

 # 7 confounders
(which(runif(300,0,1)>0.9))->confund_list1
rep(0,300)-> confound1
confound1[confund_list1]=1
(which(runif(300,0,1)>0.9))->confund_list2
rep(0,300)-> confound2
confound2[confund_list2]=1
rep(0,300)-> confound3
(which(runif(300,0,1)>0.9))->confund_list3
confound3[confund_list3]=1
rep(0,300)-> confound4
(which(runif(300,0,1)>0.9))->confund_list4
confound4[confund_list4]=1
rep(0,300)-> confound5
(which(runif(300,0,1)>0.9))->confund_list5
confound5[confund_list5]=1
rep(0,300)-> confound6
(which(runif(300,0,1)>0.9))->confund_list6
confound6[confund_list6]=1
rep(0,300)-> confound7
(which(runif(300,0,1)>0.9))->confund_list7
confound7[confund_list7]=1

drug+ rowSums(cbind(confound1+confound2+confound3+confound4+confound5+confound6+confound7)*5)->drug_con

drug_con-> drug_confound_last
cbind(confound1, confound2,confound3,confound4,confound5,confound6,confound7)-> confound_last

drug_confound<-drug_confound_last
confounding_dummy<-confound_last

sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm

pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.5257 0.5508
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.9659 0.97765
#GSEA
pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.80725 0.844
##logsitc
pred4<-prediction(log(c(res_perm[9,],res_perm[10,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf4 <- performance(pred4, "tpr", "fpr")
auc4 <- performance(pred4, "auc") #0.5247 0.5074



drug_confound<-drug
confounding_dummy<-rep(0,300)
#sapply(1:100,function(x) perm_confound(drug, drug_confound= drug, confounding_dummy = confounding_dummy,r_effic=0.3,sd_class=2,chip_sd=4,target_freq=0.05))->res_perm
res_perm->res_perm_found0



####zero confounders compare the performence of GSEA and RePhine with noise of ChIP-seq


drug_confound<-drug
confounding_dummy<-rep(0,300)

sapply(1:100,function(x) perm_confound(drug, drug_confound= drug, confounding_dummy = confounding_dummy,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm


apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.9761
#GSEA
pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") # 0.95295





sapply(1:100,function(x) perm_confound(drug, drug_confound= drug, confounding_dummy = confounding_dummy,r_effic=0.3,sd_class=2,chip_sd=3,target_freq=0.05))->res_perm

apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
auc2 <- performance(pred2, "auc") #0.9119

pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.84395


sapply(1:100,function(x) perm_confound(drug, drug_confound= drug, confounding_dummy = confounding_dummy,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm

apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
auc2 <- performance(pred2, "auc") # 0.9761 0.9656  twice with similar AUC

pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") # 0.95295 0.95195 twice with similar AUC


 sapply(1:100,function(x) perm_confound(drug, drug_confound= drug, confounding_dummy = confounding_dummy,r_effic=0.3,sd_class=2,chip_sd=4,target_freq=0.05))->res_perm

apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
auc2 <- performance(pred2, "auc") #0.8118

pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.72665



sapply(1:100,function(x) perm_confound(drug, drug_confound= drug, confounding_dummy = confounding_dummy,r_effic=0.3,sd_class=2,chip_sd=5,target_freq=0.05))->res_perm

apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
#plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.7594
#GSEA
pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.68555


sapply(1:100,function(x) perm_confound(drug, drug_confound= drug, confounding_dummy = confounding_dummy,r_effic=0.3,sd_class=2,chip_sd=6,target_freq=0.05))->res_perm
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
auc2 <- performance(pred2, "auc") #0.6494

pred3<-prediction(log(c(res_perm[7,],res_perm[8,])+0.0001)*(-1),rep(c(1,0),c(100,100)))
perf3 <- performance(pred3, "tpr", "fpr")
auc3 <- performance(pred3, "auc") #0.55





#confounders #1,3,5,7 plots based on the AUC above.
rbind(c(0.9047,0.9761,0.95295,0.8308),c(0.7056,0.95905,0.94605,0.747),c(0.6401,0.9514,0.89,0.7375),c(0.5191,0.9819,0.84505,0.6562),c(0.5257,0.9659,0.80725,0.5247))->conf_res
bp<-barplot(t(conf_res[,c(1,4,3,2)]),beside=T,ylim=c(0,1.3),col=c("#AED6F1","#B8FDAF","#FEFDC5","#FAD7A0"),border="grey")

#2,4,5,6,7
rbind(c(0.95195,0.9656),c(0.72665,0.8118),c(0.68555,0.7594),c(.55,0.6494),c())->conf_res
bp<-barplot(t(conf_res[]),beside=T,ylim=c(0,1.3),col=c("#FEFDC5","#FAD7A0"),border="grey")




perm_sd_no_GSEA<-function(drug,drug_confound, r_effic, sd_class, chip_sd, target_freq) {
    TF = correlatedValue(vari=drug, r= r_effic) ###generate TF exp
    non_TF = correlatedValue(vari=drug, r=0)
    c(runif(target_freq*1e4,0.000001,10),rep(0,1e4-target_freq*1e4))+rnorm(10000,0, chip_sd)->RP_scores ## ChIP_rp generation
    cor_noise<-as.vector(c(fisherz2r(rnorm(target_freq*1e4,fisherz(r_effic),1)), fisherz2r (rnorm(1e4-target_freq*1e4,fisherz(0),1))))   ###generate target_correlation with noise
    nonTF_cor_noise<-c(fisherz2r(rnorm(10000,fisherz(0),1)))  ###generate nonTFtarget_correlation with noise
    t(sapply(cor_noise,function(x) correlatedValue(x, drug)))-> target_matrix #generate TF target expression
    t(sapply(nonTF_cor_noise,function(x) correlatedValue(x, drug)))-> nonTF_all_matrix #generate nonTF target expression
    #adding noise
    target_matrix2<-t(apply(target_matrix[,],1,function(x) as.numeric(x+rnorm(300,0, sd_class)))) ##add noise
    
    noneTF_target_matrix2<-t(apply(nonTF_all_matrix,1,function(x) as.numeric(x+rnorm(300,0, sd_class))))##add noise
    
    #calculation
    TF+rnorm(300,0, sd_class)->TF_exp
    cor.test(drug_confound,TF_exp)->TF_P
    non_TF +rnorm(300,0, sd_class)->non_TF_exp
    cor.test(drug_confound, non_TF_exp)->nonTF_P
    
    #####partial correlation without CF
    as.numeric(cor(drug_confound,t(target_matrix2)))->cor_res
    
    
    as.numeric(cor(drug_confound,t(noneTF_target_matrix2[,])))->noncor_res
    
    
    summary(glm(cor_res ~RP_scores))$coefficients[2,4]->target_P
    summary(glm(noncor_res ~RP_scores))$coefficients[2,4]->non_target_P
    
    t(sapply(1:1000,function(x) sample(RP_scores)))->randim_RP
    
    
    
    apply(randim_RP,1,function(x) fun_random_1000(1, cor_res,x))->TF_maxmin
    apply(randim_RP,1,function(x) fun_random_1000(1, noncor_res,x))->nonTF_maxmin
    
    cumsum(fun3(1, cor_res, RP_scores))->trueT
    max(trueT)->max_ture
    min(trueT)->min_ture
    cumsum(fun3(1, noncor_res, RP_scores))->falseF
    max(falseF)->max_false
    min(falseF)->min_false
    
    min(1-sum(max_ture>TF_maxmin[1,])/1000,1-sum(min_ture<TF_maxmin[2,])/1000)->target_P_en
    min(1-sum(max_false> nonTF_maxmin[1,])/1000,1-sum(min_false <nonTF_maxmin[2,])/1000)->nontarget_P_en
    
 
    return(c(TF_P$p.val, nonTF_P$p.val, target_P, non_target_P, target_P_en, nontarget_P_en))
}
          


      
sapply(1:100,function(x )perm_sd_no_GSEA(drug=drug,drug_confound= drug,r_effic=0.3, sd_class=5, chip_sd=1, target_freq=0.05))->rr1
sapply(1:100,function(x )perm_sd_no_GSEA(drug=drug,drug_confound= drug,r_effic=0.3, sd_class=5, chip_sd=2, target_freq=0.05))->rr2
sapply(1:100,function(x )perm_sd_no_GSEA(drug=drug,drug_confound= drug,r_effic=0.3, sd_class=5, chip_sd=4, target_freq=0.05))->rr3
sapply(1:100,function(x )perm_sd_no_GSEA(drug=drug,drug_confound= drug,r_effic=0.3, sd_class=5, chip_sd=6, target_freq=0.05))->rr4
sapply(1:100,function(x )perm_sd_no_GSEA(drug=drug,drug_confound= drug,r_effic=0.3, sd_class=5, chip_sd=6, target_freq=0.1))->rr5

#FF0000 #6560FF #b6d600 #0C9508 #FF9301 #A9A9A9

apply(rr1,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,col="#FF0000",add=F,lty=1,lwd=2)
auc2 <- performance(pred2, "auc") #0.9988

apply(rr2,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
auc2 <- performance(pred2, "auc") #0.9398
plot(perf2,col="#FF9301",add=T,lty=1,lwd=2)


apply(rr3,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
auc2 <- performance(pred2, "auc") #0.7469
plot(perf2,col="#b6d600",add=T,lty=1,lwd=2)

apply(rr4,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
auc2 <- performance(pred2, "auc") #0.6305
plot(perf2,col="#6560FF",add=T,lty=1,lwd=2)


apply(rr5,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
auc2 <- performance(pred2, "auc") #0.8527
plot(perf2,col="#A9A9A9",add=T,lty=2,lwd=2)


