
##########simulate data
## require packages
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
fun3<-function(mm, Par_Score, RP_score)
{

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

####define functions(3) permutation input RP_matrix
fun_random_1000 <-function(mm, Par_Score, RP_score)
{

sort_RP_score_adj<-RP_score/max(RP_score)
sort(Par_Score,decreasing=T)->sort_Par_Score
sum(abs(sort_Par_Score)* sort_RP_score_adj)->NR_sum
P_hit<-sort_RP_score_adj*(abs(sort_Par_Score))/NR_sum
Non_Sum<-length(sort_RP_score_adj)-sum(sort_RP_score_adj)
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
cor.test(drug_confound,TF+rnorm(300,0, sd_class))->TF_P
cor.test(drug_confound, non_TF +rnorm(300,0, sd_class))->nonTF_P

#sapply(1:dim(target_matrix2)[1],function(x) pcor.test(drug_confound_last,target_matrix2[x,], confounding_dummy)[1]$estimate)->cor_res
#####partial correlation Use partial correlation to evalute the effect of the confounders and remove such effects
apply(target_matrix2,1,function(x) pcor.test(drug_confound_last,x, confounding_dummy)[1]$estimate)->cor_res
sapply(1:dim(target_matrix2)[1],function(x) pcor.test(drug_confound_last, noneTF_target_matrix2[x,], confounding_dummy)[1]$estimate)-> noncor_res
summary(glm(cor_res ~RP_scores))$coefficients[2,4]->target_P
summary(glm(noncor_res ~RP_scores))$coefficients[2,4]->non_target_P
t(sapply(1:1000,function(x) sample(RP_scores)))->randim_RP
as.numeric(cor(drug,t(target_matrix2)))->cor_res
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

###one confounders
(which(runif(300,0,1)>0.9))->confund_list1
rep(0,300)-> confound1
confound1[confund_list1]=1

drug+ confound1*2->drug_con1

drug_con1-> drug_confound_last
confound1-> confound_last

drug_confound<-drug_confound_last
confounding_dummy<-confound_last


sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm

pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.8241
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.9804



###two confounders
(which(runif(300,0,1)>0.9))->confund_list1
rep(0,300)-> confound1
confound1[confund_list1]=1
(which(runif(300,0,1)>0.9))->confund_list2
rep(0,300)-> confound2
confound2[confund_list2]=1

drug+ confound1*2+confound2*2->drug_con

drug_con-> drug_confound_last
cbind(confound1, confound2)-> confound_last

drug_confound<-drug_confound_last
confounding_dummy<-confound_last


#sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm


####zero confounders


drug_confound<-drug
confounding_dummy<-rep(0,300)
#sapply(1:100,function(x) perm_confound(drug, drug_confound= drug, confounding_dummy = confounding_dummy,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm
res_perm->res_perm_found0


pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf <- performance(pred, "tpr", "fpr")
plot(perf,col="#FFC125",lwd=2,add=F,lty=2)
auc <- performance(pred, "auc") #0.8241
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2,add=T,col="#b6d600",lwd=2)
auc2 <- performance(pred2, "auc") #0.9804


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

drug+ confound1*2+confound2*2+confound3*2->drug_con

drug_con-> drug_confound_last
cbind(confound1, confound2)-> confound_last

drug_confound<-drug_confound_last
confounding_dummy<-confound_last


sapply(1:100,function(x) perm_confound(drug, drug_confound= drug_confound_last, confounding_dummy = confound_last,r_effic=0.3,sd_class=2,chip_sd=2,target_freq=0.05))->res_perm
res_perm->res_perm_found3



res_perm_found3-> res_perm
pred<-prediction(log(c(res_perm[1,],res_perm[2,]))*(-1),rep(c(1,0),c(100,100)))
auc <- performance(pred, "auc")
apply(res_perm,2,function(x) c(max(x[c(3,5)]),max(x[c(4,6)])))->target_perm
pred2<-prediction(log(c(target_perm[1,], target_perm[2,]))*(-1),rep(c(1,0),c(100,100)))

auc2 <- performance(pred2, "auc")


rbind(c(0.9504,0.8969),c(0.9804,0.8241),c(0.9682,0.7858),c(0.9339,0.7184))->conf_res
bp<-barplot(t(conf_res),beside=T,ylim=c(0,1.3),col=c("#B8FDAF","#FEFDC5"),border="grey")



