#Geographical Random Forest 


library(lctools)
library(SpatialML)
library(caret)
library(GWmodel)
library(sp)
library(randomForest)
library(readr)
library(rgdal)
library(ranger)



#import csv for GRF and OLS models
data <- read_csv("FinalAstma.csv")
#import shapefile for GWR models
spatialdataset <- readOGR("Asthma.shp")
data=as.data.frame(data)


##caluclate bw for GWR
#108 is found as optimal
gwr_bw=bw.gwr(asthma~insurance+cancer+UrbanInten+NoInternet+UNEMP+singlepare+LimEnglish+Housing_T+NDVI+VacantH, spatialdataset, approach="AICc", kernel="bisquare",
              adaptive=TRUE,p=2, theta=0,longlat=F)




#GRF optimize mtry parameter
### cross validate random forest using grid search
control <- trainControl(method="repeatedcv",repeats=5, number=10, search="grid")
set.seed(123)
tunegrid <- expand.grid(.mtry=c(1:10),.splitrule = "variance",.min.node.size = c(5))
rf_gridsearch <- train(asthma~insurance+cancer+UrbanIntense+NoInternet+UNEMP+singleparent+LimEnglish+Housing_T+NDVI+VacantH, data= data, method="ranger",num.trees = 500, importance = "permutation",tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)
#7 is best value for mtry parameter

###Detect bw for GRF
eval_BW_GRF<- data.frame(Local=double(),
                         stringsAsFactors=FALSE
)

set.seed(12345)
for(bw in seq(from=60, to=120, by=2)){
  #set.seed(bw)
  grf16.a <- grf(asthma~insurance+cancer+UrbanIntense+NoInternet+UNEMP+singleparent+LimEnglish+Housing_T+NDVI+VacantH, dframe=data,bw=bw, kernel="adaptive",coords=data[,13:14],weighted=TRUE,ntree=500,mtry=7, nthreads = 6)
  eval_BW_GRF[bw,1]=grf16.a$LocalModelSummary$l.r.OOB
}
#bw with highest r-squared = 108   



eval_grf <- data.frame(R_squared=double(),
                       MAE=double(),
                       RMSE=double(),
                       stringsAsFactors=FALSE
)


eval_ols = data.frame(R_squared=double(),
                      MAE=double(),
                      RMSE=double(),
                      stringsAsFactors=FALSE)
eval_gwr = data.frame(R_squared=double(),
                      MAE=double(),
                      RMSE=double(),
                      stringsAsFactors=FALSE)



## random testing and splitting 10 repeats ##
for (i in 1:10){
  set.seed(i)
  smp_size <- floor(0.80 * nrow(data))
  train_ind <- sample(seq_len(nrow(data)), size = smp_size)
  train <- data[train_ind, ]
  test <- data[-train_ind, ]
  gwr_train=subset(spatialdataset, as.numeric(spatialdataset$GEOID) %in% train$FIPS )
  gwr_test=subset(spatialdataset, as.numeric(spatialdataset$GEOID) %in% test$FIPS )
  gwr.basic=gwr.predict(asthma~insurance+cancer+UrbanInten+NoInternet+UNEMP+singlepare+LimEnglish+Housing_T+NDVI+VacantH, bw=108,data=gwr_train, predictdata=gwr_test, kernel="bisquare",adaptive=TRUE, p=2,
                        theta=0, longlat=F)
  grf <- grf(asthma~insurance+cancer+UrbanIntense+NoInternet+UNEMP+singleparent+LimEnglish+Housing_T+NDVI+VacantH, importance="permutation",weighted=TRUE,dframe=train, bw=108,nthreads=6,ntree=500,kernel="adaptive",mtry=7, coords=train[,13:14])
  glm=lm(asthma~insurance+cancer+UrbanIntense+NoInternet+UNEMP+singleparent+LimEnglish+Housing_T+NDVI+VacantH, data=train)
  lm_pred= as.data.frame(predict(glm,newdata=test))
  grf_pred=predict.grf(grf, test, x.var.name="X", y.var.name="Y", local.w=1, global.w=0)
  metrics_grf=postResample(pred = grf_pred, test$asthma)
  metrics_lm=postResample(pred = lm_pred[,1],  test$asthma)
  metrics_gwr=postResample(pred=gwr.basic[["SDF"]]@data[["prediction"]],gwr_test$asthma)
  eval_grf[i,1]=metrics_grf[2]
  eval_grf[i,2]=metrics_grf[3]
  eval_grf[i,3]=metrics_grf[1]
  eval_ols[i,1]=metrics_lm[2]
  eval_ols[i,2]=metrics_lm[3]
  eval_ols[i,3]=metrics_lm[1]
  eval_gwr[i,1]=metrics_gwr[2]
  eval_gwr[i,2]=metrics_gwr[3]
  eval_gwr[i,3]=metrics_gwr[1]
}

##extract final results
mean(eval_ols$R_squared)
mean( eval_grf$R_squared)
mean( eval_grf_lm_050$R_squared)
mean( eval_gwr$R_squared)
mean(eval_ols$RMSE)
mean( eval_grf$RMSE)
mean( eval_grf_lm_050$RMSE)
mean( eval_gwr$RMSE)
mean(eval_ols$MAE)
mean( eval_grf$MAE)
mean( eval_grf_lm_050$MAE)
mean( eval_gwr$MAE)

