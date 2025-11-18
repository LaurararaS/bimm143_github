source("http://thegrantlab.org/misc/cdc.R")
head(cdc$height)
tail(cdc$weight, n=20)

install.packages("tideverse")
library(ggplot2)
ggplot(data=cdc, aes(height, weight)) + geom_point() + labs(x="Height (inches)",y="Weight (pounds)")
cor(cdc$height, cdc$weight)
hist(cdc$height)

height_m <- cdc$height * 0.0254
weight_kg <- cdc$weight * 0.454
bmi = (weight_kg)/(height_m^2)
plot(cdc$height, bmi)
cor(cdc$height, bmi)

sum(bmi >= 30)
round( (sum(bmi >= 30)/length(bmi)) * 100, 1)
plot(cdc$height[1:100], cdc$weight[1:100])
plot(cdc[1:100,5], cdc[1:100,6])


table( c(cdc[,9]))

num=0
gender = cdc[,9]
for (i in 1:20000){
  if(gender[i]=="m" & bmi[i]>=30){
    num=num+1
  }
}
num