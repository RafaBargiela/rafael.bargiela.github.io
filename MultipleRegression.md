# MULTIPLE REGRESSION IN R

Description about how to perform a multiple regression in R to analyse the effect of two or more independent variables over our dataset (dependent variable). Mainly, this will be perform by Multiple Linear Regression, but we will see also how to perform non-Linear Multiple Regression like Multiple Polynomial Regression.

## 1. Multiple Linear Regression (MLR)

The formula for MLR is as follows:

$$
	Y= \beta_0 + \beta_1X_1 + ... + \beta_nX_n + \epsilon 
$$

Where _Y_ is the predicted value, \\(\beta_0\\) is the Intercept value, meaning the value of _Y_ when the other parameters are 0. \\(\beta_1X_1\\) is the regression coefficient of the first independent variable, and \\(\beta_nX_n\\) the regression coefficient of the last independent variable. \\(\epsilon\\) is the model error, or how much variation there is from our stimation of _Y_ from the real data.

To perform a MLR in are we first will calculate a linear model using all the independent variables which we want to initially measure over our dataset. Then, we will estimate which of these varaibles are really important and get a final model including only these variables.

### 1.1 Developing the model

Our initial model will include all variables in out dataset. To do that we will use the funcion _lm_ from base R. We will use as example data with relative abundance of Bacteria and concentration of different metals. To simply add all variables in the table, we can use \\(.\\) after adding one or two variables:

```R
Initial.mod<-lm(Bacteria~V+., data=dataset.DF)
```

With _stepAIC_ from package _MASS_ we can evaluate the acuracy of our model and fetch the best fitting model:

```R
step<-stepAIC(Initial.mod,direction="both")
final.mod<-lm(step$model)
```

We can check parameters of our final model with _summary_ and also fetch coefficients values, \\(R²\\) and p-value:

```R
	summary(final.mod)
	# Coefficients
	coef<-summary(final.fit)$coefficients
	# R-squared
	R<-format(as.numeric(summary(final.fit)$adj.r.squared),nsmall=3,digits=2)
	# P-value
	f<-summary(final.fit)$fstatistic
	pvmod <- format(as.numeric(pf(f[1],f[2],f[3],lower.tail=F)),nsmall=3,digits=2)
```

### 1.2 Evaluating the model through plots
There are different plots which could help to asses the adjustment of our model. On one hand, we have the Q-Q plot, to evaluate the general fit of the whole model. Then, if we want to evaluate the influence of each of the variables included on it, we can make a plot with te estimates values and its Standard Error, and also que can perform Added-plots for each variable, based on the residuals of each independent Linear Regression over the residuals of the Linear regression of each specific variable over the rest variables.

#### 1.3 Q-Q plot
Quantile-Quantile plot is a good way to evaluate how the model fits to the data. Is a good place to print the R² and p-value. First, we need to get the standarized residuals of our final model and sort them on increasing order:

```R
res.std<-rstandard(final.fit)
res.std.sort<-sort(res.std)
```
Then, we create normal distributed random data of same length than our standarized residuals vector:
```R
i<-1:length(rst.sort)
fi<-(i-0.5)/length(rst.sort)
x.norm<-qnorm(fi) 
```
