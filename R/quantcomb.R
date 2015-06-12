

### Code for calculating quantile frequencies ###

#Count frequencies for a given group
bin.exprs = function(grp.exprs,quantiles,labels) {
	#print(grp.exprs)
	#print(quantiles)

	#intialize matrix to store adjusted expression values
	values = matrix(0,ncol=1,nrow=length(grp.exprs))

	num.quantiles = length(quantiles)

	values[grp.exprs < quantiles[1]] = labels[1]

	for(index in c(2:num.quantiles)) {
		values[grp.exprs >= quantiles[index-1] & grp.exprs < quantiles[index] ] = labels[index]
	}

	values[grp.exprs >= quantiles[num.quantiles]] = labels[(num.quantiles+1)]

	#print(values)

	#function to calculate frequencies of quantile labels
	freqs = function(quantile.label) { sum(values == quantile.label)}

	#return frequencies as a vector
	#return(sapply(-2:2,freqs))
	return(sapply(labels,freqs))
}

#Setup quantiles and groups and get group frequencies
get.quantile.freqs = function(exprs,grp1,grp2,labels) {

	#add check for number of quantiles

	num.quantiles = length(labels) - 1

	quantile.unit = 1/(num.quantiles+1)
	quantile.thresholds = sapply(1:num.quantiles,function(x){quantile.unit*x})
	#print(quantile.thresholds)

	#calculate quantiles for gene
	quantiles = quantile(exprs,quantile.thresholds)
	#print(quantiles)

	grp1.values = bin.exprs(exprs[grp1],quantiles,labels)
	grp2.values = bin.exprs(exprs[grp2],quantiles,labels)


	return(cbind(grp1.values,grp2.values))
	#return(list(grp1=grp1.values,grp2=grp2.values))
}

### ###


### Code for performing Disco-Normal Analysis ###

#Calculate first and second derivative
derivatives = function(fr,x,beta,k) {

	nrowx = length(x[,1])

 	sampln = length(fr[1,])

	nrowb = length(beta)

	n = colSums(fr)

 	xb1 = fr*x

 	xb = colSums(xb1)

 	xs1 = fr*(x*x)

 	xs = colSums(xs1)

 	nrowxs = length(xs)



 	nr=matrix(0,ncol=sampln, nrow=nrowx)
 	ar11= matrix(0,ncol=sampln, nrow=nrowx)
 	ar21= matrix(0,ncol=sampln, nrow=nrowx)
 	ar31= matrix(0,ncol=sampln, nrow=nrowx)
 	ar41= matrix(0,ncol=sampln, nrow=nrowx)

 	anr=matrix(0,ncol=sampln, nrow=nrowx)
 	aar11= matrix(0,ncol=sampln, nrow=nrowx)
 	aar21= matrix(0,ncol=sampln, nrow=nrowx)
 	aar31= matrix(0,ncol=sampln, nrow=nrowx)
 	aar41= matrix(0,ncol=sampln, nrow=nrowx)

 	f1= matrix(0,ncol=1, nrow=nrowb)
 	f2= matrix(0,ncol=1, nrow=3)
 	d1= matrix(0,ncol=nrowb, nrow=nrowb)
 	d2= matrix(0,ncol=3, nrow=3)


	for (j in 1:nrowx) {

		for (i in 1:sampln) {

  			a12= cbind(x[j,i], x[j,i]^2)

			if (k==1) {

 				beta2= rbind(beta[i], beta[sampln+i])

 				nr[j,i]=exp(a12%*%beta2)

 				ar11[j,i]=x[j,i]*nr[j,i]

 				ar21[j,i]=(x[j,i]^2)*nr[j,i]

 				ar31[j,i]=(x[j,i]^3)*nr[j,i]

 				ar41[j,i]=(x[j,i]^4)*nr[j,i]

			}

			if (k==2) {

				beta2=rbind(beta[1], beta[i+1])

				anr[j,i]=exp(a12%*%beta2)

				aar11[j,i]=x[j,i]*anr[j,i]

				aar21[j,i]=(x[j,i]^2)*anr[j,i]

				aar31[j,i]=(x[j,i]^3)*anr[j,i]

				aar41[j,i]=(x[j,i]^4)*anr[j,i]

			}
 		}
 	}

 	ar= colSums (nr)
 	aar= colSums (anr)

 	ar1= colSums(ar11)
 	ar2= colSums(ar21)
 	ar3= colSums(ar31)
 	ar4= colSums(ar41)

 	aar1= colSums(aar11)
 	aar2= colSums(aar21)
 	aar3= colSums(aar31)
 	aar4= colSums(aar41)

 	ac1=matrix(0, ncol=1, nrow=sampln)
 	ac2=matrix(0, ncol=1, nrow=sampln)

 	acc1=0
 	acc2=0

 	xbb=0

 	for (l in 1:sampln) {

 		if (k==1) {

 			f1[l]=xb[l]-n[l]*(ar1[l]/ar[l])

			f1[l+sampln]=xs[l]-n[l]*(ar2[l]/ar[l])

			d1[l, l]=-n[l]*(((-ar1[l]^2)+(ar[l]*ar2[l]))/(ar[l]^2))

			d1[l, sampln+l]=-n[l]*(((-ar2[l]*ar1[l])+(ar[l]*ar3[l]))/(ar[l]^2))

 			d1[sampln+l,l]=d1[l,sampln+l]

 			d1[l+sampln, l+sampln]=-n[l]*(((-ar2[l]^2)+(ar[l]*ar4[l]))/(ar[l]^2))
 		}

		if (k==2) {

			ac1[l]=n[l]*(aar1[l]/aar[l])

			ac2[l]=n[l]*(((-aar1[l]^2)+(aar[l]*aar2[l])))/(aar[l]^2)

			acc1=acc1+ac1[l]

 			acc2=acc2+ac2[l]

 			f2[l+1]=xs[l]-n[l]*(aar2[l]/aar[l])

 			d2[l+1,l+1]=-n[l]*(((-aar2[l]^2)+(aar[l]*aar4[l]))/(aar[l]^2))

 			d2[1, l+1]=-n[l]*(((-aar2[l]*aar1[l])+(aar[l]*aar3[l]))/(aar[l]^2))

 			d2[l+1,1]=d2[1, l+1]

 			xbb=xbb+xb[l]
 		}
	}

	if(k==1) {
		f=f1
		d=d1
	}


	if (k==2) {

		f2[1]=xbb-acc1

		d2[1, 1] = -acc2

		f=f2

		d=d2
	}

	return(list(f,d))
}

#Optimize function
newton=function(f,d,beta, tol,fr,x,k) {

	if (k==2) {
		beta=matrix(0, ncol=1, nrow=3)
	}

	while(abs(f) > tol ) {

		f3=derivatives(fr=fr,x=x,beta=beta,k=k)

		f=f3[[1]]

		d=f3[[2]]
		#print(d)
		#print(f)
		delta=solve(-1*d)%*%f
		#print(delta)
		#print(beta)
		beta=beta+delta
	}

	return(beta)
}

#Computer Chi-Square and p-value for data
disco.chisq = function(fr,x) {

	beta=matrix(0, ncol=1, nrow=4)

	f3=derivatives(fr,x,beta,k=1)

	f=f3[[1]]

	d=f3[[2]]
        #print(f)
	beta1=newton(f,d, beta, tol=0.000000000001,fr,x,k=1)

	#print(beta1)
	beta2=newton(f,d, beta, tol=0.000000000001,fr,x,k=2)
        #print(beta2)
	logl=matrix(0,ncol=2, nrow=1)

	sampln=length(fr[1,])

	n1=colSums(fr)

	xb1=fr*x

	xb= colSums (xb1)

	xs1=fr*(x*x)

	xs=colSums(xs1)

	nrowx=length(x[,1])

	nr=matrix(0,ncol=sampln, nrow=nrowx)

	larp1=matrix(0,ncol=1, nrow=2)


	for (k in 1:2) {
		if (k==1) {
			alpha=beta1[1:2]
			bes=beta1[3:4]
		}

		if (k==2) {
			alpha=matrix(beta2[1], ncol=1, nrow=2)
			bes=beta2[2:3]
		}

		for (j in 1:nrowx) {

			for (i in 1:sampln) {

 				a12= cbind(x[j,i], x[j,i]^2)

 				beta2a= rbind(alpha[i], bes[i])

 				nr[j,i]=exp(a12%*%beta2a)
	 		}
	 	}

		br=colSums(nr)

		larp1[k]=n1%*%log(br)

		logl[k]=xb%*%alpha+xs%*%bes-larp1[k]
	}

	chisq=-2*(logl[2]-logl[1])

	pvalue=1- pchisq(chisq, 1)

	return(list("beta1"=beta1,"beta2"=beta2,"chisq"=chisq,"pvalue"=pvalue))
	#c(beta1,beta2,chisq,pvalue)
}

### ###


### Code for performing analysis on a data set ###

#Calculate p-value using t-test
get.pval.ttest = function(dataf,index1,index2,datafilter=as.numeric) {

	f = function(i) {
		return(t.test(datafilter(dataf[i,index1]),datafilter(dataf[i,index2]))$p.value)
	}

	return(sapply(1:length(dataf[,1]),f))
}

#Calculate p-value using t-test
get.tstat.ttest = function(dataf,index1,index2,datafilter=as.numeric) {

	f = function(i) {
		return(as.numeric(t.test(datafilter(dataf[i,index1]),datafilter(dataf[i,index2]))$statistic))
	}

	return(sapply(1:length(dataf[,1]),f))
}

#Calculate p-value using Wilcoxon
get.pval.wilcox = function(dataf,index1,index2,datafilter=as.numeric) {

	f = function(i) {
		return(wilcox.test(datafilter(dataf[i,index1]),datafilter(dataf[i,index2]))$p.value)
	}

	return(sapply(1:length(dataf[,1]),f))
}

#Calculate fold change
fold.change = function(dataf,grp1,grp2) {
library(gtools)

	f = function(i) {
		#return(mean(dataf[i,grp1]) - mean(dataf[i,grp2]))
		a = mean(unlist(dataf[i,grp1]))
		b = mean(unlist(dataf[i,grp2]))
		#a = mean(dataf[i,grp1])
		#b = mean(dataf[i,grp2])

		return(foldchange(a,b))
	}
	return(sapply(1:nrow(dataf),f))
}

#Obtain Chi-Square and p-value for each gene
get.disco.pvalues = function(gene.exprs,grp1,grp2) {
	#gene.exprs = as.matrix(gene.exprs)



	f = function(i) {
		x=cbind(c(-2,-1,0,1,2), c(-2,-1,0,1,2))
		#k=2

		fr = get.quantile.freqs(as.numeric(gene.exprs[i,]),grp1,grp2)
		#fr
		disco.chisq(fr,x,length(c(grp1,grp2)))$pvalue
	}
	return(sapply(1:nrow(gene.exprs),f))
}


#Obtain Chi-Square and p-value for each gene
get.disco.res = function(gene.exprs,num.quantiles,grp1,grp2) {


	num.labels = num.quantiles + 1
	seq = seq(length.out = num.labels)
	labels = seq - median(seq)

	#if num.labels is odd, will get half values, so add 0.5
	if(round(labels[1]) != labels[1]) {
		labels = labels + 0.5
	}


	f = function(i) {
		#x=cbind(c(-2,-1,0,1,2), c(-2,-1,0,1,2))
		#k=2

		fr = get.quantile.freqs(exprs=as.numeric(gene.exprs[i,]),grp1=grp1,grp2=grp2,labels=labels)
		#fr = get.quantile.freqs(as.numeric(gene.exprs[i,]),grp1,grp2)
		#fr
		disco.chisq(fr=fr[,],x=cbind(labels,labels))
		#disco.chisq(fr,x,length(c(grp1,grp2)))
		#return(1)
	}
	res = lapply(1:nrow(gene.exprs),f)

        data.frame(
                   t(sapply(1:length(res),
                          function(g) {
                              c(
                                "p.value"=res[[g]]$pvalue,
                                "chisq"=res[[g]]$chisq,
                                "beta1.1"=res[[g]]$beta1[1],
                                "beta1.2"=res[[g]]$beta1[2],
                                "beta1.3"=res[[g]]$beta1[3],
                                "beta1.4"=res[[g]]$beta1[4],
                                "beta2.1"=res[[g]]$beta2[1],
                                "beta2.2"=res[[g]]$beta2[2],
                                "beta2.3"=res[[g]]$beta2[3]
                                )
                          }
                   )),
                   row.names=row.names(gene.exprs)
                   )

}

#Calculate p-value for Shapiro-Wilk test
get.pval.shapiro1 = function(dataf,index,datafilter=as.numeric) {

	f = function(i) {
		return(shapiro.test(datafilter(dataf[i,index]))$p.value)
	}

	return(sapply(1:length(dataf[,1]),f))
}



