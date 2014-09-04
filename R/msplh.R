# [This is from Mary Meyer's website: http://www.stat.colostate.edu/~meyer/srrs.htm]

# monotone regression spline code: hinge algorithm in R
# must have x=predictor,y=response.
#  can specify either k=#interior knots
#  (will place at equal x-quantiles)
#  or knots themselves.  If latter, must have knots=TRUE
# quadratic version 
# returns predicted values in yhat, knots, and vectors 
# for plotting 

mspl=function(x,y,kt,knots){
	n=length(x)
	if(knots){
		t=kt
		k=length(t)-2
	}else{
		xs=sort(x)
		k=kt
		knots=round((0:(k+1))*n/(k+1))
		knots[1]=1
		t=xs[knots]
	}
	m=k+3
	ans1=monincr(x,t)
	xpl=ans1$xpl
	spl=ans1$spl
	edges=ans1$sigma
	amat=matrix(0,nrow=m,ncol=m)
	for(i in 2:m){amat[i,i]=1}
	qmat=edges%*%t(edges)
	cvec=edges%*%y
	coef=quadprog(cvec,qmat,amat)
	yhat=t(edges)%*%coef
	yhatp=t(spl)%*%coef
	ans=new.env()
	ans$yhat=yhat
	ans$ypl=yhatp
	ans$xpl=xpl
	ans$knots=t
	ans
}

##############################################################
# Monotone increasing
##############################################################
monincr=function(x,t){
	n=length(x)
	k=length(t)-2
	m=k+3
	sigma=matrix(1:m*n,nrow=m,ncol=n)
	sigma[1,]=1:n*0+1
	for(j in 1:(k-1)){
	 	i1=x<=t[j]
	 	sigma[j+1,i1] = 0
	 	i2=x>t[j]&x<=t[j+1]
		sigma[j+1,i2] = (x[i2]-t[j])^2 / (t[j+2]-t[j]) / (t[j+1]-t[j])
	    i3=x>t[j+1]&x<=t[j+2]
		sigma[j+1,i3] = 1-(x[i3]-t[j+2])^2/(t[j+2]-t[j+1])/(t[j+2]-t[j])
	    i4=x>t[j+2]
		sigma[j+1,i4]=1
	}

	i1=x<=t[k]
	sigma[k+1,i1] = 0
	i2=x>t[k]&x<=t[k+1]
	sigma[k+1,i2] = (x[i2]-t[k])^2 / (t[k+2]-t[k]) / (t[k+1]-t[k])
	i3=x>t[k+1]&x<=t[k+2]
	sigma[k+1,i3] = 1- (x[i3]-t[k+2])^2/(t[k+2]-t[k+1])/(t[k+2]-t[k])
	i4=x>t[k+2]
	sigma[k+1,i4]=1
	
	i1=x<=t[2]
	sigma[k+2,i1]=1-(t[2]-x[i1])^2/(t[2]-t[1])^2
	i2=x>t[2]
	sigma[k+2,i2]=1
	
	i1=x<=t[k+1]
	sigma[k+3,i1]=0
	i2=x>t[k+1]&x<=t[k+2]
	sigma[k+3,i2]=(x[i2]-t[k+1])^2/(t[k+2]-t[k+1])^2
	i3=x>t[k+2]
	sigma[k+3,i3]=1

	spl=matrix(1:m*n,nrow=m,ncol=101)
	xpl=0:100/100*(max(x)-min(x))+min(x)
	spl[1,]=1:101*0+1
	for(j in 1:(k-1)){
	 	i1=xpl<=t[j]
	 	spl[j+1,i1] = 0
	 	i2=xpl>t[j]&xpl<=t[j+1]
		spl[j+1,i2] = (xpl[i2]-t[j])^2 / (t[j+2]-t[j]) / (t[j+1]-t[j])
	    i3=xpl>t[j+1]&xpl<=t[j+2]
		spl[j+1,i3] = 1-(xpl[i3]-t[j+2])^2/(t[j+2]-t[j+1])/(t[j+2]-t[j])
	    i4=xpl>t[j+2]
		spl[j+1,i4]=1
	}

	i1=xpl<=t[k]
	spl[k+1,i1] = 0
	i2=xpl>t[k]&xpl<=t[k+1]
	spl[k+1,i2] = (xpl[i2]-t[k])^2 / (t[k+2]-t[k]) / (t[k+1]-t[k])
	i3=xpl>t[k+1]&xpl<=t[k+2]
	spl[k+1,i3] = 1- (xpl[i3]-t[k+2])^2/(t[k+2]-t[k+1])/(t[k+2]-t[k])
	i4=xpl>t[k+2]
	spl[k+1,i4]=1
	
	i1=xpl<=t[2]
	spl[k+2,i1]=1-(t[2]-xpl[i1])^2/(t[2]-t[1])^2
	i2=xpl>t[2]
	spl[k+2,i2]=1
	
	i1=xpl<=t[k+1]
	spl[k+3,i1]=0
	i2=xpl>t[k+1]&xpl<=t[k+2]
	spl[k+3,i2]=(xpl[i2]-t[k+1])^2/(t[k+2]-t[k+1])^2
	i3=xpl>t[k+2]
	spl[k+3,i3]=1
	for(i in 2:m){
		ms=mean(sigma[i,])
		spl[i,]=spl[i,]-ms
		sigma[i,]=sigma[i,]-ms
	}
	ans=new.env()
	ans$sigma=sigma
	ans$xpl=xpl
	ans$spl=spl
	ans
}
##############################################################
# quadratic programming code: hinge algorithm in R
#  find theta to minimize t(theta)%*%qmat%*%theta-2*theta%*%cvec
#  subject to amat%*%theta >= 0
#
# qmat must be positive definite, amat must be irredicible
#
quadprog=function(cvec,qmat,amat){
	n=length(cvec)
	m=length(amat)/n
	sm=1e-10;h=1:m<0;obs=1:m;check=0
	umat=chol(qmat);uinv=solve(umat)
	delta=-amat%*%uinv
	y=t(uinv)%*%cvec
	b2=delta%*%y
	if(max(b2)>sm){
		i=min(obs[b2==max(b2)])
		h[i]=TRUE
	}else{check=1;theta=1:n*0}
	while(check==0){
		xmat=matrix(delta[h,],ncol=n)
		a=solve(xmat%*%t(xmat))%*%xmat%*%y
		if(min(a)<(-sm)){
			avec=1:m*0;avec[h]=a
			i=min(obs[avec==min(avec)])
			h[i]=FALSE;check=0
		}else{
			check=1
			theta=t(xmat)%*%a
			b2=delta%*%(y-theta)
			if(max(b2)>sm){
				i=min(obs[b2==max(b2)])		
				h[i]=TRUE;check=0
			}
		}
	}
	bhat=uinv%*%(y-theta)
	bhat
}
