
ssq<-function(b){sum(b**2)} #ssq

#projorth
projorth<-function(x, M)
{
  Mtx<-t(M)%*%x
  MMtx<-M%*%Mtx
  x-MMtx
}

#projL1_2L2.orth
projL1_2L2.orth<-function(x, tau, alpha=0.0000000001, M, itermax=1000,eps=1e-16){
  xold <- xnew <- x
  
  for (k in 1:itermax) {
    if (is.null(M)) 
    {
      xnew <- projL1_2L2(as.numeric(xold), tau=tau, alpha = alpha)$x
      
      while (length(xnew)<1){ xnew <- projL1_2L2(as.numeric(xold), tau=tau, alpha = alpha)$x}
      # xnew <- projl1(normalize(xold), a=a)
      # xnew <- 1/2 * ( projl1(xold, a=a) + normalize(xold) )
    } 
    else 
    {
      xnew <- projL1_2L2(as.numeric(projorth(xold, M)), tau=tau, alpha=alpha)$x
      
      while (length(xnew)<1){ xnew <- projL1_2L2(as.numeric(projorth(xold, M)), tau=tau, alpha=alpha)$x}
      
      # xnew <- projl1(projorth(normalize(xold), M), a=a)
      # xnew <- 1/3 * ( projl1(xold, a=a) + normalize(xold) + projorth(xold, M) )
    }
    
    if ( norm2(xnew - xold) < eps ) 
    {
      xold <- xnew
      return(list(x=xnew, k=k))
      break
    }
    xold <- xnew
  }
  
  return(list(x=xnew, k=k))
  
}



#projL1_2L2
projL1_2L2<-function(x, tau, alpha=0.0000000001){
  
  norm2_x <- norm2(x)
  if( norm2_x < 1e-32 ) {return(list(x=x))} #Comprobaci???n de que el denominador no es 0
  
  if( (((1-alpha)*norm1(x)+(alpha)*norm2_x^2)/norm2_x) <= tau ) 
  {
    return(list(x=x/norm2_x))
  }
  
  uneq<-x != 0
  p <- abs(x[x != 0])
  
  
  MAX=max(p)
  bMAX<-p==MAX
  nMAX<-sum(bMAX)
  
  
  
  if(tau>((1-alpha)*sqrt(length(x))+alpha))
    stop("Impossible to project, maximum ratio is : ", ((1-alpha)*sqrt(length(x))+alpha*norm2(x)))
  if(tau<=0) stop("Impossible to project, minimum ratio is : ",0)
  
  s1 <-0
  s2 <-0
  nb <- 0
  
  while (T) {
    N <- length(p)
    
    
    a_k <- p[sample(1:N,1)]
    while(a_k == MAX){
      a_k <- p[sample(1:N,1)]
    }
    
    
    p_inf_ak <- p < a_k
    p_sup_ak <- p > a_k
    p_high <- p[p_inf_ak]
    p_low <- p[p_sup_ak]
    
    
    nb_a_k <- sum(p == a_k)
    k <- nb + sum(p_sup_ak) + nb_a_k
    
    
    aksq <- a_k^2
    slow_1 <- sum(p_low) + nb_a_k*a_k
    slow_2 <- ssq(p_low) + nb_a_k*aksq
    omega_a_k<-(alpha*(s2 + slow_2) +(1-alpha)*(s1 + slow_1) - k*a_k*((1-alpha)^2) - k*(a_k^2)*alpha*(1-alpha)^2)/((1+2*alpha*a_k)*sqrt(s2 + slow_2 - 2*a_k*(1-alpha)*(s1 + slow_1) + k*aksq*(1-alpha)^2))
    
    if (omega_a_k > tau){
      if ( length(p_low) <= 1 ) break #length(p_low) == 0 
      p <- p_low
    }else{
      if (length(p_high) <=1){   #length(p_low) == 0 
        break
      }else{
        a_k_1 <- max(p_high)
        omega_a_k_1<-(alpha*(s2 + slow_2) +(1-alpha)*(s1 + slow_1) - k*a_k_1*((1-alpha)^2) - k*(a_k_1^2)*alpha*(1-alpha)^2)/((1+2*alpha*a_k_1)*sqrt(s2 + slow_2 - 2*a_k_1*(1-alpha)*(s1 + slow_1) + k*a_k_1^2*(1-alpha)^2))
        if (omega_a_k_1 >tau){
          break
        }
        p <- p_high
        nb <- k
        s1 <- s1 + slow_1
        s2 <- s2 + slow_2
      }
    }
  }
  
  
  l1<-s1 + slow_1
  l2<-s2 + slow_2
  a<-(tau^2)*l2-(alpha^2)*(l2^2)-2*alpha*l2*(1-alpha)*l1-(l1^2)*((1-alpha)^2) #termino independiente
  b<- 4*alpha*l2*(tau^2)-2*(tau^2)*(1-alpha)*l1+2*l2*k*((1-alpha)^2)*alpha+2*l1*k*((1-alpha)^3) #Coeficiente grado 1
  c<- k*(tau^2)*((1-alpha)^2)+4*(alpha^2)*l2*(tau^2)-8*alpha*(tau^2)*(1-alpha)*l1+2*k*l1*alpha*((1-alpha)^3)+2*k*l2*(alpha^2)*((1-alpha)^2)-(k^2)*((1-alpha)^4) #Coeficiente grado 2
  d<-4*k*alpha*(tau^2)*((1-alpha)^2)-8*alpha^2*(tau^2)*(1-alpha)*l1-2*(k^2)*alpha*((1-alpha)^4) #Coeficiente grado 3
  e<-4*k*(alpha^2)*(tau^2)*((1-alpha)^2)-(k^2)*(alpha^2)*((1-alpha)^4) #Coeficiente grado 4
  sol<-Re(polyroot(c(a,b,c,d,e))) #Resolucion de la ecuacion de grado 4
  sol<-sol[which(sol>=0)]
  
  
  lambda<-sol[which(sol>=min(abs(x[uneq]))&sol<=max(abs(x[uneq])))] #Seleccion de la solucion valida. Pertenece al intervalo [abs(x1),abs(xJ)]
  if(length(lambda)>1){lambda<-min(lambda)} #Si existe mas de una valida, seleccion de la minima
  if(length(lambda)<1){lambda=0}
  
  x_en <- (sign(x)*pmax(0, abs(x) - lambda*(1-alpha)))/(1+2*lambda*alpha) 
  return( list(x=x_en / norm2(x_en) , penalizacion=lambda, sol=sol) ) 
  
}


#powit.enet2
powit.enet2<-function(X,U0,V0,Uorth,Vorth,tau.u=1.4,tau.v=1.4,
         alpha.u=0.000001, alpha.v=0.000001,
         eps.pi=0.00001,eps.pocs=0.00001, 
         itermax.pi=1000, itermax.pocs=1000)
{
  #alpha.v<-alpha.v
  #alpha.u<-alpha.u
  tau.u<-tau.u
  tau.v<-tau.v
  uold <- unew <- U0
  vold <- vnew <- V0
  for (iter in 1:itermax.pi) {
    
    vnew <- projL1_2L2.orth(x=(t(X) %*% uold)[,1], tau=tau.v, alpha=alpha.v,M=Vorth,
                            itermax = itermax.pocs, eps = eps.pocs)$x
    
    #while(length(vnew)<1){vnew <- projL1_2L2.orth(x=(t(X) %*% uold)[,1], tau=tau.v, alpha=alpha.v,M=Vorth,itermax = itermax.pocs, eps = eps.pocs)$x}
    unew <- projL1_2L2.orth(x=(X%*%vnew)[,1], tau=tau.u, alpha=alpha.u, M=Uorth,
                            itermax = itermax.pocs, eps = eps.pocs)$x
    #while(length(unew)<1){ unew <- projL1_2L2.orth(x=c(X %*% vnew)[,1], tau=tau.u, alpha=alpha.u, M=Uorth, itermax = itermax.pocs, eps = eps.pocs)$x}
    if ( norm2(vnew - vold) < eps.pi && norm2(unew - uold) < eps.pi ) 
    {
      vold <- vnew
      uold <- unew
      return(list(U=unew, V=vnew, iter=iter))
      break
    }
    vold <- vnew
    uold <- unew
    
  }
  return(list(U=unew, V=vnew, iter=iter))
}

#powit.enet1
powit.enet1<-function(X,U0,V0,tau.u,tau.v,alpha.u, alpha.v,eps.pi, itermax.pi)
{
  uold <- unew <- U0
  vold <- vnew <- V0
  
  for (iter in 1:itermax.pi) {
    vnew <- projL1_2L2(x=(t(X)%*% uold)[,1], tau=tau.v, alpha=alpha.v)$x
    #while(length(vnew)<1){vnew <- projL1_2L2(x=(t(X)%*% uold)[,1], tau=tau.v, alpha=alpha.v)$x}
    
    unew <- projL1_2L2(x=(X %*% vnew)[,1], tau=tau.u, alpha=alpha.u)$x
    #while(length(unew)<1){unew <- projL1_2L2(x=(X %*% vnew)[,1], tau=tau.u, alpha=alpha.u)$x}
    if ( norm2(vnew - vold) < eps.pi && norm2(unew - uold) < eps.pi ) 
    {
      
      vold <- vnew
      uold <- unew
      return(list(U=unew, V=vnew, iter=iter))
      break
    }
    
    vold <- vnew
    uold <- unew
    
  }
  return(list(U=unew, V=vnew, iter=iter))
}

#pca.enet
pca.enet<-function(X, Q=2, tau.u = 1.4, tau.v = 1.4,alpha.u=1e-16, alpha.v=1e-16,
         itermax.pi=1000, itermax.pocs=1000,eps.pi=1e-16, eps.pocs=1e-16, 
         init.svd="svd", init.transf=1,obs.names=FALSE,plot.axis=c(1,2))
{
  library(ggplot2)
  library(reshape2)
  varnames<-colnames(X)
  obsnames<-row.names(X)
  
  X<-init.transformation(X,init.transf)
  n.comp<-Q
  res.csvd<-csvd.enet(X, Q=n.comp, tau.u = rep(tau.u,n.comp), tau.v = rep(tau.v,n.comp),alpha.u = alpha.u, alpha.v=alpha.v,
                      itermax.pi=itermax.pi, itermax.pocs=itermax.pocs,eps.pi=eps.pi, eps.pocs=eps.pocs, 
                      init.svd=init.svd)
  
  V.matrix<-(res.csvd$V)[,1:n.comp]
  row.names(V.matrix)<-varnames
  char_vector <- character(dim(V.matrix)[2])
  for(i in 1:dim(V.matrix)[2]){char_vector[i]<-paste("PC",i,sep="")}
  colnames(V.matrix)<-char_vector
  D.matrix<-res.csvd$D
  variance.proportion<-D.matrix[1:Q]^2/(norm(X,"f")^2)*100
  scores.matrix<-X%*%V.matrix
  colnames(scores.matrix)<-char_vector
  row.names(scores.matrix)<-obsnames
  
  #Loadings plot
  
  df.loadings<-matrix(0,nrow=dim(V.matrix)[1]*dim(V.matrix)[2],ncol=3)
  df.loadings[,2]<-as.numeric(V.matrix)
  df.loadings<-as.data.frame(df.loadings)
  df.loadings[,1]<-c(rep(row.names(V.matrix),dim(V.matrix)[2]))
  labels.pcs<-c()
  
  for(k in 1:dim(V.matrix)[2])
  {
    labels.pcs<-c(labels.pcs,rep(colnames(V.matrix)[k],dim(V.matrix)[1]))
  }
  df.loadings$V3<-labels.pcs
  class(df.loadings)
  p=ggplot(df.loadings, aes(V3, V1)) +
    geom_tile(aes(fill = V2), color = "white") +
    scale_fill_gradient2(low = ("violetred4"), mid = "white",
                         high = ("steelblue"), midpoint = 0)+
    ylab("variables ") +
    xlab("PCs") +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 12),
          plot.title = element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(fill = "Loadings")
  
  #Scores plot
  
  i=plot.axis[1]
  j=plot.axis[2]
  df.scores<- as.data.frame(scores.matrix)
  
  colnames(df.scores)<-colnames(scores.matrix)
  theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))
  if(obs.names)
  {
    s<-ggplot(df.scores,aes(x=df.scores[,i],y=df.scores[,j],label=row.names(df.scores) ))
    s<-s+geom_point()+ geom_text(size=2.5)+theme+
      geom_hline(yintercept=0, color="gray", size=0.5)+
      geom_vline(xintercept=0, color="gray", size=0.5)+
      xlab(colnames(df.scores)[i]) + ylab(colnames(df.scores)[j])
  }
  else
  {
    s<-ggplot(df.scores,aes(x=df.scores[,i],y=df.scores[,j]))+
      geom_point()+theme+
      geom_hline(yintercept=0, color="gray", size=0.5)+
      geom_vline(xintercept=0, color="gray", size=0.5)+
      xlab(colnames(df.scores)[i]) + ylab(colnames(df.scores)[j])
    
  }
  
  return(list(Scores=scores.matrix,
              Loadings=V.matrix,
              Variance=variance.proportion,
              Loadings.plot=p, Scores.plot=s))
}


#normalize
normalize<-function(u) 
{
  norm.u<-norm2(u)
  if(norm.u==0) {return(u)}
  else return(u/norm.u)
}


#norm1, norm2
norm2<-function(b){sqrt(sum(b^2))}
norm1<-function(b){sum(abs(b))}



#init.transformation
init.transformation<-function(X,t)
{
  varnames<-colnames(X)
  obsnames<-row.names(X)
  
  X<-as.matrix(X)
  row.names(X)<-obsnames
  colnames(X)<-varnames
  if(t==1)return(X)#Raw Data: When no transformation is required.
  if(t==2)return(X-mean(X))#Substract the global mean: Eliminate an eefect common to all the observations
  #Column centering (Remove the column means):
  if(t==3){
    Y<-matrix(0,nrow=dim(X)[1],ncol=dim(X)[2])
    for(j in 1:(dim(X)[2])){
      Y[,j]<-X[,j]-mean(X[,j])
    }
    row.names(Y)<-obsnames
    colnames(Y)<-varnames
    return(Y)}
  #Standardize columns (Remove the column means and divide by its standard deviation):
  if(t==4){
    Y<-matrix(0,nrow=dim(X)[1],ncol=dim(X)[2])
    for(j in 1:(dim(X)[2])){
      Y[,j]<-(X[,j]-mean(X[,j]))/sd(X[,j])
    }
    row.names(Y)<-obsnames
    colnames(Y)<-varnames
    return(Y)}
  #Row centering (Remove the row means)
  if(t==5){
    Y<-matrix(0,nrow=dim(X)[1],ncol=dim(X)[2])
    for(i in 1:(dim(X)[1])){
      Y[i,]<-X[i,]-mean(X[i,])
    }
    row.names(Y)<-obsnames
    colnames(Y)<-varnames
    return(Y)}
  #Standardize rows (Divide each row by its standard deviation):
  if(t==6){
    Y<-matrix(0,nrow=dim(X)[1],ncol=dim(X)[2])
    for(i in 1:(dim(X)[1])){
      Y[i,]<-(X[i,]-mean(X[i,]))/sd(X[i,])
    }
    row.names(Y)<-obsnames
    colnames(Y)<-varnames
    return(Y)} 
  #Double centering (Interaction residuals): When all the elements of the table are comparable. Useful for AMMI models
  if(t==7){
    Y<-matrix(0,nrow=dim(X)[1],ncol=dim(X)[2])
    #meanR<-rowMeans(X)
    for(j in 1:(dim(X)[2]))
    {
      for(i in 1:(dim(X)[1]))
      {
        Y[i,j]<-X[i,j]-(rowMeans(X)[i])-(colMeans(X)[j])+mean(X)
      }
    }
    row.names(Y)<-obsnames
    colnames(Y)<-varnames
    return(Y)}
}

#cv.alpha_enet
cv.alpha_enet<-function(data, Q, set.alpha.v=c(),nfolds=10,parallel=FALSE, type.measure="MSE", 
         tau.u=sqrt(dim(data)[2]), tau.v=1.4,alpha.u=0.000001,ntau=5,
         itermax.pi=500, itermax.pocs=500,eps.pi=1e-16, eps.pocs=1e-16, 
         init.svd="svd")
{
  library(ggplot2)
  library(caret)#CreateFolds
  
  data<-as.matrix(data)
  ncol<-dim(data)[2]
  
  MSE<-c()
  
  folds <- createFolds(row.names(data), k=nfolds)
  if(length(set.alpha.v)==0){set.alpha.v=seq(from=0,to=0.9,0.1)}
  for(a in 1:length(set.alpha.v))
  {
    
    if(tau.v==1.4){taus<-runif(ntau,1,(1-set.alpha.v[a])*sqrt(ncol)+set.alpha.v[a])}
    else{taus<-rep(tau.v, ntau)}
    
    alpha.v<-set.alpha.v[a]
    CV.error<-c()
    
    for (k in 1:nfolds)
    {
      x.training<-data[folds[[k]],]
      x.test<-data[-folds[[k]],]
      CV.error.taus<-c()
      for(t in 1:length(taus))
      {
        res.csvd<-csvd.enet(as.matrix(x.training), Q=Q, 
                            tau.u = rep(sqrt(dim(x.training)[1]),Q), tau.v = rep(taus[t],Q),
                            alpha.u=alpha.u, alpha.v=alpha.v,
                            itermax.pi=itermax.pi, itermax.pocs=itermax.pocs,eps.pi=eps.pi, eps.pocs=eps.pocs, 
                            init.svd=init.svd)
        V.matrix<-res.csvd$V[,1:Q]
        Vt<-t(V.matrix)
        x.aprox.test<-x.test%*%V.matrix%*%Vt
        e<-norm((x.test-x.aprox.test),type="f")^2
        CV.error.taus<-c(CV.error.taus,e)
      }
      CV.error<-c(CV.error, mean(CV.error.taus))
    }
    MSE<-c(MSE,mean(CV.error))
    
  }
  
  min=min(MSE)+0.01
  df.u<-data.frame("Alpha"=set.alpha.v,"MSE"=MSE)
  
  r<-ggplot(df.u, aes(x=Alpha, y=MSE))+
    geom_point(data=subset(df.u,MSE<min), colour="red", size=3)+
    geom_point()+
    geom_line()
  
  return(list(MSE=MSE, MSE.plot=r))
}

#csvd.enet
csvd.enet<-function(X, Q=2, tau.u = rep(1.4, Q), tau.v = rep(1.4, Q),alpha.u=1e-16, alpha.v=1e-16,
         itermax.pi=1000, itermax.pocs=1000,
         eps.pi=1e-16, eps.pocs=1e-16, init.svd="svd") {
  
  if(length(tau.v)<=1){tau.v<-rep(tau.v,Q)}
  if(length(tau.u)<=1){tau.u<-rep(tau.u,Q)}
  X<-as.matrix(X)
  n.comp<-Q
  I <- nrow(X)
  J <- ncol(X)
  if (I==1 & J==1) stop("Are you sure you want to perform the SVD of a scalar?")
  nas <- is.na(X)
  if (sum(nas) > 0) {X<-na.omit(X) #Se omiten los valores perdidos
  print(sprintf("number of NAs:%i", sum(nas), "Se omiten las observaciones con valores perdidos"))}
  
  if (init.svd=="svd") {
    svdx <- svd(X, nu=n.comp, nv=n.comp)
    U0 <- svdx$u
    V0 <- svdx$v
  } else if( init=="rand") {
    U0 <- 1/(I-1) * mvrnorm(n = I, mu = rep(0,n.comp),
                            Sigma = diag(n.comp), empirical = TRUE)
    V0 <- 1/(J-1) * mvrnorm(n = J, mu = rep(0,n.comp),
                            Sigma = diag(n.comp), empirical = TRUE)
  } else {
    stop("init should be either svd or rand.")
  }
  
  U <- matrix(0, I, n.comp)
  V <- matrix(0, J, n.comp)
  iter <- rep(NA, n.comp)
  
  
  res.powit1 <- powit.enet1(X,
                            #U0[,1], V0[,1],
                            U0[,1,drop=FALSE], V0[,1,drop=FALSE],
                            tau.u[1], tau.v[1],alpha.u, alpha.v, 
                            eps.pi,
                            itermax.pi)
  U[,1] <- res.powit1$U
  V[,1] <- res.powit1$V
  iter[1] <- res.powit1$iter
  
  if (n.comp > 1) {
    for (r in 2:n.comp) {
      
      res.powit2 <- powit.enet2(X, # Matriz original
                                U0[,r], V0[,r], # vectores iniciales
                                U[,1:(r-1),drop=FALSE], V[,1:(r-1),drop=FALSE], # restriccion ortogonalidad. No se quita lo de DROP porque necesito que sean matrices
                                tau.u[r], tau.v[r], # restriccion sparse -> Elastic-net (Lasso=Tau, 1-alpha, Ridge=alpha)
                                alpha.u, alpha.v,
                                eps.pi, eps.pocs, #Precision
                                itermax.pi, itermax.pocs) # Iteraciones maximas
      U[,r] <- res.powit2$U
      V[,r] <- res.powit2$V
      iter[r] <- res.powit2$iter
    }
  }
  D <- diag(t(U) %*% X %*% V)
  oD <- order(D, decreasing = TRUE)
  res <- list(U=U[,oD], V=V[,oD], D=D[oD], iter=iter)
  return(res)
}

#Biplot.enet
Biplot.enet<-function(X, Q=2, tau.u = 1.4, tau.v = 1.4,alpha.u=1e-16, alpha.v=1e-16,biplot.type=2,
         itermax.pi=1000, itermax.pocs=1000,eps.pi=1e-16, eps.pocs=1e-16, 
         init.svd="svd", init.transf=1,  plot.axis=c(1,2),names.obs=FALSE, 
         select.cur=FALSE, variables.cur=1, weighted.cur=FALSE, method.cur="top.scores")
{
  
  library(ggcorrplot)
  library(ggplot2)
  library(reshape2)
  
  cur.res<-c("No results")
  plot.lev<-c("No results")
  
  if(select.cur)
  {
    library(rCUR)
    cur.res<-CUR(X,c=(dim(X)[2]-1),method=method.cur,weighted = weighted.cur)
    plot.lev<-plotLeverage(cur.res)
    #abline(v=sort(cur.res$leverage)[variables.cur]*1000,col="red")
    X<-X[,cur.res@C.index[1:variables.cur]]
  }
  n.comp<-Q
  varnames<-colnames(X)
  obsnames<-row.names(X)
  
  Corr.matrix=cor(X)
  row.names(Corr.matrix)<-colnames(X)
  colnames(Corr.matrix)<-colnames(X)
  
  cor.plot<-ggcorrplot(cor(as.matrix(X)), hc.order = FALSE, outline.col = "white",
                       ggtheme = ggplot2::theme_gray,
                       colors = c("#E46726", "white", "#6D9EC1"))
  
  #Scaling data:
  X<-init.transformation(X,init.transf)
  
  #Constrained SVD 
  res.csvd<-csvd.enet(X, Q=n.comp, tau.u = rep(tau.u,Q), tau.v = rep(tau.v,Q),alpha.u=alpha.u, alpha.v=alpha.v,
                      itermax.pi=itermax.pi, itermax.pocs=itermax.pocs,eps.pi=eps.pi, eps.pocs=eps.pocs, 
                      init.svd=init.svd)
  
  D.matrix<-res.csvd$D
  B.matrix<-(res.csvd$V)[,1:Q]%*%diag(D.matrix[1:Q])
  A.matrix<-(res.csvd$U)[,1:Q]%*%diag(D.matrix[1:Q])
  
  if(biplot.type==2)
  {
    B.matrix=B.matrix
    A.matrix=A.matrix
  }
  #Biplot type (biplot.type=0 -> GH-Biplot; biplot.type=1 -> JK-Biplot -> biplot.type en (1,inf] HJ-Biplot)
  if(biplot.type<=1)
  {
    A.matrix = A.matrix %*% diag((1/D.matrix[1:Q])^(1 - biplot.type))
    B.matrix = B.matrix %*% diag((1/D.matrix[1:Q])^biplot.type)
  }
  
  variance.proportion<-D.matrix[1:Q]^2/(norm(X,"f")^2)*100
  
  Amax<-A.matrix
  Borden<-B.matrix
  sumaA<-sum(Amax^2)
  sumaB<-sum(Borden^2)
  sA<-sumaA/(dim(Amax)[1])
  sB<-sumaB/(dim(Borden)[1])
  scf<-((sB/sA)^(1/2))^(1/2)
  Amax<-Amax*scf
  Bmax<-Borden/scf
  
  dimnames <- paste("PC",1:Q,sep="")
  row.names(Amax)<-obsnames
  colnames(Amax)<-dimnames
  row.names(Bmax)<-varnames
  colnames(Bmax)<-dimnames
  
  library(ggrepel)
  library(ggplot2)
  i=plot.axis[1]
  j=plot.axis[2]
  df.u<-as.data.frame(Amax)
  df.v<-as.data.frame(Bmax)
  row.names(df.u)<-row.names(Amax)
  colnames(df.u)<-colnames(Amax)
  row.names(df.v)<-row.names(Bmax)
  colnames(df.v)<-colnames(Bmax)
  df.v<-df.v[rowSums(df.v[,plot.axis])!=0,]
  
  axis.labels<-paste(colnames(df.u),sprintf("(%0.2f%%)", variance.proportion))
  
  p = ggplot()+
    geom_point(data = df.u, aes(x=df.u[,i], y=df.u[,j]),color="blue",size=2)+
    geom_segment(data = df.v, aes(x=0,y=0,xend=df.v[,i],yend=df.v[,j]),
                 arrow=arrow(length=unit(2/3,"picas"), type="closed"),size=1)+ 
    theme_bw()+
    geom_hline(yintercept=0, color="gray", size=0.5)+
    geom_vline(xintercept=0, color="gray", size=0.5)+
    geom_label_repel(data=df.v,aes(label=row.names(df.v),x=df.v[,i],y=df.v[,j],fontface=2),hjust=1,
                     color="black",size=2.5) +
    labs(x = axis.labels[i], y=axis.labels[j])
  
  if(names.obs)
  {
    p=p+geom_text(data=df.u, aes(label=row.names(df.u),x=df.u[,i],y=df.u[,j]),size=2.5)
  }
  
  
  return(list(Corr.matrix=Corr.matrix,
              Row.coordinates=Amax,
              Col.coordinates=Bmax,
              Pseudo.singular.values=res.csvd$D,
              Variance=variance.proportion,
              Cum.Variance=cumsum(variance.proportion), 
              corrplot=cor.plot,
              biplot.plot=p,
              results.cur=cur.res))
}

#bic.tau.enet
bic.tau.enet<-function(data, Q, ntau=3,tau.u=sqrt(dim(data)[1]), alpha.u=0.000001,alpha.v=0.5,
         itermax.pi=500, itermax.pocs=500,eps.pi=1e-16, eps.pocs=1e-16, 
         init.svd="svd")
{
  library(ggplot2)
  alpha.u<-alpha.u
  alpha.v<-alpha.v
  
  data<-as.matrix(data)
  nrow<-dim(data)[1]
  ncol<-dim(data)[2]
  
  BIC.error<-c()
  set.tau.v<-seq(from=1, to=((1-alpha.v)*sqrt(ncol)+alpha.v), length.out=ntau)
  
  for(a in 1:length(set.tau.v))
  {
    
    V.csvd<-csvd.enet(data, Q=Q, tau.u = rep(nrow,Q), tau.v = rep(set.tau.v[a],Q),
                      alpha.u=alpha.u, alpha.v=alpha.v,
                      itermax.pi=itermax.pi, itermax.pocs=itermax.pocs,
                      eps.pi=eps.pi, eps.pocs=eps.pocs, 
                      init.svd=init.svd)$V[,1:Q]
    
    Vt.csvd<-t(V.csvd)
    x.aprox.csvd<-data%*%V.csvd%*%Vt.csvd
    error.csvd<-norm((data-x.aprox.csvd),type="f")^2
    
    V.svd<-svd(data)$v[,1:Q]
    Vt.svd<-t(V.svd)
    x.aprox.svd<-data%*%V.svd%*%Vt.svd
    error.svd<-norm((data-x.aprox.svd),type="f")^2
    
    
    df<-length(which(V.svd!=0))
    bic<-(error.csvd/error.svd)-(df)*(log(nrow)/nrow)
    BIC.error<-c(BIC.error,bic)
    
  }
  
  df.u<-data.frame("Tau"=set.tau.v,"BIC"=BIC.error)
  
  b<-ggplot(df.u, aes(x=Tau, y=BIC))+
    geom_point()+
    geom_line()
  
  return(list(tau=set.tau.v, BIC=BIC.error, BIC.plot=b))
}

