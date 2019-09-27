'
int D=as<int>(D_);
int K=as<int>(K_);
int W=as<int>(W_);
IntegerVector  N_d(Nd_);
IntegerMatrix  subj_species_topic_matrix(data);
IntegerMatrix  subj_topic_matrix(D,K); 
IntegerMatrix  topic_species_matrix(K,W); 
IntegerVector  topic_vector(K);
NumericVector alpha(alpha_); 
NumericVector beta(beta_);
double sum_beta;
double sum_alpha;
NumericMatrix pi(D,K);
NumericMatrix phi(K,W);
int S=as<int>(S_);

for (int d=0; d<D; ++d){
  for (int k=0; k<K; ++k){
    subj_topic_matrix(d,k)=0; 
  }
}
for (int k=0; k<K; ++k){
  for (int w=0; w<W; ++w){
    topic_species_matrix(k,w)=0; 
  }
}
for (int k=0; k<K; ++k){
  topic_vector[k]=0;
}

NumericMatrix dataK(K,D);
NumericMatrix dataK2(K,D);
NumericVector dataK2_plus(K);
NumericVector dataK_plus(K);

NumericMatrix dataW(W,K);
NumericMatrix dataW2(W,K);
NumericVector dataW2_plus(W);
NumericVector dataW_plus(W);

NumericMatrix  pi_plus(D,K);
NumericMatrix  phi_plus(K,W);

int L=subj_species_topic_matrix.nrow();

for (int s=0; s<S; ++s){
  sum_alpha=0;
  sum_beta=0;
  for (int w=0; w<W; ++w){
    sum_beta+=beta[w];
  }
  for (int k=0; k<K; ++k){
    sum_alpha+=alpha[k];
  }
  for (int n=0; n<L; ++n){
    int h=subj_species_topic_matrix(n,0);
    int l=subj_species_topic_matrix(n,1);
    int assign=subj_species_topic_matrix(n,2);
    if (assign>0){
      subj_topic_matrix((h-1),(assign-1))-=1;
      topic_species_matrix((assign-1),(l-1))-=1;
      topic_vector[assign-1]-=1;
    }
    
    double prob(K); 
    NumericVector p(K);
    for (int k=0; k<K; ++k){
      prob=subj_topic_matrix((h-1),k)+alpha[k];
      prob=prob*(topic_species_matrix(k,(l-1))+beta[l-1]);
      prob=prob/(topic_vector[k]+sum_beta);
      p[k]=prob;
      if (k!=0){ p[k]+=p[k-1];}
    }
    double u=R::runif(0,1);
    u=p[K-1]*u;
    for (int k=0; k<K; ++k){
      assign= k;
      if (u<p[k]){ break ;}
    }
    
    subj_topic_matrix((h-1),(assign))+=1;
    topic_species_matrix((assign),(l-1))+=1; 
    topic_vector[assign]+=1;
    subj_species_topic_matrix(n,2)=(assign+1);
  }
  
  for ( int k=0; k<K; ++k){
    dataK_plus[k]=0;
    dataK2_plus[k]=0;
    for ( int d=0; d<D; ++d){
      dataK(k,d)=(R::digamma(subj_topic_matrix(d,k)+alpha[k])-R::digamma(alpha[k]));
      dataK2(k,d)=(R::digamma(N_d[d]+sum_alpha)-R::digamma(sum_alpha));
      dataK_plus[k]+=dataK(k,d);
      dataK2_plus[k]+=dataK2(k,d);
    }
    alpha[k]=alpha[k]*dataK_plus[k]/dataK2_plus[k];
  }
  
  for ( int w=0; w<W; ++w){
    dataW_plus[w]=0;
    dataW2_plus[w]=0;
    for ( int k=0; k<K; ++k){
      dataW(w,k)=(R::digamma(topic_species_matrix(k,w)+beta[w])-R::digamma(beta[w]));
      dataW2(w,k)=(R::digamma(topic_vector[k]+sum_beta)-R::digamma(sum_beta));
      dataW_plus[w]+=dataW(w,k);
      dataW2_plus[w]+=dataW2(w,k);
    }
    beta[w]=beta[w]*dataW_plus[w]/dataW2_plus[w];
  }
  
  for ( int k=0; k<K; ++k){
    for ( int w=0; w<W; ++w){
      phi(k,w)=(topic_species_matrix(k,w)+beta[w])/(topic_vector[k]+sum_beta);
      if (s>=1000){ phi_plus(k,w)=phi_plus(k,w)+phi(k,w);}
    }
  }
  
  for ( int d=0; d<D; ++d){
    for ( int k=0; k<K; ++k){
      pi(d,k)=(subj_topic_matrix(d,k)+alpha[k])/(N_d[d]+sum_alpha);
      if (s>=1000){  pi_plus(d,k)=pi_plus(d,k)+pi(d,k);}
    }
  }
  
}

List ret;
ret["phi_plus"]=phi_plus; ret["pi_plus"]=pi_plus;
return wrap(ret);
'