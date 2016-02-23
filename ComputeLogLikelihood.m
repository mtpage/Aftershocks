function logL = ComputeLogLikelihood(times, a, p, b, c, Mmain, Mmin, t1, t2)

  N=length(times);
  k=10^(a+b*(Mmain-Mmin));
  
  if p==1
    logL=N*log(k)-p*sum(log(times+c))-k*(log(t2+c)-log(t1+c));
  else
    logL=N*log(k)-p*sum(log(times+c))-k*((t2+c)^(1-p)-(t1+c)^(1-p))/(1-p);
  end
  
end
  
  
