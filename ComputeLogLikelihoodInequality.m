function logL = ComputeLogLikelihoodInequality(times, a, p, b, c, mainshockMags, Mmin, G, t1, t2, a0, p0)

  % Basic log likelihood with time-varying Mc
  logL = sum(log(Lambda(times, a, p, b, c, mainshockMags, Mmin, G))) - integral(@(t) Lambda(t, a, p, b, c, mainshockMags, Mmin, G),t1,t2);

  % Now to add inequality constraint penalty...
  % (Comment this section out to remove inequality constraint)
  penalty = findPenalty(a, p, b, c, mainshockMags, Mmin, t1, t2, a0,  p0);
  logL = logL + penalty;

end
  
  
function lambda = Lambda(t, a, p, b, c, Mmain, Mmin, G)
% Omori rate of aftershocks above time-dependent Mc
  
  Mmain2=Mmain'*ones(1,length(t));
  t2=ones(length(Mmain),1)*t;
  
  Mc = max(Mmain2/2 - G - log10(t2), Mmin);
  
  lambda = 10.^(a+b*(Mmain2-Mc));
  if size(lambda,1)>1 % If there are multiple mainshocks in stack
    lambda=sum(lambda);
  end
  lambda = lambda .* (t+c).^-p;
  
end


function penalty = findPenalty(a, p, b, c, mainshockMags, Mmin, t1, t2, a0, p0)  
% Penalty for inequality constraint
% Keep time-dependent Mc fit above Mcat fit at all points in time interval
  
     % Equivalent mainshock magnitude for stacked sequence
     Mequiv=(1/b)*log10(sum(10.^(b*mainshockMags)));
 
     % Rate vs. Time for all data (down to catalog Mc)
     Lambda0=@(t) Lambda(t, a0, p0, b, c, Mequiv, Mmin, Inf);
     
     % Rate vs. Time for time-varying Mc
     Lambda1=@(t) Lambda(t, a, p, b, c, Mequiv, Mmin, Inf);
   
     tMaxDiff = fminbnd(@(t) Lambda1(t) - Lambda0(t), t1, t2);
     maxDiff = Lambda1(tMaxDiff) - Lambda0(tMaxDiff);     
     % fminbnd() can get stuck in local minima - so also check endpoints
     d1=Lambda1(t1) - Lambda0(t1);
     d2=Lambda1(t2) - Lambda0(t2);
     maxDiff = min(maxDiff,d1);
     maxDiff = min(maxDiff,d2);
     
     if maxDiff < 0
       % Only penalize if Lambda0 is greater than Lambda1
       penalty = 10^6*maxDiff; 
     else
       penalty = 0;
     end
     
end
   
