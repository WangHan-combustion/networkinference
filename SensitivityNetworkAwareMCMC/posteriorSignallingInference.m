function [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posteriorSignallingInference(key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector)
                                                                
if modelProblem==7
        
    % Calculates the logarithm of the numerator of the Hastings ratio    
    logp1=0.0;
        
    % Parameter in allowed prior range
    if logp1~=-inf
       p1=0.1*(valArray(1)~=delta)*(valArray(2)~=delta)*normal2(valArray(1),valArray(2),[1;2],[2 0;0 2])+0.9*(valArray(1)~=delta)*(valArray(2)==delta)*normal(valArray(1),1,2);
      
       logp1=log(p1); 
       
       lp1=logp1;
       ll1=logp1;
    end    

elseif modelProblem==8    
    
     [ll1,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,input,data,valArray,pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);
     lp1=logprior(N,valArray,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateProposedInfluentialReactions,meanVector,valVector);
     logp1=ll1+lp1;    

elseif modelProblem==9
    

        
        lp1=logprior(N,valArray,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateProposedInfluentialReactions,meanVector,valVector);
        
        [ll1,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,input,data,valArray,pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);
            
        logp1=ll1+lp1;
    
elseif modelProblem==10     

     lp1=logprior(N,valArray,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateProposedInfluentialReactions,meanVector,valVector);
                               
     [ll1,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,input,data,valArray,pos,modelProblem,dataVariance,nReactions,nSpecies,shoudlDisplay);
     logp1=ll1+lp1;
     
     
elseif modelProblem==5
        
             
    logp1=0.0;
    
    % Check if the parameter is within the allowed prior range
    for n=1:N

               
        if valArray(n)~=delta && (valArray(n)<priorRange(1,n) || valArray(n)>priorRange(2,n))                         
            lp1=-inf;
            logp1=-inf;
            break;
        end
    end
    
    ll1=0.0;
    lp1=0.0;
    
    % Parameter in allowed prior range
    if logp1~=-inf
            
        valPosterior=valArray;        
                      
        numberOfVariables=0;
        reactionVector=[];
        x=[];
        normalMean=[];
        normalSigma=[];
        
        for n=1:N
            
            if valPosterior(n)~=delta
                
                numberOfVariables=numberOfVariables+1;
                reactionVector(numberOfVariables)=n;    % Here I am just using it as a reference for variables/not reactions
                x(numberOfVariables)=valPosterior(n);
                normalMean(numberOfVariables)=(priorRange(2,n)+priorRange(1,n))/2;
                normalSigma(numberOfVariables)=(priorRange(2,n)-priorRange(1,n))/5;
                
            end
            
            if valArray(n)~=delta && valPosterior(n)==delta
                
                logp1=logp1-log(priorRange(2,n)-priorRange(1,n));
                
            end
            
        end
                
        logp1=logp1+logTruncatedNormalPosterior(x,normalMean,normalSigma,N,reactionType,priorRange,reactionVector);
        
        
        if typeSampling==1
            
            for i=1:size(key,1)
                if key(i,:)==intermediateProposedInfluentialReactions'
                    index=i;
                    break;
                end
            end
            
            
            logp1=log(keyCount(index))+logp1;
            
        end
        
    end
    
    
elseif modelProblem==6    

        logp1=0.0;
    
    % Check if the parameter is within the allowed prior range
    for n=1:N

               
        if valArray(n)~=delta && (valArray(n)<priorRange(1,n) || valArray(n)>priorRange(2,n))                         
            lp1=-inf;
            logp1=-inf;
            break;
        end
    end
    
    ll1=0.0;
    lp1=0.0;
    
    % Parameter in allowed prior range
    if logp1~=-inf
            
        valPosterior=valArray;        
                     
        numberOfVariables=0;
        reactionVector=[];
        x=[];
        normalMean=[];
        normalSigma=[];
        
        for n=1:N
            
            if valPosterior(n)~=delta
                
                numberOfVariables=numberOfVariables+1;
                reactionVector(numberOfVariables)=n;    % Here I am just using it as a reference for variables/not reactions
                x(numberOfVariables)=valPosterior(n);
                normalMean(numberOfVariables)=(priorRange(2,n)+priorRange(1,n))/2;
                normalSigma(numberOfVariables)=(priorRange(2,n)-priorRange(1,n))/5;
                
            end
            
            if valArray(n)~=delta && valPosterior(n)==delta
                
                logp1=logp1-log(priorRange(2,n)-priorRange(1,n));
                
            end
            
        end
                

        logp1=logp1+logTruncatedNormalPosterior(x,normalMean,normalSigma,N,reactionType,priorRange,reactionVector);
        
        
        if typeSampling==1
            
            for i=1:size(key,1)
                if key(i,:)==intermediateProposedInfluentialReactions'
                    index=i;
                    break;
                end
            end
            
            
            logp1=log(keyCount(index))+logp1;
            
        end
        
    end
    
elseif modelProblem==11
            
        lp1=logprior(N,valArray,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateProposedInfluentialReactions,meanVector,valVector);
        
        [ll1,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,input,data,valArray,pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);
        logp1=ll1+lp1;    
    
elseif modelProblem==14
    

        
        lp1=logprior(N,valArray,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateProposedInfluentialReactions,meanVector,valVector);
        
        [ll1,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,input,data,valArray,pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);
        logp1=ll1+lp1;        
        
    
else
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calclulates the logarithm of the numerator of the Hastings ratio
    
    logp1=0.0;
    ll1=0.0;
    
    % Check if the parameter is within the allowed prior range
    for n=1:N
               
        if valArray(n)~=delta && (valArray(n)<priorRange(1,n) || valArray(n)>priorRange(2,n))                         
            lp1=-inf;
            logp1=-inf;
            break;
        end
    end
    
   
   
    % Parameter in allowed prior range
    if logp1~=-inf || typeOfCall==2

       [ll1,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,input,data,valArray,pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);       
        
        
        % Count the number of forward model solves        
        lp1=logprior(N,valArray,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateProposedInfluentialReactions,meanVector,valVector);
                
        logp1=ll1+lp1;
    end

end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ll2=leadLogLikelihood;
 lp2=logPos;
 logp2=logPos;    
    