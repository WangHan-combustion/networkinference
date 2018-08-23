function [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posteriorSignallingInference(RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector)

if modelProblem==11
            
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