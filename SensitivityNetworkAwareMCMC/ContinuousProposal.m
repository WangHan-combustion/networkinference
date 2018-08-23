function [loglik]=ContinuousProposal(optimizedParameters,modelProblem,currentInfluentialReactions,proposedInfluentialReactions,nReactions,algorithm,X,priorRange,reactionType,AddDeleteReaction,K1,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables)

loglik=0;

variableNumber=0;

% Add/delete reaction (Includes both add/delete move and swap move)

% Create a vector of proposed values : Ysamples

if algorithm==2
    
    loglik=0;
    reactionVector=[];
    
    
    for n=1:nReactions
        
        if any(n==proposedInfluentialReactions) && any(n==currentInfluentialReactions)
            
        elseif any(n==proposedInfluentialReactions) %|| any(n==currentInfluentialReactions)
                            
                % Using the regular proposal
                reactionVector=[reactionVector n];
                                
                if reactionType(n)==1
                    
                    m1=centeringLocation(variableNumber+1);
                    m2=centeringLocation(variableNumber+2);
                                                           

                       
                       Ysamples(variableNumber+1)=X(n+sum(reactionType(1:n-1)));
                       Ysamples(variableNumber+2)=X(1+n+sum(reactionType(1:n-1)));
                                              
  
                    
                    variableNumber=variableNumber+2;
                    
                else
                    
                    m1=centeringLocation(variableNumber+1);
                    

                        
                       Ysamples(variableNumber+1)=X(n+sum(reactionType(1:n-1)));
                      
  
                    
                    variableNumber=variableNumber+1;
                    
                end
            
        end
        

    end

    
    if variableNumber~=0
     if proposalType==1
       loglik=logMultiVariateNormal(Ysamples,centeringLocation,eigenVectorMatrix,eigenValueMatrix,numberOfVariables);        
     elseif proposalType==2
       loglik=logMultiVariateBeta(X,centeringLocation,alphaParameters,betaParameters,gammaParameters,reactionVector,reactionType,nReactions,priorRange,optimizedParameters);            
     elseif proposalType==3    
       loglik=logMultiVariateTruncatedNormal(Ysamples,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,reactionVector,reactionType,nReactions,priorRange,optimizedParameters);  
     end        
    end 
    
elseif algorithm==3
    
    loglik=0;
    
    for n=1:nReactions
        
        if any(n==proposedInfluentialReactions) && any(n==currentInfluentialReactions)
            
        elseif any(n==proposedInfluentialReactions) || any(n==currentInfluentialReactions)
            
            if reactionType(n)==1
                loglik=loglik-1/2*log(2*22/7*varVector(K1+1,n+sum(reactionType(1:n-1))))-((X(n+sum(reactionType(1:n-1)))-mVector(K1+1,n+sum(reactionType(1:n-1))))^2)/2/varVector(K1+1,n+sum(reactionType(1:n-1)))-1/2*log(2*22/7*varVector(K1+1,1+n+sum(reactionType(1:n-1))))-((X(1+n+sum(reactionType(1:n-1)))-mVector(K1+1,1+n+sum(reactionType(1:n-1))))^2)/2/varVector(K1+1,1+n+sum(reactionType(1:n-1)));
            else
                loglik=loglik-1/2*log(2*22/7*varVector(K1+1,n+sum(reactionType(1:n-1))))-((X(n+sum(reactionType(1:n-1)))-mVector(K1+1,n+sum(reactionType(1:n-1))))^2)/2/varVector(K1+1,n+sum(reactionType(1:n-1)));
            end
            
        end
    end
end
       


