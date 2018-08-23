function [logN1,logD1,ll1,logp1,logp2,logProp1,logProp2,numForwardModelSolves,AdditionalLnPrior]=calcN1D1(RJMCMCType,key,keyCount,typeSampling,algorithm,N,K,K1,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables,yParameters,dParameters,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,currentModelUpdateWeight,reverseModelUpdateWeight,currentModelWeightWithinModel,reverseModelWeightWithinModel,currentNumTReactions,reverseNumTReactions,currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,currentMoveNumber,reverseMoveNumber,currentWithinModelMoveNumber,reverseWithinModelMoveNumber,AddDeleteReaction,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,AddDeleteReactionSwap,shouldDisplay,numForwardModelSolves,optimizedParameters,meanVector,valVector)
  
% typeOfCall 1 indicates call to posterior for sampling
typeOfCall=1;

if shouldDisplay==1 
      disp('posterior')
end 
              


if typeSampling==1
           
      if ~all(intermediateCurrentInfluentialReactions==intermediateProposedInfluentialReactions) || AddDeleteReaction==2
           
          diffReactions=setdiff(proposedInfluentialReactions,intermediateProposedInfluentialReactions);
          
          [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posterior(RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,yParameters,dParameters,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector);
               
          AdditionalLnPrior=logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions,meanVector,valVector);

               
      else
               
          logp1=logPos;
          logp2=logPos;
          ll1=ll;
           
          
          diffReactions=setdiff(proposedInfluentialReactions,intermediateProposedInfluentialReactions);
          
          AdditionalLnPrior=logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions,meanVector,valVector);          
          
          
          %{
          for n=1:nReactions
             if currentInfluentialReactions(n)==n && proposedInfluentialReactions(n)~=n
               
                 AdditionalLnPrior=AdditionalLnPrior-logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,n);                
                 break;
             elseif currentInfluentialReactions(n)~=n && proposedInfluentialReactions(n)==n
                                
                 AdditionalLnPrior=AdditionalLnPrior+logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,n);
                 break;
             end    
          end    
          %}
                          
      end
           
else
           
      [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posterior(RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,yParameters,dParameters,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector);
           
      AdditionalLnPrior=0;
           
end


        if shouldDisplay==1
            disp('proposal')
        end  



if AddDeleteReaction==2
 % within model move   
 [logProp1]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,dParameters,N,K,K1,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables,priorRange,reverseModelUpdateWeight,reverseModelWeightWithinModel,0,reverseMoveNumber,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,nReactions,AddDeleteReaction,reactionType,yParameters,reverseWithinModelMoveNumber); 
elseif AddDeleteReaction==3

elseif AddDeleteReaction==1            
 % add/delete move  
 [logProp1]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,dParameters,N,K,K1,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables,priorRange,reverseModelUpdateWeight,reverseModelWeightWithinModel,reverseNumTReactions(reverseMoveNumber),reverseMoveNumber,intermediateProposedInfluentialReactions,intermediateCurrentInfluentialReactions,nReactions,AddDeleteReaction,reactionType,yParameters);  
end

logN1=logp1+logProp1;

if AddDeleteReaction==2
 % within model move   
 [logProp2]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,yParameters,N,K,K1,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables,priorRange,currentModelUpdateWeight,currentModelWeightWithinModel,0,currentMoveNumber,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,nReactions,AddDeleteReaction,reactionType,dParameters,currentWithinModelMoveNumber);
  
elseif AddDeleteReaction==3

elseif AddDeleteReaction==1
 % add/delete move   
 [logProp2]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,yParameters,N,K,K1,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables,priorRange,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions(currentMoveNumber),currentMoveNumber,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,nReactions,AddDeleteReaction,reactionType,dParameters); 
end    

logD1=logp2+logProp2;

