function [logN1,logD1,ll1,logp1,logp2,logProp1,logProp2,numForwardModelSolves,AdditionalLnPrior]=calcN1D1(RJMCMCType,key,keyCount,typeSampling,algorithm,N,K,K1,proposalType,centeringLocationForward,centeringLocationOnBoundaryForward,eigenVectorMatrixForward,eigenValueMatrixForward,gradientVectorForward,HessianMatrixForward,numberOfVariablesForward,centeringLocationBackward,centeringLocationOnBoundaryBackward,eigenVectorMatrixBackward,eigenValueMatrixBackward,gradientVectorBackward,HessianMatrixBackward,numberOfVariablesBackward,yParameters,dParameters,dParametersForward,dParametersBackward,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,currentModelUpdateWeight,reverseModelUpdateWeight,currentModelWeightWithinModel,reverseModelWeightWithinModel,currentNumTReactions,reverseNumTReactions,currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,intermediateCurrentInfluentialReactionsForward,intermediateProposedInfluentialReactionsForward,intermediateCurrentInfluentialReactionsBackward,intermediateProposedInfluentialReactionsBackward,currentMoveNumber,reverseMoveNumber,currentWithinModelMoveNumber,reverseWithinModelMoveNumber,AddDeleteReaction,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,AddDeleteReactionSwap,shouldDisplay,numForwardModelSolves,optimizedParameters,meanVector,valVector)

% typeOfCall 1 indicates call to posterior for sampling
typeOfCall=1;

if shouldDisplay==1 
      disp('posterior')
end 
              
if typeSampling==1
           
      %if ~all(intermediateCurrentInfluentialReactions==intermediateProposedInfluentialReactions) || AddDeleteReaction==2
           
          diffReactions=setdiff(proposedInfluentialReactions,intermediateProposedInfluentialReactions);
          
          [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posterior(key,keyCount,intermediateProposedInfluentialReactions,N,yParameters,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector);
               
          AdditionalLnPrior=logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions,meanVector,valVector);

      %{         
      else
               
          logp1=logPos;
          logp2=logPos;
          ll1=ll;
           
          
          diffReactions=setdiff(proposedInfluentialReactions,intermediateProposedInfluentialReactions);
          
          AdditionalLnPrior=logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions,meanVector,valVector);          
          
                          
      end
      %}     
else
           
      [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posterior(key,keyCount,intermediateProposedInfluentialReactions,N,yParameters,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector);
           
      AdditionalLnPrior=0;
           
end


        if shouldDisplay==1
            disp('proposal')
        end  



if AddDeleteReaction==2
 % within model move   
 [logProp1]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,dParameters,N,K,K1,proposalType,centeringLocationBackward,centeringLocationOnBoundaryBackward,eigenVectorMatrixBackward,eigenValueMatrixBackward,gradientVectorBackward,HessianMatrixBackward,numberOfVariablesBackward,priorRange,reverseModelUpdateWeight,reverseModelWeightWithinModel,0,reverseMoveNumber,intermediateCurrentInfluentialReactionsBackward,intermediateProposedInfluentialReactionsBackward,nReactions,AddDeleteReaction,reactionType,dParametersBackward,reverseWithinModelMoveNumber); 
elseif AddDeleteReaction==3

elseif AddDeleteReaction==1            
 % add/delete move  
 [logProp1]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,dParameters,N,K,K1,proposalType,centeringLocationBackward,centeringLocationOnBoundaryBackward,eigenVectorMatrixBackward,eigenValueMatrixBackward,gradientVectorBackward,HessianMatrixBackward,numberOfVariablesBackward,priorRange,reverseModelUpdateWeight,reverseModelWeightWithinModel,reverseNumTReactions(reverseMoveNumber),reverseMoveNumber,intermediateCurrentInfluentialReactionsBackward,intermediateProposedInfluentialReactionsBackward,nReactions,AddDeleteReaction,reactionType,dParametersBackward);  
end

logN1=logp1+logProp1;

if AddDeleteReaction==2
 % within model move   
 [logProp2]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,yParameters,N,K,K1,proposalType,centeringLocationForward,centeringLocationOnBoundaryForward,eigenVectorMatrixForward,eigenValueMatrixForward,gradientVectorForward,HessianMatrixForward,numberOfVariablesForward,priorRange,currentModelUpdateWeight,currentModelWeightWithinModel,0,currentMoveNumber,intermediateCurrentInfluentialReactionsForward,intermediateProposedInfluentialReactionsForward,nReactions,AddDeleteReaction,reactionType,dParametersForward,currentWithinModelMoveNumber);
  
elseif AddDeleteReaction==3

elseif AddDeleteReaction==1
 % add/delete move   
 [logProp2]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,yParameters,N,K,K1,proposalType,centeringLocationForward,centeringLocationOnBoundaryForward,eigenVectorMatrixForward,eigenValueMatrixForward,gradientVectorForward,HessianMatrixForward,numberOfVariablesForward,priorRange,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions(currentMoveNumber),currentMoveNumber,intermediateCurrentInfluentialReactionsForward,intermediateProposedInfluentialReactionsForward,nReactions,AddDeleteReaction,reactionType,dParametersForward); 
end    

logD1=logp2+logProp2;

