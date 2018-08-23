function [loglikelihood]=proposal(optimizedParameters,RJMCMCType,modelProblem,shouldDisplay,algorithm,valArray,N,K,K1,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables,priorRange,modelUpdateWeight,modelWeightWithinModel,modelNumTReactions,moveNumber,currentInfluentialReactions,proposedInfluentialReactions,nReactions,AddDeleteReaction,reactionType,val,withinModelMoveNumber)
 
if algorithm==2
    
    if AddDeleteReaction==1
        
        % proposal for the direction in which the reaction is being added
        
        if RJMCMCType~=4
            % Evaluate the continuous component of the proposal
            [loglikelihood1]=ContinuousProposal(optimizedParameters,modelProblem,currentInfluentialReactions,proposedInfluentialReactions,nReactions,algorithm,valArray,priorRange,reactionType,AddDeleteReaction,K1,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables);
            
                                 
            loglikelihood=log(modelUpdateWeight(moveNumber))+loglikelihood1;

        else
                        
            loglikelihood=log(modelUpdateWeight(moveNumber));
            
        end
        
    elseif AddDeleteReaction==2
        
        % Within model random-walk
        
        % The continuous part of the proposal neednt be computed since the
        % proposal is symmetric
        loglikelihood=1.0;%log(modelWeightWithinModel(withinModelMoveNumber));
        
    end
        
end
    

