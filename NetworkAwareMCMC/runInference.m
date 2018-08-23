function [AcceptanceRate,numForwardModelSolves,N,K,numSamples,formatSamples,priorRange,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,typeMoves] = runInference(N,K,numSamples,inputX,data,pos,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,dataVariance,modelWeightScalingRegular,modelWeightScalingSwap,delta,RJMCMCType,algorithm,typeSampling,shouldDisplay,OnTheFlyModelDetermination,currentRemovedReactions,AddDelMoveProb,WithinMoveProb,SwapMoveProb,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,optimizedParameters,derivativeBased,proposalType,meanVector,valVector)

tic

for i=1
    
 % Timer set

 rng(i,'twister')
 rng(i,'twister')    
        
 currentRemovedReactions=sort([currentRemovedReactions; sort(randsample([3 4 5 6 9],unidrnd(5)))']);

 [AcceptanceRate,numForwardModelSolves,N,K,numSamples,formatSamples,priorRange,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,typeMoves] = RJMCMCSolver(i,N,K,numSamples,inputX,data,pos,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,dataVariance,modelWeightScalingRegular,modelWeightScalingSwap,delta,RJMCMCType,algorithm,typeSampling,shouldDisplay,OnTheFlyModelDetermination,currentRemovedReactions,AddDelMoveProb,WithinMoveProb,SwapMoveProb,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,optimizedParameters,derivativeBased,proposalType,meanVector,valVector);
     
end    
