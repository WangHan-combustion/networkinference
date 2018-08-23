function [AcceptanceRate,numForwardModelSolves,N,K,numSamples,formatSamples,priorRange,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,typeMoves] = multipleRealizations(N,K,numSamples,inputX,data,pos,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,dataVariance,modelWeightScalingRegular,modelWeightScalingSwap,delta,RJMCMCType,algorithm,typeSampling,shouldDisplay,OnTheFlyModelDetermination,currentRemovedReactions,AddDelMoveProb,WithinMoveProb,SwapMoveProb,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,optimizedParameters,derivativeBased,proposalType,meanVector,valVector)

tic

for i=1
   
 % Timer set

 rng(i,'twister')
 rng(i,'twister')    

 currentRemovedReactions=sort([currentRemovedReactions; sort(randsample([3 5 6 13 11 8],unidrnd(6)))']);
 
 [AcceptanceRate,numForwardModelSolves,N,K,numSamples,formatSamples,priorRange,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,typeMoves] = RJMCMCSolver(i,N,K,numSamples,inputX,data,pos,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,dataVariance,modelWeightScalingRegular,modelWeightScalingSwap,delta,RJMCMCType,algorithm,typeSampling,shouldDisplay,OnTheFlyModelDetermination,currentRemovedReactions,AddDelMoveProb,WithinMoveProb,SwapMoveProb,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,optimizedParameters,derivativeBased,proposalType,meanVector,valVector);
     
end    
