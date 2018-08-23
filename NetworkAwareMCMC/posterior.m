function [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posterior(RJMCMCType,key,keyCount,currentInfluentialReactions,proposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector)

% MAPK signalling pathway
if pos==8
    [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posteriorSignallingInference(RJMCMCType,key,keyCount,currentInfluentialReactions,proposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector);
end


