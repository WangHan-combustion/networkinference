function [f]=optimizationObjective(x,updateReactions,RJMCMCType,key,keyCount,currentInfluentialReactions,proposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector)

variableNumber=0;

for i=1:length(updateReactions)
        
    if reactionType(updateReactions(i))==1
        valArray(updateReactions(i)+sum(reactionType(1:(updateReactions(i)-1))))=x(variableNumber+1);
        valArray(1+updateReactions(i)+sum(reactionType(1:(updateReactions(i)-1))))=x(variableNumber+2);
        
        variableNumber=variableNumber+2;
    else
        valArray(updateReactions(i)+sum(reactionType(1:(updateReactions(i)-1))))=x(variableNumber+1);
                
        variableNumber=variableNumber+1;
        
    end    
    
end

% typeOfCall =2 call to posterior for proposal optimization
typeOfCall=2;
[logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posterior(RJMCMCType,key,keyCount,currentInfluentialReactions,proposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector);

if proposalType==1
 f=-logp1;
else
 f=-ll1;
end