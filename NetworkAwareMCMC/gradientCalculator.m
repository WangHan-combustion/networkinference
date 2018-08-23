function [gradientVector]=gradientCalculator(x,updateReactions,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector)

% Forward difference gradient

l=length(x);

gradientVector=zeros(l,1);

stepSize=sqrt(eps/2);

for i=1:l

 xTemp=x;   
    
 xTemp(i)=x(i);
    
 y1=optimizationObjective(xTemp,updateReactions,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);

 xTemp(i)=x(i)+stepSize;
 
 y2=optimizationObjective(xTemp,updateReactions,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);

 gradientVector(i)=(y2-y1)/stepSize;
 
end


% Objective value at boundary


% This function returns the gradient of the loglikelihood or the
% logposterior 

% G to be negated since -logposterior is the objective
%gradientVector=-gradientVector;
