function [HessianMatrix]=hessianCalculator(x,updateReactions,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector)

% Diagonal elements of the hessian matrix by central differences

l=length(x);

HessianMatrix=zeros(l,l);

stepSize=nthroot(eps/2,4);

for i=1:l
    
    for j=1:i
        
        if i==j
            
            xTemp=x;
            
            xTemp(i)=x(i)+stepSize;
                        
            y1=optimizationObjective(xTemp,updateReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
            
            xTemp(i)=x(i);
            
            y2=optimizationObjective(xTemp,updateReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
            
            xTemp(i)=x(i)-stepSize;
            
            y3=optimizationObjective(xTemp,updateReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
            
            HessianMatrix(i,i)=(y1-2*y2+y3)/stepSize^2;               
        
        else    
            
           
            xTemp1=x;
            
            xTemp1(i)=xTemp1(i)+stepSize;
            xTemp1(j)=xTemp1(j)+stepSize;                        
            
            y1=optimizationObjective(xTemp1,updateReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
            
            xTemp2=x;
            
            xTemp2(i)=xTemp2(i)+stepSize;
            xTemp2(j)=xTemp2(j)-stepSize; 
            
            y2=optimizationObjective(xTemp2,updateReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
                        
            xTemp3=x;
            
            xTemp3(i)=xTemp3(i)-stepSize;
            xTemp3(j)=xTemp3(j)+stepSize;           
            
            y3=optimizationObjective(xTemp3,updateReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
                        
            xTemp4=x;
                        
            xTemp4(i)=xTemp4(i)-stepSize;
            xTemp4(j)=xTemp4(j)-stepSize;            
            
            y4=optimizationObjective(xTemp4,updateReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
              
            HessianMatrix(i,j)=(y1-y2-y3+y4)/4/stepSize/stepSize;
            
            HessianMatrix(j,i)=HessianMatrix(i,j);
                        
        end
    end
end

% This function returns the Hessian of -logposterior or -loglikelihood,
% thus it is positive semi-definite at the optimal point
