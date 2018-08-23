function [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posterior(key,keyCount,proposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector)



% MAPK signalling pathway
if pos==8
    [logp1,logp2,ll1,ll2,lp1,lp2,numForwardModelSolves]=posteriorSignallingInference(key,keyCount,proposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,typeOfCall,shouldDisplay,meanVector,valVector);       
    
end

%{


% Catalysis example/Stagnation flow reactor
if pos==7
   [logp1,logp2,ll1,ll2]=posteriorSteamReforming(N,valArray,val,leadLogLikelihood,pos,priorRange,input,output,delta,dataVariance);
end

% 2-D example simulated
if pos==1
   [logp1,logp2,ll1,ll2]=posterior2Dexample(N,valArray,val);
end

% 5-D example simulated
if pos==2
    [logp1,logp2,ll1,ll2]=posterior5Dsimulated(N,valArray,val);
end

% 5-D example easy posterior
if pos==3
    [logp1,logp2,ll1,ll2]=posterior5DsimulatedSecond(N,valArray,val);
end

% 8-D linear regression
if pos==4
    [logp1,logp2,ll1,ll2]=posterior8DlinearRegression(N,valArray,val,leadLogLikelihood,pos,priorRange,input,output,delta,dataVariance);
end

% 10-D logistic regression
if pos==5
    [logp1,logp2,ll1,ll2]=posterior10DlogisticRegression(N,valArray,val,leadLogLikelihood,pos,priorRange,input,output,delta,dataVariance);
end

if pos==6
    [logp1,logp2,ll1,ll2]=posteriorNeuralNetwork(N,valArray,val,leadLogLikelihood,pos,priorRange,input,output,delta,dataVariance);
end

%}


