function [AcceptanceRate,numForwardModelSolves,N,K,numSamples,formatSamples,priorRange,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,typeMoves] = RJMCMCSolver(Seed,N,K,numSamples,inputX,data,pos,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,dataVariance,modelWeightScalingRegular,modelWeightScalingSwap,delta,RJMCMCType,algorithm,typeSampling,shouldDisplay,OnTheFlyModelDetermination,currentRemovedReactions,AddDelMoveProb,WithinMoveProb,SwapMoveProb,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,optimizedParameters,derivativeBased,proposalType,meanVector,valVector)

numForwardModelSolves=0;

% Determine the initial set of influential and uninfluential reactions

if typeSampling==0

    % For this case the current InfluentialReactions should not be sorted
    
    currentInfluentialReactions=zeros(nReactions,1);
    
    for i=1:nReactions
        if ~any(i==currentRemovedReactions)
            currentInfluentialReactions(i)=i;
        end
    end
    
    currentKeyLocation=0;
    currentModelIndex=0;   
    
    intermediateCurrentInfluentialReactions=currentInfluentialReactions;

elseif typeSampling==1
    
    % For this case the current InfluentialReactions should not be sorted
    
    currentInfluentialReactions=zeros(nReactions,1);
    
    for i=1:nReactions
        if ~any(i==currentRemovedReactions)
            currentInfluentialReactions(i)=i;
        end
    end
    
    currentKeyLocation=0;
    currentModelIndex=0;        
    
    [intermediateCurrentInfluentialReactions,currentUninfluentialReactions] = identifyInfluentialReactions(initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,currentRemovedReactions);

    
    if ~any(intermediateCurrentInfluentialReactions==1)
        intermediateCurrentInfluentialReactions(1)=1;
        intermediateCurrentInfluentialReactions=sort(intermediateCurrentInfluentialReactions);
    end     
        
    if ~any(intermediateCurrentInfluentialReactions==2)
        intermediateCurrentInfluentialReactions(1)=2;
        intermediateCurrentInfluentialReactions=sort(intermediateCurrentInfluentialReactions);
    end    
 
    %{
    if ~any(intermediateCurrentInfluentialReactions==9)
        intermediateCurrentInfluentialReactions(1)=9;
        intermediateCurrentInfluentialReactions=sort(intermediateCurrentInfluentialReactions);
    end        
 
    if ~any(intermediateCurrentInfluentialReactions==4)
        intermediateCurrentInfluentialReactions(1)=4;
        intermediateCurrentInfluentialReactions=sort(intermediateCurrentInfluentialReactions);
    end     
    
    
    if ~any(intermediateCurrentInfluentialReactions==10)
        intermediateCurrentInfluentialReactions(1)=10;
        intermediateCurrentInfluentialReactions=sort(intermediateCurrentInfluentialReactions);
    end   
    
    if ~any(intermediateCurrentInfluentialReactions==12)
        intermediateCurrentInfluentialReactions(1)=12;
        intermediateCurrentInfluentialReactions=sort(intermediateCurrentInfluentialReactions);
    end     
    %}
    
    
end

% Set the prior probability
[priorRange]=priorSetup(pos,modelProblem,N,data,valVector);

% Initialize chain
[formatSamples,ll,logPos,logPosTrue,numForwardModelSolves,logPosterior] = initializeChain(typeSampling,inputX,data,modelProblem,N,numSamples,pos,priorRange,numForwardModelSolves,dataVariance,nReactions,reactionType,nSpecies,currentInfluentialReactions,intermediateCurrentInfluentialReactions,delta,shouldDisplay,key,keyCount,meanVector,valVector);

totalAcceptance = 0;
totalSamples    = 0;

typeMoves=zeros(20,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define all the sufficient statistics (to be used in the adaptive MCMC algorithm)

% Parameters used in the adaptive MCMC algorithm

adaptationBatchSize=1000;
burnSamples=10*adaptationBatchSize;
batchSamples=zeros(N,adaptationBatchSize);

% Default is adaptive parameters
flagAdap=1;

% Number of zero components
K1=1;

[~,~,wVector]=initializeProposalParameters(priorRange,delta,N,nReactions,reactionType,K,K1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the proposal probabilities
[currentPossibleMoves,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,reactionRatesDRGReactions,production,consumption,rDRGEPReactions,rDRGEPRemovedReactions] = initialProposalProbabilities(currentInfluentialReactions,intermediateCurrentInfluentialReactions,formatSamples(1:N,numSamples),priorRange,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,modelWeightScalingRegular,RJMCMCType,algorithm,typeSampling,AddDelMoveProb,WithinMoveProb,SwapMoveProb,inputX,key,keyAllModels,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,nonZeroPriorModel,wVector,K,K1,modelIndex,currentKeyLocation);

i=1;

currentReactionsSamples=ones(nReactions,numSamples);
proposedReactionsSamples=ones(nReactions,numSamples);


tic
totalTime=0;

interval=1000;  

while (i<=numSamples)      
   
  if mod(i,interval)==0
      
      i
      totalTime=totalTime+toc;
      toc
      tic
      
  end
     
  
  if i==1
    d1=formatSamples(1:N,numSamples);
  else
    d1=formatSamples(1:N,i-1); 
  end  
   
  if flagAdap==0
   [samples,acc1,num1,ll,logPos,logPosTrue,currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,currentPossibleMoves,currentModelIndex,currentKeyLocation,reactionRatesDRGReactions,production,consumption,AMoveNum,WMoveNum,SMoveNum,AMoveAcc,WMoveAcc,SMoveAcc,rDRGEPReactions,rDRGEPRemovedReactions,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,typeMoves,numKeys,numForwardModelSolves] = generateIndSampler(i,d1,N,K,K1,delta,ll,logPos,logPosTrue,pos,priorRange,inputX,data,dataVariance,currentInfluentialReactions,intermediateCurrentInfluentialReactions,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,currentPossibleMoves,currentModelIndex,currentKeyLocation,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,modelWeightScalingRegular,modelWeightScalingSwap,RJMCMCType,algorithm,typeSampling,OnTheFlyModelDetermination,reactionRatesDRGReactions,production,consumption,AMoveNum,WMoveNum,SMoveNum,AMoveAcc,WMoveAcc,SMoveAcc,AddDelMoveProb,WithinMoveProb,SwapMoveProb,rDRGEPReactions,rDRGEPRemovedReactions,shouldDisplay,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,wVector,typeMoves,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,numForwardModelSolves,optimizedParameters,derivativeBased,proposalType,meanVector,valVector);
  else
   [samples,acc1,num1,ll,logPos,logPosTrue,currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,currentPossibleMoves,currentModelIndex,currentKeyLocation,reactionRatesDRGReactions,production,consumption,AMoveNum,WMoveNum,SMoveNum,AMoveAcc,WMoveAcc,SMoveAcc,rDRGEPReactions,rDRGEPRemovedReactions,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,typeMoves,numKeys,numForwardModelSolves] = generateIndSampler(i,d1,N,K,K1,delta,ll,logPos,logPosTrue,pos,priorRange,inputX,data,dataVariance,currentInfluentialReactions,intermediateCurrentInfluentialReactions,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,currentPossibleMoves,currentModelIndex,currentKeyLocation,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,modelWeightScalingRegular,modelWeightScalingSwap,RJMCMCType,algorithm,typeSampling,OnTheFlyModelDetermination,reactionRatesDRGReactions,production,consumption,AMoveNum,WMoveNum,SMoveNum,AMoveAcc,WMoveAcc,SMoveAcc,AddDelMoveProb,WithinMoveProb,SwapMoveProb,rDRGEPReactions,rDRGEPRemovedReactions,shouldDisplay,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,wVector,typeMoves,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,numForwardModelSolves,optimizedParameters,derivativeBased,proposalType,meanVector,valVector);      
  end    
  
  currentReactionsSamples(1:nReactions,i)=currentInfluentialReactions;
  proposedReactionsSamples(1:nReactions,i)=proposedInfluentialReactions;  
    
  formatSamples(1:N,i)=samples;
  logPosterior(i)=logPosTrue;
     
  totalAcceptance = totalAcceptance+acc1;
  totalSamples = totalSamples+num1;
     
  i=i+1;
  
end

disp('Moves of different types')

disp('Number of add/delete move')
disp(AMoveNum)

disp('Number of add/delete acceptances')
disp(AMoveAcc)

disp('Number of within model move')
disp(WMoveNum)

disp('Number of within model acceptances')
disp(WMoveAcc)

disp('Number of swap reaction moves')
disp(SMoveNum)

disp('Number of swap reaction acceptances')
disp(SMoveAcc)

disp('acceptance ratio');
totalAcceptance
totalSamples
AcceptanceRate=totalAcceptance/totalSamples;

fname1 = sprintf('MP11Optimization%sSampling3000KTruncatedNormalTest%s%d.mat',typeOfSampling,attempt,Seed);
save(fname1);
