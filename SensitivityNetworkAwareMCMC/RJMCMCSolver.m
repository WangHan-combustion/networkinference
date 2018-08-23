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
    
    
    %{
    if ~any(intermediateCurrentInfluentialReactions==11)
        intermediateCurrentInfluentialReactions(1)=11;
        intermediateCurrentInfluentialReactions=sort(intermediateCurrentInfluentialReactions);
    end     
    %}
    %{
    if ~any(intermediateCurrentInfluentialReactions==8)
        intermediateCurrentInfluentialReactions(1)=8;
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

S_O=zeros(K,N);
S_R=zeros(K,N);
S_M=zeros(K,N);
S_Q=zeros(K,N);

adaptationBatchSize=1000;
burnSamples=10*adaptationBatchSize;
batchSamples=zeros(N,adaptationBatchSize);

% Default is adaptive parameters
flagAdap=1;

% Number of zero components
K1=1;

[mVector,varVector,wVector]=initializeProposalParameters(priorRange,delta,N,nReactions,reactionType,K,K1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the proposal probabilities
[currentPossibleMoves,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,reactionRatesDRGReactions,production,consumption,rDRGEPReactions,rDRGEPRemovedReactions] = initialProposalProbabilities(currentInfluentialReactions,intermediateCurrentInfluentialReactions,formatSamples(1:N,numSamples),priorRange,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,modelWeightScalingRegular,RJMCMCType,algorithm,typeSampling,AddDelMoveProb,WithinMoveProb,SwapMoveProb,inputX,key,keyAllModels,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,nonZeroPriorModel,wVector,K,K1,modelIndex,currentKeyLocation);

i=1;

AMoveNum=0;
WMoveNum=0;
SMoveNum=0;

AMoveAcc=0;
WMoveAcc=0;
SMoveAcc=0;

% Prepare file for writing

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
  
  formatSamples(1:N,i)=samples;
  logPosterior(i)=logPosTrue;
     
  currentReactionsSamples(1:nReactions,i)=currentInfluentialReactions;
  proposedReactionsSamples(1:nReactions,i)=proposedInfluentialReactions;
  
  
  totalAcceptance = totalAcceptance+acc1;
  totalSamples = totalSamples+num1;
   
  typeOfSampling='MC';
  
  attempt='Test309d';

  if i==1000
  
   fname1 = sprintf('MP11Optimization%sSampling1KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end     

  if i==2000
  
   fname1 = sprintf('MP11Optimization%sSampling2KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end 
  
  if i==3000
  
   fname1 = sprintf('MP11Optimization%sSampling3KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  
  if i==5000
  
   fname1 = sprintf('MP11Optimization%sSampling5KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
    
  if i==10000
  
   fname1 = sprintf('MP11Optimization%sSampling10KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==20000
  
   fname1 = sprintf('MP11Optimization%sSampling20KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==30000
  
   fname1 = sprintf('MP11Optimization%sSampling30KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end 
  
  
  if i==50000
  
   fname1 = sprintf('MP11Optimization%sSampling50KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==60000
  
   fname1 = sprintf('MP11Optimization%sSampling60KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end     

  if i==70000
  
   fname1 = sprintf('MP11Optimization%sSampling70KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==75000
  
   fname1 = sprintf('MP11Optimization%sSampling75KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    
  
  if i==100000
  
   fname1 = sprintf('MP11Optimization%sSampling100KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==125000
  
   fname1 = sprintf('MP11Optimization%sSampling125KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  
  
  
  if i==150000
  
   fname1 = sprintf('MP11Optimization%sSampling150KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==175000
  
   fname1 = sprintf('MP11Optimization%sSampling175KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    
    
  if i==200000
  
   fname1 = sprintf('MP11Optimization%sSampling200KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==225000
  
   fname1 = sprintf('MP11Optimization%sSampling225KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==250000
  
   fname1 = sprintf('MP11Optimization%sSampling250KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  
 
  if i==275000
  
   fname1 = sprintf('MP11Optimization%sSampling275KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    
  
  if i==300000
  
   fname1 = sprintf('MP11Optimization%sSampling300KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    

  if i==325000
  
   fname1 = sprintf('MP11Optimization%sSampling325KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end      
  
  if i==350000
  
   fname1 = sprintf('MP11Optimization%sSampling350KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==375000
  
   fname1 = sprintf('MP11Optimization%sSampling375KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end     
  
  if i==400000
  
   fname1 = sprintf('MP11Optimization%sSampling400KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end     

  if i==425000
  
   fname1 = sprintf('MP11Optimization%sSampling425KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==450000
  
   fname1 = sprintf('MP11Optimization%sSampling450KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==475000
  
   fname1 = sprintf('MP11Optimization%sSampling475KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==500000
  
   fname1 = sprintf('MP11Optimization%sSampling500KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    

  if i==525000
  
   fname1 = sprintf('MP11Optimization%sSampling525KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==550000
  
   fname1 = sprintf('MP11Optimization%sSampling550KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==575000
  
   fname1 = sprintf('MP11Optimization%sSampling575KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==600000
  
   fname1 = sprintf('MP11Optimization%sSampling600KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    

  if i==625000
  
   fname1 = sprintf('MP11Optimization%sSampling625KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==650000
  
   fname1 = sprintf('MP11Optimization%sSampling650KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==675000
  
   fname1 = sprintf('MP11Optimization%sSampling675KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==700000
  
   fname1 = sprintf('MP11Optimization%sSampling700KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==725000
  
   fname1 = sprintf('MP11Optimization%sSampling725KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==750000
  
   fname1 = sprintf('MP11Optimization%sSampling750KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==775000
  
   fname1 = sprintf('MP11Optimization%sSampling775KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==800000
  
   fname1 = sprintf('MP11Optimization%sSampling800KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==825000
  
   fname1 = sprintf('MP11Optimization%sSampling825KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==850000
  
   fname1 = sprintf('MP11Optimization%sSampling850KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==875000
  
   fname1 = sprintf('MP11Optimization%sSampling875KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==900000
  
   fname1 = sprintf('MP11Optimization%sSampling900KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==925000
  
   fname1 = sprintf('MP11Optimization%sSampling925KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==950000
  
   fname1 = sprintf('MP11Optimization%sSampling950KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==975000
  
   fname1 = sprintf('MP11Optimization%sSampling975KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1000000
  
   fname1 = sprintf('MP11Optimization%sSampling1000KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  
  
  if i==1100000
  
   fname1 = sprintf('MP11Optimization%sSampling1100KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1125000
  
   fname1 = sprintf('MP11Optimization%sSampling1125KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  
  
  
  if i==1150000
  
   fname1 = sprintf('MP11Optimization%sSampling1150KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1175000
  
   fname1 = sprintf('MP11Optimization%sSampling1175KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    
    
  if i==1200000
  
   fname1 = sprintf('MP11Optimization%sSampling1200KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1225000
  
   fname1 = sprintf('MP11Optimization%sSampling1225KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==1250000
  
   fname1 = sprintf('MP11Optimization%sSampling1250KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  
 
  if i==1275000
  
   fname1 = sprintf('MP11Optimization%sSampling1275KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    
  
  if i==1300000
  
   fname1 = sprintf('MP11Optimization%sSampling1300KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    

  if i==1325000
  
   fname1 = sprintf('MP11Optimization%sSampling1325KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end      
  
  if i==1350000
  
   fname1 = sprintf('MP11Optimization%sSampling1350KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1375000
  
   fname1 = sprintf('MP11Optimization%sSampling1375KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end     
  
  if i==1400000
  
   fname1 = sprintf('MP11Optimization%sSampling1400KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end     

  if i==1425000
  
   fname1 = sprintf('MP11Optimization%sSampling1425KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==1450000
  
   fname1 = sprintf('MP11Optimization%sSampling1450KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1475000
  
   fname1 = sprintf('MP11Optimization%sSampling1475KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1500000
  
   fname1 = sprintf('MP11Optimization%sSampling1500KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    

  if i==1525000
  
   fname1 = sprintf('MP11Optimization%sSampling1525KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==1550000
  
   fname1 = sprintf('MP11Optimization%sSampling1550KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1575000
  
   fname1 = sprintf('MP11Optimization%sSampling1575KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1600000
  
   fname1 = sprintf('MP11Optimization%sSampling1600KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end    

  if i==1625000
  
   fname1 = sprintf('MP11Optimization%sSampling1625KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==1650000
  
   fname1 = sprintf('MP11Optimization%sSampling1650KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1675000
  
   fname1 = sprintf('MP11Optimization%sSampling1675KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1700000
  
   fname1 = sprintf('MP11Optimization%sSampling1700KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1725000
  
   fname1 = sprintf('MP11Optimization%sSampling1725KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==1750000
  
   fname1 = sprintf('MP11Optimization%sSampling1750KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1775000
  
   fname1 = sprintf('MP11Optimization%sSampling1775KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1800000
  
   fname1 = sprintf('MP11Optimization%sSampling1800KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1825000
  
   fname1 = sprintf('MP11Optimization%sSampling1825KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==1850000
  
   fname1 = sprintf('MP11Optimization%sSampling1850KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1875000
  
   fname1 = sprintf('MP11Optimization%sSampling1875KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1900000
  
   fname1 = sprintf('MP11Optimization%sSampling1900KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end  

  if i==1925000
  
   fname1 = sprintf('MP11Optimization%sSampling1925KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   
  
  if i==1950000
  
   fname1 = sprintf('MP11Optimization%sSampling1950KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==1975000
  
   fname1 = sprintf('MP11Optimization%sSampling1975KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end   

  if i==2000000
  
   fname1 = sprintf('MP11Optimization%sSampling2000KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
   save(fname1);  
  
  end 
  
  
  
  i=i+1;
  
end


%formatSamples=samples;

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

fname1 = sprintf('MP11Optimization%sSampling2000KTruncatedNormal%s%d.mat',typeOfSampling,attempt,Seed);
save(fname1);
