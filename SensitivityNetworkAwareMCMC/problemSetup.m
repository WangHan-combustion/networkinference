function [Seed,algorithm,RJMCMCType,typeSampling,shouldDisplay,OnTheFlyModelDetermination,numSamples,K,modelProblem,initialIncludedSpeciesPool,targetNode,N,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,pos,inputX,data,nData,dataVariance,modelWeightScalingRegular,modelWeightScalingSwap,delta,currentRemovedReactions,AddDelMoveProb,WithinMoveProb,SwapMoveProb,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,optimizedParameters,derivativeBased,proposalType,meanVector,valVector] = syntheticDataGeneration()

Seed=1;

rng(Seed,'twister')
rng(Seed,'twister')
%%
% This script is used to generate synthetic data for all problems

    % Simulation settings
    pos=8;
    RJMCMCType=2;
    algorithm=2;
    K=2;
    
    AddDelMoveProb=0.8;
    WithinMoveProb=0.2;
    SwapMoveProb=0.0;
        
    modelProblem=11;
    typeSampling=1;
    
    if modelProblem==8 || modelProblem==9 || modelProblem==10 || modelProblem==11 || modelProblem==14
     proposalType=1;
    else 
     proposalType=3;  % 1. normal distribution 2. beta distribution 3. truncated normal distribution
    end
    
    optimizedParameters=1;
    derivativeBased=1;
    shouldDisplay=0;
    OnTheFlyModelDetermination=1;

    %% 13 reaction problem (model problem 11) unconstrained parameter
    if modelProblem==11
        
        numSamples=1000;
                                
        dataStandardDeviation=0.2;
        dataVariance=dataStandardDeviation*dataStandardDeviation;
        delta=-100;
                
        load timesExample2nData30Final.mat
        load dataExample2nData30Final.mat

        nData=size(data,2);
        
        modelWeightScalingRegular=[0.0    0
                                   0.1    1000
                                   0.2    100
                                   0.3    100
                                   0.4    10
                                   0.5    10
                                   0.6    10
                                   0.7    5
                                   0.8    5
                                   0.9    5
                                   1.0    5
                                   2.0    0];
        
        modelWeightScalingSwap=[0.0   0
                                   0.1   1
                                   0.2   1
                                   0.3   1
                                   0.4   1
                                   0.5   5
                                   0.6   6
                                   0.7   7
                                   0.8   8
                                   0.9   9
                                   1.0  10
                                   2.0   0];
        
        targetNode=8;
        currentRemovedReactions=[7];
        
        nReactions=13;
        nSpecies=15;
        
        initialIncludedSpeciesPool=[3 4 5 8 11 12 13 15];
        speciesReference=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
        
        reactionType=[0 1 0 0 0 0 0 0 0 0 0 0 0];
        
        N=nReactions+sum(reactionType(1:nReactions));
        
        stoichiometricMatrix= [-1  1  3  0  0  0  0  0  0  0  0  0  3
                                1  0  0  0  0  0  0  0  0  0  0  0  0
                                0 -1  0  0  0  0  0  0  0  0  0  0  0
                                0 -1  0  0  0  0  0  0  0  0  0  0  0
                                0  0 -1  1  0  0  0  0  0  0  0  0  0
                                0  0  1 -1  3  0  0  0  0  0  0  0  0
                                0  0  0  0  1  3  0  0 -1  0  0  0  0
                                0  0  0  0  0 -1  1 -1  0  0  0  0  0
                                0  0  0  0  0  1 -1  1  0  0  0  0  0
                                0  0  0  0  0  0  0  3  0 -1  1  0  0
                                0  0  0  0  0  0  0  0  3  3  0  0  0
                                0  0  0  0 -1  0  0  0  1  0  0  0  0
                                0  0  0  0  0  0  0  0  0  1 -1  0  0
                                0  0  0  0  0  0  0  0  0  0  3 -1  1
                                0  0  0  0  0  0  0  0  0  0  0  1 -1];
                               
                         % 1       2       2       3       4       5       6       7       8       9      10      11      12      13
    
        offSet=0.7;            
                
        meanVector=[0.0 1.5 0.0 0.5+offSet 2.0 2.0+offSet 0.4+offSet 3.0 0.5 1.0 0.0 0.5 4.0 2.5];
    
                      % 1     2       2       3       4       5       6        7         8       9       10       11      12      13                        
    
        valVector=[  0.1   0.1     0.1     0.1     0.1     0.1     0.1      0.1      0.1       0.01       0.01     0.1     0.01      0.1];

        
        
    end
    %%      
    

    %% Determine all possible moves for the chosen problem
                
    [key,keyCount,keyMoveDetermination,keyAllModels,modelIndex,nonZeroPriorModel,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys] = determineAllPossibleMoves(typeSampling,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix);
       
