function [currentPossibleMoves,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,reactionRatesDRGReactions,production,consumption,rDRGEPReactions,rDRGEPRemovedReactions] = initialProposalProbabilities(currentInfluentialReactions,intermediateCurrentInfluentialReactions,valVector,priorRange,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,modelWeightScalingRegular,RJMCMCType,algorithm,typeSampling,AddDelMoveProb,WithinMoveProb,SwapMoveProb,inputX,key,keyAllModels,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,nonZeroPriorModel,wVector,K,K1,modelIndex,currentKeyLocation)

% This function returns the weights for different move probabilities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the initial proposal move proabilities

% Making it +1 for within model move +1 for swap move
scaleLength=size(modelWeightScalingRegular,1);

reactionRatesDRGReactions=0;
production=0;
consumption=0;
rDRGEPReactions=0;
rDRGEPRemovedReactions=zeros(nReactions,nReactions);

% default values
currentPossibleMoves=0;

locationOfCurrent=0;

%% Chemical kinetic reversible-jump MCMC with model check (RJMCMCType 2)
if RJMCMCType==2 || RJMCMCType==4
    
    
    if typeSampling==1
 
        currentTotalPossibleMoves=nReactions+2;
        
        currentNumTReactions=zeros(currentTotalPossibleMoves,1);
        
        if algorithm==2
            
            currentModelUpdateWeight=zeros(currentTotalPossibleMoves,1);
            currentModelWeightWithinModel=zeros(nReactions,1);
            
            for i=1:nReactions
                
               if i~=7 && i~=2 && i~=1 && i~=10 && i~=12 && i~=1 && i~=8 && i~=11 && i~=13 && i~=1                currentNumTReactions(i)=1;
                currentModelUpdateWeight(i)=1;
               end 
                
            end
                        
        end
        
        % Within model move
        if nnz(intermediateCurrentInfluentialReactions)~=0
            currentModelUpdateWeight(currentTotalPossibleMoves-1)=1;
        else
            currentModelUpdateWeight(currentTotalPossibleMoves-1)=0;
        end
        
        % Swap move
        currentModelUpdateWeight(currentTotalPossibleMoves)=0;
        
        for i=1:nReactions
            if any(i==intermediateCurrentInfluentialReactions)
                currentModelWeightWithinModel(i)=1;
            else
                currentModelWeightWithinModel(i)=0;
            end
        end                
        
    elseif typeSampling==0
        
        currentTotalPossibleMoves=nReactions+2;
        
        currentNumTReactions=zeros(currentTotalPossibleMoves,1);
        
        if algorithm==2
            
            currentModelUpdateWeight=zeros(currentTotalPossibleMoves,1);
            currentModelWeightWithinModel=zeros(nReactions,1);
            
            for i=1:nReactions
                
               if i~=7 && i~=2 && i~=1 && i~=10 && i~=12 && i~=8 && i~=11 && i~=13 && i~=1%&& i~=2 && i~=12 && i~=10 && i~=7 && i~=6
                currentNumTReactions(i)=1;
                currentModelUpdateWeight(i)=1;
               end 
                
            end
                        
        end
        
        % Within model move
        if nnz(currentInfluentialReactions)~=0
            currentModelUpdateWeight(currentTotalPossibleMoves-1)=1;
        else
            currentModelUpdateWeight(currentTotalPossibleMoves-1)=0;
        end
        
        % Swap move
        currentModelUpdateWeight(currentTotalPossibleMoves)=0;
        
        for i=1:nReactions
            if any(i==currentInfluentialReactions)
                currentModelWeightWithinModel(i)=1;
            else
                currentModelWeightWithinModel(i)=0;
            end
        end
        
    elseif typeSampling==3
                
        index=currentKeyLocation;
        
        locationOfCurrent=index;
        
        currentPossibleMoves=AllPossibleMoves{index};
        numAddDeleteMoves=size(currentPossibleMoves,1);
        
        currentTotalPossibleMoves=numAddDeleteMoves+2;
        
        % +1 for within model move and +1 for swap move
        %currentModelUpdateWeight=zeros(currentTotalPossibleMoves,1);
        currentModelWeightWithinModel=zeros(nReactions,1);
                
        currentNumTReactions=ones(currentTotalPossibleMoves,1);
        
        if algorithm==2
            % This needs to be updated
            currentModelUpdateWeight=currentNumTReactions;
            
        end
        
        % Within model move
        if nnz(currentInfluentialReactions)~=0
            currentModelUpdateWeight(currentTotalPossibleMoves-1)=1;
        else
            currentModelUpdateWeight(currentTotalPossibleMoves-1)=0;
        end
        
        % Swap move
        currentModelUpdateWeight(currentTotalPossibleMoves)=0;
        
        for i=1:nReactions
            if any(i==currentInfluentialReactions)
                currentModelWeightWithinModel(i)=1;
            else
                currentModelWeightWithinModel(i)=0;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end
%%

% add/delete probabilities
currentModelUpdateWeight(1:(currentTotalPossibleMoves-2))=currentModelUpdateWeight(1:(currentTotalPossibleMoves-2))/sum(currentModelUpdateWeight(1:(currentTotalPossibleMoves-2)));

currentModelUpdateWeight(currentTotalPossibleMoves-1)=currentModelUpdateWeight(currentTotalPossibleMoves-1);

currentModelUpdateWeight(currentTotalPossibleMoves)=currentModelUpdateWeight(currentTotalPossibleMoves);

currentModelWeightWithinModel=currentModelWeightWithinModel/sum(currentModelWeightWithinModel);
