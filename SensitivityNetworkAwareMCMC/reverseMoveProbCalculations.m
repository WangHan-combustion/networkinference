function [reversePossibleMoves,reverseModelUpdateWeight,reverseModelWeightWithinModel,reverseNumTReactions,locationOfProposed,reverseMoveNumber,reverseWithinModelMoveNumber,reactionRatesDRGReactionsRev,productionRev,consumptionRev,rDRGEPReactionsRev,rDRGEPRemovedReactionsRev,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,proposedKeyLocation] = reverseMoveProbCalculations(currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,currentMoveNumber,currentWithinModelMoveNumber,withinModel,y1,priorRange,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,modelWeightScalingRegular,RJMCMCType,algorithm,typeSampling,OnTheFlyModelDetermination,AddDelMoveProb,WithinMoveProb,SwapMoveProb,inputX,key,keyAllModels,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,nonZeroPriorModel,wVector,K,K1,modelIndex,proposedKeyLocation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reverse move probability calculations

% +1 is for within model move +1 for swap model move
scaleLength=size(modelWeightScalingRegular,1);

reactionRatesDRGReactionsRev=0;
productionRev=0;
consumptionRev=0;
rDRGEPReactionsRev=0;

rDRGEPRemovedReactionsRev=zeros(nReactions,nReactions);

reversePossibleMoves=0;

%default
locationOfProposed=0;

%% Chemical kinetic reversible-jump MCMC with model check (RJMCMC Type 2)
if RJMCMCType==2 || RJMCMCType==4
    
    if typeSampling==1

        proposedTotalPossibleMoves=nReactions+2;
        
        reverseModelUpdateWeight=zeros(proposedTotalPossibleMoves,1);
        reverseModelWeightWithinModel=zeros(nReactions,1);
        
        reverseNumTReactions=zeros(proposedTotalPossibleMoves,1);
        
        if algorithm==2
            % Add/Delete move
            for i=1:nReactions
                
               if i~=7 && i~=2 && i~=1 && i~=9 && i~=4 && i~=10 && i~=12 %&& i~=8 && i~=11 && i~=13%&& i~=1 %&& i~=10 && i~=12%&& i~=1 && i~=2 && i~=12 && i~=10 && i~=7 && i~=6
               %if i~=2
                reverseNumTReactions(i)=1;
                reverseModelUpdateWeight(i)=1;
               end 
                
            end
        elseif algorithm==3
            % Add/Delete move
            for i=1:nReactions
                reverseNumTReactions(i)=1;
                
                % If reaction is currently prsent, it will be removed else it
                % will be added
                if any(proposedInfluentialReactions==i)
                    reverseModelUpdateWeight(i)=wVector(K1,i);
                else
                    reverseModelUpdateWeight(i)=wVector(K1+1,i);
                end
                
            end
        end
        
        % Within model move
        if nnz(intermediateProposedInfluentialReactions)~=0
            reverseModelUpdateWeight(nReactions+1)=1;
        else
            reverseModelUpdateWeight(nReactions+1)=0;
        end
        
        % Swap move
        reverseModelUpdateWeight(nReactions+2)=0;
        
        for i=1:nReactions
            if any(i==intermediateProposedInfluentialReactions)
                reverseModelWeightWithinModel(i)=1;
            else
                reverseModelWeightWithinModel(i)=0;
            end
        end
        
        reverseMoveNumber=currentMoveNumber;
        reverseWithinModelMoveNumber=currentWithinModelMoveNumber;        
        
    elseif typeSampling==0
        
        proposedTotalPossibleMoves=nReactions+2;
        
        reverseModelUpdateWeight=zeros(proposedTotalPossibleMoves,1);
        reverseModelWeightWithinModel=zeros(nReactions,1);
        
        reverseNumTReactions=zeros(proposedTotalPossibleMoves,1);
        
        if algorithm==2
            % Add/Delete move
            for i=1:nReactions
                
               if i~=7 && i~=2 && i~=1 && i~=9 && i~=4 && i~=10 && i~=12 %&& i~=8 && i~=11 && i~=13%&& i~=1 %&& i~=10 && i~=12%&& i~=1 && i~=2 && i~=12 && i~=10 && i~=7 && i~=6
               %if i~=2
                reverseNumTReactions(i)=1;
                reverseModelUpdateWeight(i)=1;
               end 
                
            end
        elseif algorithm==3
            % Add/Delete move
            for i=1:nReactions
                reverseNumTReactions(i)=1;
                
                % If reaction is currently prsent, it will be removed else it
                % will be added
                if any(proposedInfluentialReactions==i)
                    reverseModelUpdateWeight(i)=wVector(K1,i);
                else
                    reverseModelUpdateWeight(i)=wVector(K1+1,i);
                end
                
            end
        end
        
        % Within model move
        if nnz(proposedInfluentialReactions)~=0
            reverseModelUpdateWeight(nReactions+1)=1;
        else
            reverseModelUpdateWeight(nReactions+1)=0;
        end
        
        % Swap move
        reverseModelUpdateWeight(nReactions+2)=0;
        
        for i=1:nReactions
            if any(i==proposedInfluentialReactions)
                reverseModelWeightWithinModel(i)=1;
            else
                reverseModelWeightWithinModel(i)=0;
            end
        end
        
        reverseMoveNumber=currentMoveNumber;
        reverseWithinModelMoveNumber=currentWithinModelMoveNumber;
        
    elseif typeSampling==3
        
        index=proposedKeyLocation;
               
        % addition/deletion move
        % determine the move
        
        if withinModel==0
                                    
            proposedUninfluentialReactions=(1:nReactions)';
            
            for j=1:nReactions
               if any(proposedInfluentialReactions==j)
                   proposedUninfluentialReactions(j)=0;
               end                    
            end 
                        
            if OnTheFlyModelDetermination
                if keyMoveDetermination(index)==0
                    
                    [reactionSet,AllowedMoves,numPossiblePaths,numOfModelsForEachMove,~,transDimReactions,numTransDimReactions] = determineMoves(typeSampling,key,keyCount,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,proposedInfluentialReactions,proposedUninfluentialReactions);
                    
                    AllPossibleMoves(index)={AllowedMoves};
                    AllPossibleNumPaths(index)={numPossiblePaths};
                    AllPossibleMovesTReactions(index)={transDimReactions};
                    AllPossibleMovesNumTReactions(index)={numTransDimReactions};
                    AllPossibleNumOfModelsForEachMove(index)={numOfModelsForEachMove};                    
                    keyMoveDetermination(index)=1;
                end
            end
            
        end
        
        
        locationOfProposed=index;
        
        reversePossibleMoves=AllPossibleMoves{index};
        numAddDeleteMoves=size(reversePossibleMoves,1);
        
        proposedTotalPossibleMoves=numAddDeleteMoves+2;
        
        % +1 for within model move and +1 for swap move
        reverseModelWeightWithinModel=zeros(nReactions,1);
        
        reverseNumTReactions=AllPossibleNumPaths{index};
        
        reverseNumTReactions=ones(proposedTotalPossibleMoves,1);
        
        if algorithm==2
            
            reverseModelUpdateWeight=reverseNumTReactions;
            
        elseif algorithm==3
            
            reverseModelUpdateWeight=ones(1,proposedTotalPossibleMoves);
            
            % Update this using the adaptive component weights
            for i=1:numAddDeleteMoves
                for j=1:nReactions
                    if any(reversePossibleMoves(i,:)==j)
                        reverseModelUpdateWeight(i)=reverseModelUpdateWeight(i)*wVector(K1+1,j);
                    else
                        reverseModelUpdateWeight(i)=reverseModelUpdateWeight(i)*wVector(K1,j);
                    end
                end
            end
            
        end       
                                
        % Within model move
        if nnz(proposedInfluentialReactions)~=0
            reverseModelUpdateWeight(numAddDeleteMoves+1)=1;
        else
            reverseModelUpdateWeight(numAddDeleteMoves+1)=0;
        end
        
        % Swap move
        reverseModelUpdateWeight(numAddDeleteMoves+2)=0;
        
        for i=1:nReactions
            if any(i==proposedInfluentialReactions)
                reverseModelWeightWithinModel(i)=1;
            else
                reverseModelWeightWithinModel(i)=0;
            end
        end
        
        if withinModel==0
           
           for i=1:numAddDeleteMoves 
             if reversePossibleMoves(i,:)==currentInfluentialReactions'
                 reverseMoveNumber=i;
                 break;
             end
           end    
               
        elseif withinModel==1
           reverseMoveNumber=numAddDeleteMoves+1; 
        elseif withinModel==2
           reverseMoveNumber=numAddDeleteMoves+2; 
        end    
            
        
        reverseWithinModelMoveNumber=currentWithinModelMoveNumber;
        
        
    end
    
end

% add/delete probabilities
reverseModelUpdateWeight(1:(proposedTotalPossibleMoves-2))=reverseModelUpdateWeight(1:(proposedTotalPossibleMoves-2))/sum(reverseModelUpdateWeight(1:(proposedTotalPossibleMoves-2)));

reverseModelUpdateWeight(proposedTotalPossibleMoves-1)=reverseModelUpdateWeight(proposedTotalPossibleMoves-1);

reverseModelUpdateWeight(proposedTotalPossibleMoves)=reverseModelUpdateWeight(proposedTotalPossibleMoves);

reverseModelWeightWithinModel=reverseModelWeightWithinModel/sum(reverseModelWeightWithinModel);

    