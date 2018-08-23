function [currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,currentMoveNumber,currentWithinModelMoveNumber,AMoveNum,WMoveNum,SMoveNum,withinModel,AddDeleteReaction,AddDeleteReactionSwap,proposedModelIndex,proposedKeyLocation]=sampleProposalMove(typeSampling,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,AddDelMoveProb,WithinMoveProb,SwapMoveProb,modelWeightScalingRegular,modelWeightScalingSwap,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,currentInfluentialReactions,intermediateCurrentInfluentialReactions,reactionRatesDRGReactions,production,consumption,AMoveNum,WMoveNum,SMoveNum,RJMCMCType,rDRGEPReactions,rDRGEPRemovedReactions,currentPossibleMoves,keyAllModels,modelIndex,keyCount,key,currentModelIndex,currentKeyLocation)

AddDeleteReactionSwap=0;

totalPossibleMoves=length(currentModelUpdateWeight);

% +1 for within model move +1 for swap move
cumulativeModelWeightVector=zeros(totalPossibleMoves,1);
cumulativeModelWeightWithinModel=zeros(nReactions,1);

proposedModelIndex=0;
proposedKeyLocation=0;

currentMoveNumber=0;
currentWithinModelMoveNumber=0;

% Cumulating probabilities of add/del moves
for nAll=1:totalPossibleMoves-2
    if nAll==1
        cumulativeModelWeightVector(nAll)=currentModelUpdateWeight(nAll);
    else
        cumulativeModelWeightVector(nAll)=cumulativeModelWeightVector(nAll-1)+currentModelUpdateWeight(nAll);
    end
end

% Cumulating probabilities all within model move types
for nOnlyR=1:nReactions
    if nOnlyR==1
        cumulativeModelWeightWithinModel(nOnlyR)=currentModelWeightWithinModel(nOnlyR);
    else
        cumulativeModelWeightWithinModel(nOnlyR)=cumulativeModelWeightWithinModel(nOnlyR-1)+currentModelWeightWithinModel(nOnlyR);
    end
end


for n=1:nReactions
    if n<=2
        cumulativeModelWeightVector(n)=-10;
    end
end


moveType=rand(1);

% Model Check algorithm
if typeSampling==1

    proposedInfluentialReactions=currentInfluentialReactions;
    
    if moveType<=AddDelMoveProb
        
        a=rand(1);
        
        for i=1:totalPossibleMoves-2
            if i==1
                % add/delete move
                if a<=cumulativeModelWeightVector(i)
                 %if i==1   
                    currentMoveNumber=i;
                                        
                    currentTReaction=currentMoveNumber;
                    
                    withinModel=0;
                    AddDeleteReaction=1;
                    AMoveNum=AMoveNum+1;
                    
                    if any(currentTReaction==currentInfluentialReactions)
                        
                        % delete reaction
                        proposedInfluentialReactions(currentTReaction)=0;
                        
                    else
                        
                        % add reaction
                        proposedInfluentialReactions(currentTReaction)=currentTReaction;
                        
                    end
                    
                  %  break;
                    
                end
            elseif i<=(totalPossibleMoves-2)                
                % add/delete move
                if a>cumulativeModelWeightVector(i-1) && a<=cumulativeModelWeightVector(i)
                    
                    currentMoveNumber=i;
                                        
                    currentTReaction=currentMoveNumber;
                    
                    withinModel=0;
                    AddDeleteReaction=1;
                    AMoveNum=AMoveNum+1;
                    
                    if any(currentTReaction==currentInfluentialReactions)
                        
                        % delete reaction
                        proposedInfluentialReactions(currentTReaction)=0;
                        
                    else
                        
                        % add reaction
                        proposedInfluentialReactions(currentTReaction)=currentTReaction;
                        
                    end
                    
                end
                
            end
        end

        proposedKeyLocation=1;
        for j=nReactions:-1:1
            proposedKeyLocation=proposedKeyLocation+(proposedInfluentialReactions(j)~=0)*2^(nReactions-j);
        end    
        
        %% Determine the location of the key
        %[~,proposedKeyLocation]=ismember(proposedInfluentialReactions',keyAllModels,'rows');
        
        intermediateProposedInfluentialReactions=key(modelIndex(proposedKeyLocation),:)';        
        
    else
        
        proposedInfluentialReactions=currentInfluentialReactions;
        intermediateProposedInfluentialReactions=intermediateCurrentInfluentialReactions;    

        currentMoveNumber=totalPossibleMoves-1;
        % within model move
        bb=rand(1);
        
        for j=1:nReactions
            if j==1
                if bb<=cumulativeModelWeightWithinModel(j)
                    currentWithinModelMoveNumber=j;
                end
            else
                if bb>cumulativeModelWeightWithinModel(j-1) && bb<=cumulativeModelWeightWithinModel(j)
                    currentWithinModelMoveNumber=j;
                end
            end
        end
        
        withinModel=1;
        WMoveNum=WMoveNum+1;
        
        AddDeleteReaction=2;
        
        
    end
    
elseif typeSampling==0
    
    proposedInfluentialReactions=currentInfluentialReactions;
    
    if moveType<=AddDelMoveProb
        
        a=rand(1);
        
        for i=1:totalPossibleMoves-2
            if i==1
                % add/delete move
                if a<=cumulativeModelWeightVector(i)
                    
                    currentMoveNumber=i;
                                        
                    currentTReaction=currentMoveNumber;
                    
                    withinModel=0;
                    AddDeleteReaction=1;
                    AMoveNum=AMoveNum+1;
                    
                    if any(currentTReaction==currentInfluentialReactions)
                        
                        % delete reaction
                        proposedInfluentialReactions(currentTReaction)=0;
                        
                    else
                        
                        % add reaction
                        proposedInfluentialReactions(currentTReaction)=currentTReaction;
                        
                    end
                    
                end
            elseif i<=(totalPossibleMoves-2)
                
                % add/delete move
                if a>cumulativeModelWeightVector(i-1) && a<=cumulativeModelWeightVector(i)
                    
                    currentMoveNumber=i;
                                        
                    currentTReaction=currentMoveNumber;
                    
                    withinModel=0;
                    AddDeleteReaction=1;
                    AMoveNum=AMoveNum+1;
                    
                    if any(currentTReaction==currentInfluentialReactions)
                        
                        % delete reaction
                        proposedInfluentialReactions(currentTReaction)=0;
                        
                    else
                        
                        % add reaction
                        proposedInfluentialReactions(currentTReaction)=currentTReaction;
                        
                    end
                    
                end
                
            end
        end
    else
        
        proposedInfluentialReactions=currentInfluentialReactions;
        currentMoveNumber=totalPossibleMoves-1;
        % within model move
        bb=rand(1);
        
        for j=1:nReactions
            if j==1
                if bb<=cumulativeModelWeightWithinModel(j)
                    currentWithinModelMoveNumber=j;
                end
            else
                if bb>cumulativeModelWeightWithinModel(j-1) && bb<=cumulativeModelWeightWithinModel(j)
                    currentWithinModelMoveNumber=j;
                end
            end
        end
        
        withinModel=1;
        WMoveNum=WMoveNum+1;
        
        AddDeleteReaction=2;
        
        
    end

    % In this case the intermediate current and proposed influential
    % reactions are same as the overall current and proposed influential
    % reactions
    intermediateCurrentInfluentialReactions=currentInfluentialReactions;
    intermediateProposedInfluentialReactions=proposedInfluentialReactions;
    
    
elseif typeSampling==3
    
    if moveType<=AddDelMoveProb
                
        a=rand(1);
        
        for i=1:totalPossibleMoves-2
            if i==1
                % add/delete move
                if a<=cumulativeModelWeightVector(i)
                    
                    currentMoveNumber=i;
                    
                    % current T reaction is the move number
                                        
                    withinModel=0;
                    AddDeleteReaction=1;
                    AMoveNum=AMoveNum+1;
                    
                    proposedInfluentialReactions=currentPossibleMoves(currentMoveNumber,:)';
                    
                    
                end
            elseif i<=(totalPossibleMoves-2)
                
                % add/delete move
                if a>cumulativeModelWeightVector(i-1) && a<=cumulativeModelWeightVector(i)
                    
                    currentMoveNumber=i;
                                        
                    withinModel=0;
                    AddDeleteReaction=1;
                    AMoveNum=AMoveNum+1;
                    
                    proposedInfluentialReactions=currentPossibleMoves(currentMoveNumber,:)';                    
                end
                
            end
        end
        
        % Determine the model indices of the current and proposed
        % influential reactions
        
        %{
        
        for j=1:size(keyAllModels,1)
            
            if currentInfluentialReactions'==keyAllModels(j,:)                
                currentModelIndex=modelIndex(j);
                currentKeyLocation=j;
                break;                
            end
            
        end
 
        %}
        
        
        for j=1:size(keyAllModels,1)
            
            if proposedInfluentialReactions'==keyAllModels(j,:)                
                proposedModelIndex=modelIndex(j);
                proposedKeyLocation=j;
                break;                
            end
            
        end 

        
        intermediateCurrentInfluentialReactions=currentInfluentialReactions;
        intermediateProposedInfluentialReactions=proposedInfluentialReactions;
        
        
        if currentModelIndex~=proposedModelIndex
            
            
                intermediateCurrentInfluentialReactions=key(currentModelIndex,:)';
                intermediateProposedInfluentialReactions=key(proposedModelIndex,:)';
            
            %{
            for j=1:nReactions
               
               if any(j==currentInfluentialReactions)
                   currentUnInfluentialReactions(j)=0;
               else
                   currentUnInfluentialReactions(j)=j;
               end
               
               if any(j==proposedInfluentialReactions)
                   proposedUnInfluentialReactions(j)=0;
               else
                   proposedUnInfluentialReactions(j)=j;
               end
               
            end
            
            currentUnInfluentialReactions=sort(currentUnInfluentialReactions);
            currentUnInfluentialReactions=currentUnInfluentialReactions(currentUnInfluentialReactions~=0);
            
            proposedUnInfluentialReactions=sort(proposedUnInfluentialReactions);
            proposedUnInfluentialReactions=proposedUnInfluentialReactions(proposedUnInfluentialReactions~=0);
            
            [intermediateCurrentInfluentialReactions,~] = identifyInfluentialReactions(modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,currentUnInfluentialReactions);
            [intermediateProposedInfluentialReactions,~] = identifyInfluentialReactions(modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,proposedUnInfluentialReactions);
            %}
        end

        
    else
        
        
        proposedInfluentialReactions=currentInfluentialReactions;
        proposedModelIndex=currentModelIndex;
        proposedKeyLocation=currentKeyLocation;
        
        
        currentMoveNumber=totalPossibleMoves-1;
        % within model move
        bb=rand(1);
        
        for j=1:nReactions
            if j==1
                if bb<=cumulativeModelWeightWithinModel(j)
                    currentWithinModelMoveNumber=j;
                end
            else
                if bb>cumulativeModelWeightWithinModel(j-1) && bb<=cumulativeModelWeightWithinModel(j)
                    currentWithinModelMoveNumber=j;
                end
            end
        end
        

        
        withinModel=1;
        WMoveNum=WMoveNum+1;
        
        AddDeleteReaction=2;

        intermediateCurrentInfluentialReactions=currentInfluentialReactions;
        intermediateProposedInfluentialReactions=proposedInfluentialReactions;    
        
    end   
        
elseif typeSampling==4    

     %{
     if moveType<=AddDelMoveProb
        
        % Randomly generate the proposed influential reactions vector
                                
        for j=1:size(key,1)
            
            if currentInfluentialReactions'==key(j,:)                
                currentModelIndex=j;
                break;                
            end
            
        end
                
        % Number of models with current model index
        NumberCurrentModel=keyCount(currentModelIndex);
        
        % Location of all models with current model index
        ModelsCurrentModelIndex=find(modelIndex==currentModelIndex);
                        
        % Randomly pick one these models
        selectModel=unidrnd(NumberCurrentModel);
                
        % Moved to the regular space
        intermediateCurrentInfluentialReactions=keyAllModels(ModelsCurrentModelIndex(selectModel),:)';
                                
        intermediateProposedInfluentialReactions=intermediateCurrentInfluentialReactions;
                
        a=rand(1);
        
        for i=1:totalPossibleMoves-2
            if i==1
                % add/delete move
                if a<=cumulativeModelWeightVector(i)
                    
                    currentMoveNumber=i;
                    
                    % current T reaction is the move number
                    
                    currentTReaction=currentMoveNumber;
                    
                    withinModel=0;
                    AddDeleteReaction=1;
                    AMoveNum=AMoveNum+1;
                    
                    if any(currentTReaction==intermediateProposedInfluentialReactions)
                        
                        % delete reaction
                        intermediateProposedInfluentialReactions(currentTReaction)=0;
                        
                    else
                        
                        % add reaction
                        intermediateProposedInfluentialReactions(currentTReaction)=currentTReaction;
                        
                    end
                    
                end
            elseif i<=(totalPossibleMoves-2)
                
                % add/delete move
                if a>cumulativeModelWeightVector(i-1) && a<=cumulativeModelWeightVector(i)
                    
                    currentMoveNumber=i;
                    
                    % current T reaction is the move number
                    
                    currentTReaction=currentMoveNumber;
                    
                    withinModel=0;
                    AddDeleteReaction=1;
                    AMoveNum=AMoveNum+1;
                    
                    if any(currentTReaction==intermediateProposedInfluentialReactions)
                        
                        % delete reaction
                        intermediateProposedInfluentialReactions(currentTReaction)=0;
                        
                    else
                        
                        % add reaction
                        intermediateProposedInfluentialReactions(currentTReaction)=currentTReaction;
                        
                    end
                    
                end
                
            end
        end
        
        
        % Randomly generate the proposed influential reactions vector
        
        for j=1:size(keyAllModels,1)
            
            if intermediateProposedInfluentialReactions'==keyAllModels(j,:)
                
                currentModelIndex=modelIndex(j);
                break;
                
            end
            
        end 
        
        proposedInfluentialReactions=key(currentModelIndex,:)';
                           
    else
                
        for j=1:size(key,1)
            
            if currentInfluentialReactions'==key(j,:)                
                currentModelIndex=j;
                break;                
            end
            
        end
                
        % Number of models with current model index
        NumberCurrentModel=keyCount(currentModelIndex);
        
        % Location of all models with current model index
        ModelsCurrentModelIndex=find(modelIndex==currentModelIndex);
                        
        % Randomly pick one these models
        selectModel=unidrnd(NumberCurrentModel);
                
        intermediateCurrentInfluentialReactions=keyAllModels(ModelsCurrentModelIndex(selectModel),:)';

        intermediateProposedInfluentialReactions=intermediateCurrentInfluentialReactions;
        
        proposedInfluentialReactions=currentInfluentialReactions;
        
        currentMoveNumber=totalPossibleMoves-1;
        % within model move
        bb=rand(1);
        
        for j=1:nReactions
            if j==1
                if bb<=cumulativeModelWeightWithinModel(j)
                    currentWithinModelMoveNumber=j;
                end
            else
                if bb>cumulativeModelWeightWithinModel(j-1) && bb<=cumulativeModelWeightWithinModel(j)
                    currentWithinModelMoveNumber=j;
                end
            end
        end
        
        withinModel=1;
        WMoveNum=WMoveNum+1;
        
        AddDeleteReaction=2;
      
        if nnz(intermediateCurrentInfluentialReactions)==0
            currentWithinModelMoveNumber=0;
        end    
        
             
     end   
     %}
    

    
end

