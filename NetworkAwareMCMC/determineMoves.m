function [reactionSet,AllowedMoves,numPossiblePaths,numOfModelsForEachMove,countAllowedMoves,transDimReactions,numTransDimReactions] = determineMoves(typeSampling,key,keyCount,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,influentialReactions,uninfluentialReactions)

influentialReactions'

uninfluentialReactions'

% Determine all possible legitimate moves

%% Main outputs: AllowedMoves               - all possible moves from given influential reaction set
%%               countAllowedMoves          - number of possible moves
%%               transDimReactions          - reactions that can act as the final reaction
%%               numPossiblePaths           - number of ways to get to specific cluster from the given cluster
%%               numOfModelsForEachMove     - number of models in the cluster of the proposed move

%% Main inputs:  influential reactions   - reactions that influential to the target
%%               uninfluential reactions - reactions that are not influential to the target

% addition/deletion move

if typeSampling ==1
    
    countAllowedMoves=0;
    AllowedMoves=zeros(nReactions+1,nReactions);
    numOfModelsForEachMove=zeros(nReactions+1,1);
    numPossiblePaths=zeros(2^nReactions,1);
    transDimReactions=zeros(nReactions+1,nReactions,2);
    
    numInfluentialReactions=nnz(influentialReactions);
    numUninfluentialReactions=nnz(uninfluentialReactions);
    
    % factorMult is same as keyCount
    factorMult=1;
    
    numReactionsToBeAdded=5;
    
    %{
    % Adding one possible "AllowedMoves" as the self move
    
    % Here we just declare self move as a possible move
    AllowedMoves(countAllowedMoves+1,:)=influentialReactions';
    numPossiblePaths(countAllowedMoves+1)=numPossiblePaths(countAllowedMoves+1)+1;
                    
    for keyLoopCounter=1:size(key,1)
          if key(keyLoopCounter,:)==influentialReactions'
             index=keyLoopCounter;
             break;
          end
    end
                    
    numOfModelsForEachMove(countAllowedMoves+1)=keyCount(index);
                        
    countAllowedMoves=countAllowedMoves+1;
    %}
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    reactionSize=1;
    numOfElements=zeros(nReactions,1);
    
    % (count-1) is the number of elements in the reaction set
    count=1;
    
    % Reaction set refers to the set of reactions that are not sufficient to
    % alter the model
    reactionSet={[]};
    
    if numUninfluentialReactions~=0
        for i=1:numUninfluentialReactions
            
            % create the reactionExcluded vector corresponding to each of the
            % reaction added from the uninfluential reactions
            
            reactionExcluded=uninfluentialReactions';
            reactionExcluded(reactionExcluded==uninfluentialReactions(numInfluentialReactions+i))=0;
            reactionExcluded=reactionExcluded(reactionExcluded~=0);
            
            [influentialReactionsTemp,uninfluentialReactionsTemp]=identifyInfluentialReactions(initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,reactionExcluded);
            
            if all(influentialReactions==influentialReactionsTemp)
                
                % Still in the same cluster
                factorMult=factorMult+1;
                
                reactionSet(reactionSize,count)={sort(uninfluentialReactions(numInfluentialReactions+i))};
                count=count+1;
                
                %{
                transDimReactions(1,nReactions+1,1)=nReactions+1;
                transDimReactions(1,nReactions+1,2)=transDimReactions(1,nReactions+1,2)+1;
                
                % Adding the current uninfluential reaction %
                               
                numPossiblePaths(1)=numPossiblePaths(1)+2;
                %}
                
            else
                
                flag1=0;
                
                for k=1:countAllowedMoves
                    if all(AllowedMoves(k,:)==influentialReactionsTemp')
                        flag1=1;
                        transDimReactions(k,uninfluentialReactions(numInfluentialReactions+i),1)=uninfluentialReactions(numInfluentialReactions+i);
                        transDimReactions(k,uninfluentialReactions(numInfluentialReactions+i),2)=transDimReactions(k,uninfluentialReactions(numInfluentialReactions+i),2)+1;
                        numPossiblePaths(k)=numPossiblePaths(k)+1;                        
                    end
                end
                
                if flag1==0
                    AllowedMoves(countAllowedMoves+1,:)=influentialReactionsTemp';
                    transDimReactions(countAllowedMoves+1,uninfluentialReactions(numInfluentialReactions+i),1)=uninfluentialReactions(numInfluentialReactions+i);
                    transDimReactions(countAllowedMoves+1,uninfluentialReactions(numInfluentialReactions+i),2)=transDimReactions(countAllowedMoves+1,uninfluentialReactions(numInfluentialReactions+i),2)+1;
                    numPossiblePaths(countAllowedMoves+1)=numPossiblePaths(countAllowedMoves+1)+1;
                    
                    for keyLoopCounter=1:size(key,1)
                       if key(keyLoopCounter,:)==influentialReactionsTemp'
                         index=keyLoopCounter;
                         break;
                       end
                    end
                    
                    numOfModelsForEachMove(countAllowedMoves+1)=keyCount(index);
                    
                    countAllowedMoves=countAllowedMoves+1;                    
                end
                
                
            end
            
        end
        
        numOfElements(reactionSize)=count-1;
        
        
        for reactionSize=1:numReactionsToBeAdded-1
            if count>1
                reactionSize
                [factorMult,count,reactionSet,numOfElements,AllowedMoves,numPossiblePaths,numOfModelsForEachMove,transDimReactions,countAllowedMoves]=newPoolReactions(key,keyCount,reactionSize+1,factorMult,numOfElements,AllowedMoves,numPossiblePaths,numOfModelsForEachMove,transDimReactions,countAllowedMoves,reactionSet,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,influentialReactions,uninfluentialReactions);
            end
        end
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If reaction i is part of the influential reactions
    for i=1:numInfluentialReactions
        
        % Determine the new set of uninfluentialReactions
        reactionExcluded=uninfluentialReactions';
        % This is fine because not all reactions can be
        % uninfluential
        reactionExcluded(1)=influentialReactions(numUninfluentialReactions+i);
        
        % influentialReactionsTemp tells me which cluster I have moved into
        [influentialReactionsTemp,uninfluentialReactionsTemp]=identifyInfluentialReactions(initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,reactionExcluded(reactionExcluded~=0));
        
        if (nnz(influentialReactions)-nnz(influentialReactionsTemp))<=numReactionsToBeAdded
            flag1=0;
            
            % Check if the new set of uninfluentialReactions is already
            % part of allowed moves
            
            for k=1:countAllowedMoves
                if all(AllowedMoves(k,:)==influentialReactionsTemp')
                    flag1=1;
                    transDimReactions(k,influentialReactions(numUninfluentialReactions+i),1)=influentialReactions(numUninfluentialReactions+i);
                    transDimReactions(k,influentialReactions(numUninfluentialReactions+i),2)=transDimReactions(k,influentialReactions(numUninfluentialReactions+i),2)+1;
                    numPossiblePaths(k)=numPossiblePaths(k)+1;
                end
            end
            
            % If the new influential reactions vector is not already in
            % AllowedMoves then include it
            if flag1==0
                
                AllowedMoves(countAllowedMoves+1,:)=influentialReactionsTemp';
                transDimReactions(countAllowedMoves+1,influentialReactions(numUninfluentialReactions+i),1)=influentialReactions(numUninfluentialReactions+i);
                transDimReactions(countAllowedMoves+1,influentialReactions(numUninfluentialReactions+i),2)=transDimReactions(countAllowedMoves+1,influentialReactions(numUninfluentialReactions+i),2)+1;
                numPossiblePaths(countAllowedMoves+1)=numPossiblePaths(countAllowedMoves+1)+1;
                
                for keyLoopCounter=1:size(key,1)
                   if key(keyLoopCounter,:)==influentialReactionsTemp'
                      index=keyLoopCounter;
                      break;
                   end
                end
                    
                numOfModelsForEachMove(countAllowedMoves+1)=keyCount(index);                
                                                
                countAllowedMoves=countAllowedMoves+1;
                
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    AllowedMoves=AllowedMoves(1:countAllowedMoves,:);
    
    % This is different from the function in reactionManager
    transDimReactions=transDimReactions(1:countAllowedMoves,:,1);
    
    numPossiblePaths=numPossiblePaths(1:countAllowedMoves);
    
    numOfModelsForEachMove=numOfModelsForEachMove(1:countAllowedMoves);
    
    numTransDimReactions=zeros(countAllowedMoves,1);
    
    for i=1:countAllowedMoves
        numTransDimReactions(i)=nnz(transDimReactions(i,:,1));
    end
    
elseif typeSampling==3
    
    reactionSet={[]};
    
    countAllowedMoves=0;
    AllowedMoves=zeros(nReactions+1,nReactions);
    numPossiblePaths=zeros(nReactions+1,1);
    transDimReactions=zeros(nReactions+1,nReactions+1,2);
        
    numReactionsToBeAdded=4;
    
    reactionExcluded=sort(uninfluentialReactions);
    reactionExcluded=reactionExcluded(reactionExcluded~=0);
    
    [influentialReactionsTemp1,uninfluentialReactionsTemp1]=identifyInfluentialReactions(initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,reactionExcluded);
    
    for i=1:nReactions
        
        influentialReactionsTemp=influentialReactions;
        uninfluentialReactionsTemp=uninfluentialReactions;
        
        if any(i==influentialReactions)
            influentialReactionsTemp(i)=0;
            uninfluentialReactionsTemp(i)=i;
        else
            influentialReactionsTemp(i)=i;
            uninfluentialReactionsTemp(i)=0;
        end
        
        reactionExcluded=sort(uninfluentialReactionsTemp);
        reactionExcluded=reactionExcluded(reactionExcluded~=0);
        
        [influentialReactionsTemp2,uninfluentialReactionsTemp2]=identifyInfluentialReactions(initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,reactionExcluded);
        
        if abs(nnz(influentialReactionsTemp1)-nnz(influentialReactionsTemp2))<=numReactionsToBeAdded
            
            AllowedMoves(countAllowedMoves+1,:)=influentialReactionsTemp';
            %%%%% Not used
            transDimReactions(countAllowedMoves+1,i,1)=i;
            transDimReactions(countAllowedMoves+1,i,2)=transDimReactions(countAllowedMoves+1,i,2)+1;
            %%%%%
            numPossiblePaths(countAllowedMoves+1)=numPossiblePaths(countAllowedMoves+1)+1;
            countAllowedMoves=countAllowedMoves+1;
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    AllowedMoves=AllowedMoves(1:countAllowedMoves,:);
    
    % This is different from the function in reactionManager
    transDimReactions=transDimReactions(1:countAllowedMoves,:,1);
      
    numPossiblePaths=numPossiblePaths(1:countAllowedMoves);
    
    numTransDimReactions=zeros(countAllowedMoves,1);
    
    for i=1:countAllowedMoves
        numTransDimReactions(i)=1;
    end
    
end
