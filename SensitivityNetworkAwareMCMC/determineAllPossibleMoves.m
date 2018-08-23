function [key,keyCount,keyMoveDetermination,keyAllModels,modelIndex,nonZeroPriorModel,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,count] = determineAllPossibleMoves(typeSampling,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix)

% Number of reactions

count=0;

AllPossibleMoves=cell(1,1);
AllPossibleMovesTReactions=cell(1,1);
AllPossibleMovesNumTReactions=cell(1,1);
AllPossibleNumPaths=cell(1,1);
AllPossibleNumOfModelsForEachMove=cell(1,1);
key=ones(1,nReactions);
keyCount=zeros(1);
keyMoveDetermination=zeros(1);
currRemovedReactions=zeros(nReactions,1);

keyAllModels=ones(1,nReactions);
modelIndex=zeros(1);    
nonZeroPriorModel=zeros(1);

counter=0;

addendum=1;
i=0;

tic

if typeSampling==0 || typeSampling==1 || typeSampling==3
    
    while i<=addendum
        
        removedReactions=currRemovedReactions;
        
        if i~=0
            removedReactions(i)=i;
        end
        
        [~,key,keyCount,keyMoveDetermination,keyAllModels,modelIndex,nonZeroPriorModel,AllPossibleMoves,AllPossibleNumPaths,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,count,counter]=recurseRemovedReactions(addendum+1,counter,removedReactions,key,keyCount,keyMoveDetermination,keyAllModels,modelIndex,nonZeroPriorModel,AllPossibleMoves,AllPossibleNumPaths,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,count,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix);
        
        i=i+addendum;
        
    end
    
end

if typeSampling==3
    
    keyMoveDetermination=zeros(size(keyAllModels,1),1);
    
end

keyAllModels=flipud(keyAllModels);
modelIndex=fliplr(modelIndex);
nonZeroPriorModel=fliplr(nonZeroPriorModel);
