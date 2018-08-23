function [currRemovedReactions,key,keyCount,keyMoveDetermination,keyAllModels,modelIndex,nonZeroPriorModel,AllPossibleMoves,AllPossibleNumPaths,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,count,counter]=recurseRemovedReactions(addendum,counter,currRemovedReactions,key,keyCount,keyMoveDetermination,keyAllModels,modelIndex,nonZeroPriorModel,AllPossibleMoves,AllPossibleNumPaths,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,count,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix)

i=0;

while i<=addendum
    
    removedReactions=currRemovedReactions;
    
    if i~=0
        removedReactions(i)=i;
    end
        
    if addendum<nReactions
        [~,key,keyCount,keyMoveDetermination,keyAllModels,modelIndex,nonZeroPriorModel,AllPossibleMoves,AllPossibleNumPaths,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,count,counter]=recurseRemovedReactions(addendum+1,counter,removedReactions,key,keyCount,keyMoveDetermination,keyAllModels,modelIndex,nonZeroPriorModel,AllPossibleMoves,AllPossibleNumPaths,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,count,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix);
    else
               
                
        [influentialReactions,uninfluentialReactions] = identifyInfluentialReactions(initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,removedReactions(removedReactions~=0));                
        
        
                
            [li,loc]=ismember(influentialReactions',key,'rows');
                
            if ~li
                                                           
                key(count+1,:)=influentialReactions';
                keyCount(count+1,1)=1;
                keyMoveDetermination(count+1,1)=0;
                %AllPossibleMoves(count+1)={AllowedMoves};
                %AllPossibleMovesTReactions(count+1)={transDimReactions};
                %AllPossibleMovesNumTReactions(count+1)={numTransDimReactions};
                count=count+1;
            
            else
                keyCount(loc,1)=keyCount(loc,1)+1;                
            end

        
        counter=counter+1
        
        for n=1:nReactions
         if any(removedReactions==n)
             keyAllModels(counter,n)=0;
         else
             keyAllModels(counter,n)=n;
         end    
        end
        
        if ~li
            modelIndex(counter)=count;
        else
            modelIndex(counter)=loc;
        end
        
        if nnz(influentialReactions)>nReactions/2
           nonZeroPriorModel(counter)=1;
        else
           nonZeroPriorModel(counter)=0; 
        end     
    end
    
    i=i+addendum;
end