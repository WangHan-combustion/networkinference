function [influentialReactions,uninfluentialReactions] = identifyInfluentialReactions(initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,removedReactions)

%%
% includedSpeciesPool: the set of species currently in the model

% removedReactions: the set of reactions that have been explicitly turned off

% Assume removedReactions is a set of non zero entries, no zero values in
% the vector

% includedReactionPool is determined based on the includedSpeciesPool and
% removedReactions

%%                                      
                    
%%                   

%% Definitions:

% Outputs:

% influentialReactions: Set of reactions that DO currently have an influence on the output
% uninfluentialReactions: All reactions that DO NOT currently have an influence on the output

% Inputs:

% removedReactions

%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  BASIC REACTION MANAGER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the currently included reactions

% the type of model (unsteady/steady)
operation=1;

% reactions to kept excluded
% reactionExcluded and removedReactions not padded with zeros
reactionExcluded=removedReactions;

[~,influentialReactions,uninfluentialReactions] = basicReactionManager(stoichiometricMatrix,reactionType,speciesReference,initialIncludedSpeciesPool,operation,reactionExcluded);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADVANCED REACTION MANAGER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% removing other reactions that do not matter from the set of included
% reactions (This is applicable if there are catalyst)

reactionElimination=zeros(nnz(influentialReactions),1);
numReactionElimination=0;
    
numReactionsAlreadyRemoved=nnz(uninfluentialReactions);
        
% We test here if any further reactions can be removed
        
for i=1:nnz(influentialReactions)
        
        %  This is fine because in basicManager reactions are added from
        %  the left
        reactionOfInterest=influentialReactions(i);
                
        [~,~,includedReactionPool]=find(influentialReactions);
        [pathExist] = advancedReactionManager(reactionOfInterest,stoichiometricMatrix,speciesReference,includedReactionPool,targetNode);
        
        if pathExist==0
            reactionElimination(numReactionElimination+1)=reactionOfInterest;
            numReactionElimination=numReactionElimination+1;
        end
        
end
    
for j=1:numReactionElimination
        
   influentialReactions(influentialReactions==reactionElimination(j))=0;        
   uninfluentialReactions(numReactionsAlreadyRemoved+j)=reactionElimination(j);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT INFLUENTIAL AND UNINFLUENTIAL REACTIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uninfluentialReactions=sort(uninfluentialReactions);
influentialReactions=sort(influentialReactions);


