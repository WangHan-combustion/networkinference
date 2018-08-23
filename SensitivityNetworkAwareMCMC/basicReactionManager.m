function [speciesPoolTemp,influentialReactions,uninfluentialReactions] = basicReactionManager(stoichiometricMatrix,reactionType,speciesReference,initialIncludedSpeciesPool,operation,reactionExcluded)

% This function determine the set of influential and uninfluential
% reactions by consideration of all active reactions given the
% initialIncluded speices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are constants

nReactions=size(stoichiometricMatrix,2);          % total number of reactions in the problem
nSpecies=size(stoichiometricMatrix,1);            % total number of species in the problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No need to preallocate, influentialReactions=zeros(nReactions,1);         % influential reactions (active reactions)
uninfluentialReactions=zeros(nReactions,1);       % uninfluential reactions (nonactive reactions) 

flag=1;

% Initialize all variables

nReactionsTemp=0;                                           % Number of reactions added
nSpeciesTemp=length(initialIncludedSpeciesPool);            % Number oof species added

reactionPoolTemp=zeros(nReactions,1);                       % This keeps track of all reactions added
speciesPoolTemp=zeros(nSpecies,1);                          % This keeps track of all species added

speciesPoolTemp(1:nSpeciesTemp)=initialIncludedSpeciesPool; % This keeps track of all species added

while flag
        
    sizeInitialReactions=nReactionsTemp;
        
    for r=1:nReactions
        if ~any(reactionExcluded==r) && ~any(reactionPoolTemp==r)
            
            indf=find(stoichiometricMatrix(:,r)==-1|stoichiometricMatrix(:,r)==3);
            indb=find(stoichiometricMatrix(:,r)==1|stoichiometricMatrix(:,r)==3);
            
            s1=size(indf,1);
            s2=size(indb,1);
            pf=zeros(s1,1);
            pb=zeros(s2,1);
            
            for i=1:s1
              if any(speciesReference(indf(i))==speciesPoolTemp)
                  pf(i)=1;
              end
            end
            
            
            for j=1:s2
              if any(speciesReference(indb(j))==speciesPoolTemp)
                  pb(j)=1;
              end
            end

                       
            % Update the reaction pool if the reaction found to be active
            if reactionType(r)==1
             if all(pf) || all(pb)

                reactionPoolTemp(nReactionsTemp+1)=r;
                nReactionsTemp=nReactionsTemp+1;

             end
            else
             if all(pf)

                reactionPoolTemp(nReactionsTemp+1)=r;
                nReactionsTemp=nReactionsTemp+1;
             end
            end
            
        end
    end
    
    sizeFinalReactions=nReactionsTemp;
        
    sizeInitialSpecies=nSpeciesTemp;
    
    for s=1:nReactionsTemp
        
        % ind is a column vector because stoichiometricMatrix(:,x) is a
        % column vector
        ind=find(stoichiometricMatrix(:,reactionPoolTemp(s)));
                
        % Update the species pool if species not already present
        for k=1:size(ind,1) 
           if ~any(speciesReference(ind(k))==speciesPoolTemp) 

             speciesPoolTemp(nSpeciesTemp+1)=speciesReference(ind(k));
             nSpeciesTemp=nSpeciesTemp+1;
           end 
        end
    end
    
    sizeFinalSpecies=nSpeciesTemp;
    
    if sizeFinalSpecies-sizeInitialSpecies==0 && sizeFinalReactions-sizeInitialReactions==0
        flag=0;
    end
            
end


influentialReactions=reactionPoolTemp;

% Keeping track of all excluded reactions

count=0;
for i=1:nReactions
   if ~any(influentialReactions==i) 
     uninfluentialReactions(count+1)=i;  
     count=count+1;
   end    
end



