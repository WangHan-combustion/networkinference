function [pathExist] = advancedReactionManager(reactionOfInterest,stoichiometricMatrix,speciesReference,includedReactionPool,targetNode)

pathExist=0;

nReactions=size(stoichiometricMatrix,2);       % number of reactions in the problem
nSpecies=size(stoichiometricMatrix,1);         % number of species in the problem

reactionPoolTemp=zeros(nReactions,1);
speciesPoolTemp=zeros(nSpecies,1);

nReactionsTemp=0;
nSpeciesTemp=0;


%% Include the reaction of interest in the reaction vector (reactionPoolTemp) and corresponding species in the species vector (speciesPoolTemp)

ind=find(stoichiometricMatrix(:,reactionOfInterest)==1|stoichiometricMatrix(:,reactionOfInterest)==-1);

for k=1:size(ind,1)
    speciesPoolTemp(nSpeciesTemp+1)=speciesReference(ind(k));
    nSpeciesTemp=nSpeciesTemp+1;
end

flag=1;

reactionPoolTemp(nReactionsTemp+1)=reactionOfInterest;
nReactionsTemp=nReactionsTemp+1;
    
while flag
    
    for r=1:length(includedReactionPool)
                
        if ~any(reactionPoolTemp==includedReactionPool(r))  % Look for reactions in the includedReactionPool which are not in the reactionPoolTemp
                                         
            %
            indR=find(stoichiometricMatrix(:,includedReactionPool(r))==-1 | stoichiometricMatrix(:,includedReactionPool(r))==3);
            indP=find(stoichiometricMatrix(:,includedReactionPool(r))==1);
            % 
            
            sizeReactants=length(indR);   % Number of reactants added
            sizeProducts=length(indP);    % Number of products added
            
            pReactants=zeros(sizeReactants,1);
            pProducts=zeros(sizeProducts,1);
                        
            %
             for i=1:sizeReactants
              
                if any(speciesReference(indR(i))==speciesPoolTemp)
                    pReactants(i)=1;
                end
                
             end
            
            for j=1:sizeProducts

                if any(speciesReference(indP(j))==speciesPoolTemp)
                    pProducts(j)=1;
                end                
                
            end           
            %
                               
            %
            if any(pReactants) || any(pProducts)
                reactionPoolTemp(nReactionsTemp+1)=includedReactionPool(r);
                nReactionsTemp=nReactionsTemp+1;
            end             
            % 
                      
        end
        
    end
            
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Test if the while loop should stop
    sizeInitialSpecies=nSpeciesTemp;
            
    for s=1:nReactionsTemp
        ind=find(stoichiometricMatrix(:,reactionPoolTemp(s))==1|stoichiometricMatrix(:,reactionPoolTemp(s))==-1);
        
        for k=1:size(ind,1)
           if ~any(speciesReference(ind(k))==speciesPoolTemp) 
             speciesPoolTemp(nSpeciesTemp+1)=speciesReference(ind(k));
             nSpeciesTemp=nSpeciesTemp+1;
           end 
        end
    end
    
    sizeFinalSpecies=nSpeciesTemp;
    
    if sizeFinalSpecies-sizeInitialSpecies==0
        flag=0;                       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

    includedSpeciesPool=speciesPoolTemp;

    if any(includedSpeciesPool==targetNode)
        pathExist=1;
    end
    
