function [samples,iterAcceptance,iterSamples,ll,logPos,logPosTrue,currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,currentPossibleMoves,currentModelIndex,currentKeyLocation,reactionRatesDRGReactions,production,consumption,AMoveNum,WMoveNum,SMoveNum,AMoveAcc,WMoveAcc,SMoveAcc,rDRGEPReactions,rDRGEPRemovedReactions,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,typeMoves,numKeys,numForwardModelSolves] = generateIndSampler(iterationCount,d1,N,K,K1,delta,ll,logPos,logPosTrue,pos,priorRange,inputX,data,dataVariance,currentInfluentialReactions,intermediateCurrentInfluentialReactions,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,currentPossibleMoves,currentModelIndex,currentKeyLocation,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,modelWeightScalingRegular,modelWeightScalingSwap,RJMCMCType,algorithm,typeSampling,OnTheFlyModelDetermination,reactionRatesDRGReactions,production,consumption,AMoveNum,WMoveNum,SMoveNum,AMoveAcc,WMoveAcc,SMoveAcc,AddDelMoveProb,WithinMoveProb,SwapMoveProb,rDRGEPReactions,rDRGEPRemovedReactions,shouldDisplay,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,wVector,typeMoves,numKeys,keyAllModels,modelIndex,nonZeroPriorModel,numForwardModelSolves,optimizedParameters,derivativeBased,proposalType,meanVector,valVector)

% optional return arguments

iterAcceptance = 0;
iterSamples = 0;

y1=zeros(N,1);

dParameters=d1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proposal preparation

if RJMCMCType==2 || RJMCMCType==3 || RJMCMCType==4
    
    % Bayesian inference of chemical kinetic models
    
    if algorithm==2 || algorithm==3
        
        if shouldDisplay==1
            disp('sample move')
        end    
        
        % sample new value
        [currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,currentMoveNumber,currentWithinModelMoveNumber,AMoveNum,WMoveNum,SMoveNum,withinModel,AddDeleteReaction,AddDeleteReactionSwap,proposedModelIndex,proposedKeyLocation]=sampleProposalMove(typeSampling,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,AddDelMoveProb,WithinMoveProb,SwapMoveProb,modelWeightScalingRegular,modelWeightScalingSwap,currentModelUpdateWeight,currentModelWeightWithinModel,currentNumTReactions,locationOfCurrent,currentInfluentialReactions,intermediateCurrentInfluentialReactions,reactionRatesDRGReactions,production,consumption,AMoveNum,WMoveNum,SMoveNum,RJMCMCType,rDRGEPReactions,rDRGEPRemovedReactions,currentPossibleMoves,keyAllModels,modelIndex,keyCount,key,currentModelIndex,currentKeyLocation);

       
        if ~any(intermediateProposedInfluentialReactions==1)
          intermediateProposedInfluentialReactions(1)=1;
          intermediateProposedInfluentialReactions=sort(intermediateProposedInfluentialReactions);
        end          
        
              
        if ~any(intermediateProposedInfluentialReactions==2)
          intermediateProposedInfluentialReactions(1)=2;
          intermediateProposedInfluentialReactions=sort(intermediateProposedInfluentialReactions);
        end      
         
        
        if ~any(intermediateProposedInfluentialReactions==9)
          intermediateProposedInfluentialReactions(1)=9;
          intermediateProposedInfluentialReactions=sort(intermediateProposedInfluentialReactions);
        end      
                
        
        if ~any(intermediateProposedInfluentialReactions==4)
          intermediateProposedInfluentialReactions(1)=4;
          intermediateProposedInfluentialReactions=sort(intermediateProposedInfluentialReactions);
        end
        
        
        if ~any(intermediateProposedInfluentialReactions==10)
          intermediateProposedInfluentialReactions(1)=10;
          intermediateProposedInfluentialReactions=sort(intermediateProposedInfluentialReactions);
        end        
 
        if ~any(intermediateProposedInfluentialReactions==12)
          intermediateProposedInfluentialReactions(1)=12;
          intermediateProposedInfluentialReactions=sort(intermediateProposedInfluentialReactions);
        end    
        
        
        %{
         if ~any(intermediateProposedInfluentialReactions==11)
          intermediateProposedInfluentialReactions(1)=11;
          intermediateProposedInfluentialReactions=sort(intermediateProposedInfluentialReactions);
        end        
        %}
        %{
        if ~any(intermediateProposedInfluentialReactions==8)
          intermediateProposedInfluentialReactions(1)=8;
          intermediateProposedInfluentialReactions=sort(intermediateProposedInfluentialReactions);
        end           
        %}
        %}
        if shouldDisplay==1
            disp('rev move prob')
        end        
        
        [reversePossibleMoves,reverseModelUpdateWeight,reverseModelWeightWithinModel,reverseNumTReactions,locationOfProposed,reverseMoveNumber,reverseWithinModelMoveNumber,reactionRatesDRGReactionsRev,productionRev,consumptionRev,rDRGEPReactionsRev,rDRGEPRemovedReactionsRev,key,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,proposedKeyLocation] = reverseMoveProbCalculations(currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,currentMoveNumber,currentWithinModelMoveNumber,withinModel,y1,priorRange,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,modelWeightScalingRegular,RJMCMCType,algorithm,typeSampling,OnTheFlyModelDetermination,AddDelMoveProb,WithinMoveProb,SwapMoveProb,inputX,key,keyAllModels,keyCount,keyMoveDetermination,AllPossibleMoves,AllPossibleNumPaths,AllPossibleNumOfModelsForEachMove,AllPossibleMovesTReactions,AllPossibleMovesNumTReactions,numKeys,nonZeroPriorModel,wVector,K,K1,modelIndex,proposedKeyLocation);
        
        if shouldDisplay==1
            disp('proposalParameters')
        end  
        
        if shouldDisplay==1
          %save currentParameters.mat dParameters
          withinModel
          currentInfluentialReactions
          proposedInfluentialReactions          
        end
        
        [iCIRF,iPIRF,iCIRB,iPIRB,dParametersForward,dParametersBackward] = determineForwardBackwardInfluentialReactions(withinModel,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,dParameters,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
        
        if withinModel==0
            
            % Forward
            % If Add/Delete Move
            [centeringLocationForward,centeringLocationOnBoundaryForward,eigenVectorMatrixForward,eigenValueMatrixForward,gradientVectorForward,HessianMatrixForward,numberOfVariablesForward,AreReactionAddedForward,currentModelIndex,proposedModelIndex]=proposalParameters(iterationCount,RJMCMCType,key,keyCount,iCIRF,iPIRF,N,dParametersForward,dParametersForward,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,keyAllModels,modelIndex,optimizedParameters,derivativeBased,proposalType,currentModelIndex,proposedModelIndex,shouldDisplay,meanVector,valVector);
            
            % Backward
            [centeringLocationBackward,centeringLocationOnBoundaryBackward,eigenVectorMatrixBackward,eigenValueMatrixBackward,gradientVectorBackward,HessianMatrixBackward,numberOfVariablesBackward,AreReactionAddedBackward,currentModelIndex,proposedModelIndex]=proposalParameters(iterationCount,RJMCMCType,key,keyCount,iCIRB,iPIRB,N,dParametersBackward,dParametersBackward,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,keyAllModels,modelIndex,optimizedParameters,derivativeBased,proposalType,currentModelIndex,proposedModelIndex,shouldDisplay,meanVector,valVector);
            
            
        else            
            
            centeringLocationForward=0;
            centeringLocationOnBoundaryForward=0;
            eigenVectorMatrixForward=0;
            eigenValueMatrixForward=0;
            gradientVectorForward=0;
            HessianMatrixForward=0;
            numberOfVariablesForward=0;
            AreReactionAddedForward=0;
            
            centeringLocationBackward=0;
            centeringLocationOnBoundaryBackward=0;
            eigenVectorMatrixBackward=0;
            eigenValueMatrixBackward=0;
            gradientVectorBackward=0;
            HessianMatrixBackward=0;
            numberOfVariablesBackward=0;
            AreReactionAddedBackward=0;            
                                    
        end
        
        if shouldDisplay==1
            disp('sample parameters')
        end  
        
        [y1,detVariance]=sampleProposalContinuous(typeSampling,currentModelIndex,proposedModelIndex,modelProblem,algorithm,nReactions,reactionType,N,priorRange,dParametersForward,iCIRF,iPIRF,currentWithinModelMoveNumber,delta,withinModel,proposalType,centeringLocationForward,centeringLocationOnBoundaryForward,eigenVectorMatrixForward,eigenValueMatrixForward,gradientVectorForward,HessianMatrixForward,numberOfVariablesForward,AreReactionAddedForward,K,K1,optimizedParameters);
                
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


yParameters=y1;

[logN1,logD1,ll1,logp1,logp2,logProp1,logProp2,numForwardModelSolves,AdditionalLnPrior]=calcN1D1(RJMCMCType,key,keyCount,typeSampling,algorithm,N,K,K1,proposalType,centeringLocationForward,centeringLocationOnBoundaryForward,eigenVectorMatrixForward,eigenValueMatrixForward,gradientVectorForward,HessianMatrixForward,numberOfVariablesForward,centeringLocationBackward,centeringLocationOnBoundaryBackward,eigenVectorMatrixBackward,eigenValueMatrixBackward,gradientVectorBackward,HessianMatrixBackward,numberOfVariablesBackward,yParameters,dParameters,dParametersForward,dParametersBackward,ll,logPos,pos,priorRange,inputX,data,delta,dataVariance,currentModelUpdateWeight,reverseModelUpdateWeight,currentModelWeightWithinModel,reverseModelWeightWithinModel,currentNumTReactions,reverseNumTReactions,currentInfluentialReactions,proposedInfluentialReactions,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,iCIRF,iPIRF,iCIRB,iPIRB,currentMoveNumber,reverseMoveNumber,currentWithinModelMoveNumber,reverseWithinModelMoveNumber,AddDeleteReaction,modelProblem,initialIncludedSpeciesPool,targetNode,nReactions,reactionType,nSpecies,speciesReference,stoichiometricMatrix,AddDeleteReactionSwap,shouldDisplay,numForwardModelSolves,optimizedParameters,meanVector,valVector);

[acceptance]=calcAcceptanceRatio(logN1,logD1);


b=rand(1);

if shouldDisplay
    

    currentModelUpdateWeight
    
    reverseModelUpdateWeight
    
    %{
    currentNumTReactions
    
    reverseNumTReactions
    %}
    
    logN1
    
    logD1
    
    
    ll1
    
    ll
    
    logp1
    
    logp2
    
    logProp1
    
    logProp2
    
    logPos
    
    AdditionalLnPrior
        
    centeringLocationOnBoundaryForward
    
    eigenVectorMatrixForward
    
    eigenValueMatrixForward
    
    HessianMatrixForward
    
    gradientVectorForward
    
    numberOfVariablesForward

    centeringLocationForward    
    
    centeringLocationOnBoundaryBackward
    
    eigenVectorMatrixBackward
    
    eigenValueMatrixBackward
    
    HessianMatrixBackward
    
    gradientVectorBackward
    
    numberOfVariablesBackward

    centeringLocationBackward      
            
    detVariance
    
    y1'
    
    d1'
    
    b
        
    withinModel
    
    if b<=acceptance
        display('Accept')
    end
    
    disp('intermediateCurrentInfluentialReactionsForward')
    iCIRF'
    
    disp('intermediateProposedInfluentialReactionsForward')
    iPIRF'
    
    disp('intermediateCurrentInfluentialReactionsBackward')
    iCIRB'
    
    disp('intermediateProposedInfluentialReactionsBackward')
    iPIRB'
    
    disp('currentInfluentialReactions')    
    currentInfluentialReactions'
 
    disp('proposedInfluentialReactions')    
    proposedInfluentialReactions'    

    acceptance       
    
    display('End of iteration')     
    
    pause
end

% Accept/reject proposed sample

moveNumber=length(setxor(currentInfluentialReactions,proposedInfluentialReactions));

if moveNumber~=0
    typeMoves(moveNumber,1)=typeMoves(moveNumber,1)+1;
end

if (b<=acceptance)
    
    iterAcceptance = iterAcceptance+1;
    iterSamples = iterSamples+1;
    
    samples(1:N,1) = y1;
        
    % Update the previous-loglikelihood value
    currentInfluentialReactions=proposedInfluentialReactions;
    intermediateCurrentInfluentialReactions=intermediateProposedInfluentialReactions;
    
    currentModelUpdateWeight=reverseModelUpdateWeight;
    currentModelWeightWithinModel=reverseModelWeightWithinModel;
    currentNumTReactions=reverseNumTReactions;
    
    currentPossibleMoves=reversePossibleMoves;
    
    currentModelIndex=proposedModelIndex;
    currentKeyLocation=proposedKeyLocation;
    
    locationOfCurrent=locationOfProposed;
    
    reactionRatesDRGReactions=reactionRatesDRGReactionsRev;
    production=productionRev;
    consumption=consumptionRev;
    rDRGEPReactions=rDRGEPReactionsRev;
    rDRGEPRemovedReactions=rDRGEPRemovedReactionsRev;
    
    ll = ll1;    
    logPos=logp1;
    logPosTrue=logp1+AdditionalLnPrior;
    
    % update the acceptance counts
    if withinModel==0
        AMoveAcc=AMoveAcc+1;
    elseif withinModel==1
        WMoveAcc=WMoveAcc+1;
    else
        SMoveAcc=SMoveAcc+1;
    end
    
    if moveNumber~=0
        typeMoves(moveNumber,2)=typeMoves(moveNumber,2)+1;
    end

else
    
    samples(1:N,1) = d1;
    iterSamples = iterSamples+1;
    
end
