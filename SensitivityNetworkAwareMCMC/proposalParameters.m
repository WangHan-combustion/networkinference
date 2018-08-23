function [centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables,AreReactionAdded,currentModelIndex,proposedModelIndex]=proposalParameters(iterationCount,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,keyAllModels,modelIndex,optimizedParameters,derivativeBased,proposalType,currentModelIndex,proposedModelIndex,shouldDisplay,meanVector,valVector)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the centering location, the eigenVectorMatrix, the
% eigenValueMatrix, the gradientVector, the HessianMatrix at the centering location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize eigenVectorMatrix, eigenValueMatrix, gradientVector,
% HessianMatrix
eigenVectorMatrix=0;
eigenValueMatrix=0;
eigenValueVector=0;
gradientVector=0;
HessianMatrix=0;
offSetRatio=0.01;
selfCode=1;   

centeringLocation=[];
centeringLocationOnBoundary=[];
addReactions=[];
deleteReactions=[];
numberOfVariables=0;
AreReactionAdded=0;

hessianThreshold=1E-1;     % hessian lower bound
hessianUpperBound=1E12;   % hessian uppper bound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shouldDisplay==1
   
    disp('determining the reactions')
    
end

% Determine the reactions to be added or deleted

for n=1:nReactions
    
    if any(n==intermediateProposedInfluentialReactions) && any(n==intermediateCurrentInfluentialReactions)
    elseif any(n==intermediateProposedInfluentialReactions)
        
        AreReactionAdded=1;
        
        addReactions=[addReactions n];
        
        if reactionType(n)==1
            
            numberOfVariables=numberOfVariables+2;
            
        else
            
            numberOfVariables=numberOfVariables+1;
            
        end        
        
    elseif any(n==intermediateCurrentInfluentialReactions)
        
        deleteReactions=[deleteReactions n];
        
        if reactionType(n)==1
            
            numberOfVariables=numberOfVariables+2;
            
        else
            
            numberOfVariables=numberOfVariables+1;
            
        end
        
    else
        
    end
    
   
end

% Check whether an optimization problem needs to be solved

% Initialize lower bounds, upper bounds and nominal center values
lb=zeros(1,numberOfVariables);
ub=zeros(1,numberOfVariables);
x0=zeros(1,numberOfVariables);

%%% If reactions are being added:
if isempty(deleteReactions) && numberOfVariables~=0
    
    variableNumber=0;
    
    if shouldDisplay==1
      disp('initializing parameters')
    end

   
    for i=1:length(addReactions)
        
        if reactionType(addReactions(i))==1
            
            % Initial location: mid-point of the prior
            a1=rand(1);
            a2=rand(1);
            
            lb(variableNumber+1)=priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))));
            lb(variableNumber+2)=priorRange(1,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1))));
            
            ub(variableNumber+1)=priorRange(2,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))));
            ub(variableNumber+2)=priorRange(2,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1))));
            
            %x0(variableNumber+1)=priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+a1*(-priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+priorRange(2,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
            %x0(variableNumber+2)=priorRange(1,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+a2*(-priorRange(1,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+priorRange(2,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
            
            x0(variableNumber+1)=meanVector(addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+sqrt(valVector(addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))))*randn(1);
            x0(variableNumber+2)=meanVector(1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+sqrt(valVector(1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))))*randn(1);
            
            variableNumber=variableNumber+2;
            
        else
            
            % Initial location: mid-point of the prior
            a1=rand(1);
            
            lb(variableNumber+1)=priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))));
            
            ub(variableNumber+1)=priorRange(2,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))));

            
            x0(variableNumber+1)=meanVector(addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+sqrt(valVector(addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))))*randn(1);                  
            %x0(variableNumber+1)=priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+a1*(-priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+priorRange(2,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        
            variableNumber=variableNumber+1;
            
        end
        
    end
    
    if shouldDisplay==1
   
       disp('adding reactions') 
        
       disp('optimizing parameters')
    
       %save x0.mat x0;
       
       %x0
       
       %load x0.mat
       
    end 
   
    if optimizedParameters==1
        
        %if modelProblem==8 || modelProblem==9 || modelProblem==10
        % Proposal type 1 is a normal distribution
        if proposalType==1
                                            
                % Model problems 8, 9, 10 use derivative-based methods by
                % default
                warning('off','all')                
                  
                if shouldDisplay==1
                 opt = optimset('Display','on','TolFun',1E-6,'TolX',1E-6);                    
                else
                 opt = optimset('Display','off','TolFun',1E-6,'TolX',1E-6);                    
                end    
                                
                [centeringLocation,~,~,~,~,HessianMatrix]=fminunc(@(centeringLocation)optimizationObjective(centeringLocation,addReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector),x0,opt);
 
                
                centeringLocationOnBoundary=0;
        else
            
            if derivativeBased==1
            
                if shouldDisplay==1

                  opt = optimset('Algorithm','active-set','FunValCheck','on','Display','iter-detailed','TolFun',1E-6,'TolX',1E-6);
                                       
                else
                    
                 opt = optimset('Algorithm','active-set','Display','off','TolFun',1E-6,'TolX',1E-6);
                 
                end
                                 
                % Compute forward model
                
                %{
                addReactions
                coordinate1=1;
                for x5=3:0.1:5
                    
                   coordinate2=1;
                   for x7=1:0.1:3 
                      [f(coordinate1,coordinate2)]=optimizationObjective([x5 x7],addReactions,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay);
                      coordinate2=coordinate2+1;
                   end    
                   coordinate1=coordinate1+1
                end    
                
                save forwardModel.mat f;
                %}
                
                [centeringLocation,fval,~,~,lambda,grad,H]=fmincon(@(centeringLocation)optimizationObjective(centeringLocation,addReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector),x0,[],[],[],[],lb,ub,[],opt);

                
                
            else
                
                
                opt.algorithm = NLOPT_LN_COBYLA;
                opt.lower_bounds = lb;
                opt.upper_bounds = ub;
                opt.min_objective = (@(centeringLocation)optimizationObjective(centeringLocation,addReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector));
                
                
                opt.fc_tol = 1e-5*ones(1,length(lb));
                opt.xtol_rel = 1e-6;
                opt.ftol_rel = 1e-6;
                %opt.verbose=1;
                
                [centeringLocation, fval, retcode] = nlopt_optimize(opt, x0);
                
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check if the centering location is the boundary
        
        %if modelProblem~=8 && modelProblem~=9 && modelProblem~=10
        if proposalType==3
            
            l=length(addReactions);
            
            variableNumber=0;
            
            for i=1:l
                
                if reactionType(addReactions(i))==1
                    
                    if centeringLocation(variableNumber+1)<=(priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocation(variableNumber+1)=priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+...
                            offSetRatio*(priorRange(2,addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+1)=1;
                    elseif centeringLocation(variableNumber+1)>=(priorRange(2,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocation(variableNumber+1)=priorRange(2,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))-...
                            offSetRatio*(priorRange(2,addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+1)=2;
                    else
                        centeringLocationOnBoundary(variableNumber+1)=0;
                    end                    
 
                    if centeringLocation(variableNumber+2)<=(priorRange(1,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocation(variableNumber+2)=priorRange(1,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+...
                            offSetRatio*(priorRange(2,1+addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+2)=1;
                    elseif centeringLocation(variableNumber+2)>=(priorRange(2,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocation(variableNumber+2)=priorRange(2,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))-...
                            offSetRatio*(priorRange(2,1+addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,1+addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+2)=2;
                    else
                        centeringLocationOnBoundary(variableNumber+2)=0;
                    end                     
                    
                    variableNumber=variableNumber+2;
                else

                    if centeringLocation(variableNumber+1)<=(priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocation(variableNumber+1)=priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))+...
                            offSetRatio*(priorRange(2,addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+1)=1;
                    elseif centeringLocation(variableNumber+1)>=(priorRange(2,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocation(variableNumber+1)=priorRange(2,addReactions(i)+sum(reactionType(1:(addReactions(i)-1))))-...
                            offSetRatio*(priorRange(2,addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,addReactions(i)+sum(reactionType(1:(addReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+1)=2;
                    else
                        centeringLocationOnBoundary(variableNumber+1)=0;
                    end                    
                                        
                    variableNumber=variableNumber+1;
                end
                
            end

            
        end   
                
        if proposalType==1
            
             if shouldDisplay==1
   
               disp('determining finite differences')
    
             end
                                    
            if selfCode==1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % My code for finite-difference gradient Hessian
                
                if shouldDisplay==1
                    
                    disp('determining the gradients')
                    
                end
                                
                if shouldDisplay==1
                    
                    centeringLocation
             
                    disp('determining the hessian')
                    
                end
                
                [HessianMatrix]=hessianCalculator(centeringLocation,addReactions,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
                
                gradientVector=-gradientVector;
                
            else
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DERIVEST gradient and Hessian
                
                [hessianDERIVEST,err] = hessian(@(x)optimizationObjective(x,addReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector),centeringLocation);
                                
                [gradientDERIVEST,err] = gradest(@(x)optimizationObjective(x,addReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector),centeringLocation);
                
                % Convert the gradient to the gradient of the true objective
                gradientVector=-gradientDERIVEST;
                
                % Hessian matrix should be positive definite because the
                % objective function is -logp or -ll which is minimized
                HessianMatrix=hessianDERIVEST;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
        end
        
        if shouldDisplay==1
            disp('Initial Hessian')
            x0
            centeringLocation
            HessianMatrix 
        end    
        
        
        [eigenVectorMatrix,eigenValueMatrix]=eig(HessianMatrix);
                
        % Return the main diagonal
        eigenValueVector=diag(eigenValueMatrix);       
                        
        %if modelProblem==8 || modelProblem==9 || modelProblem==10
        if proposalType==1
                                                                     
                indexLocations=find(eigenValueVector<=hessianThreshold);
                
                % Truncate eigenvalues below lower bound
                if ~isempty(indexLocations)
                    eigenValueVector(indexLocations)=hessianThreshold;
                end
                
                indexLocations=find(eigenValueVector>=hessianUpperBound);
                
                % Truncate eigenvalues above upper bound
                if ~isempty(indexLocations)
                    eigenValueVector(indexLocations)=hessianUpperBound;
                end
                
        else
            
                                                       
               for variableNumber=1:numberOfVariables
                  
                   if centeringLocationOnBoundary(variableNumber)==1
                       if gradientVector(variableNumber)>0
                          gradientVector(variableNumber)=-1E-6; 
                       end
                       if HessianMatrix(variableNumber,variableNumber)<1E-6
                          HessianMatrix(variableNumber,variableNumber)=1E-6; 
                       end   
                       if HessianMatrix(variableNumber,variableNumber)>1E6
                          HessianMatrix(variableNumber,variableNumber)=1E6; 
                       end                        
                   elseif centeringLocationOnBoundary(variableNumber)==2
                       if gradientVector(variableNumber)<0
                           gradientVector(variableNumber)=1E-6;
                       end
                       if HessianMatrix(variableNumber,variableNumber)<1E-6
                          HessianMatrix(variableNumber,variableNumber)=1E-6; 
                       end  
                       if HessianMatrix(variableNumber,variableNumber)>1E6
                          HessianMatrix(variableNumber,variableNumber)=1E6; 
                       end                        
                   elseif centeringLocationOnBoundary(variableNumber)==0
                       if HessianMatrix(variableNumber,variableNumber)<hessianThreshold
                          HessianMatrix(variableNumber,variableNumber)=hessianThreshold;
                       end                          
                       if HessianMatrix(variableNumber,variableNumber)>1E6
                          HessianMatrix(variableNumber,variableNumber)=1E6; 
                       end                                
                   end    
                   
               end    
        end
        
                
    else
        
        eigenValueVector=zeros(numberOfVariables,1);
        
        for i=1:length(addReactions)
            
            if reactionType(addReactions(i))==1
                
                eigenValueVector(i+sum(reactionType(addReactions(1:i-1))))=((priorRange(2,addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,addReactions(i)+sum(reactionType(1:addReactions(i)-1))))/1)^(-2);
                eigenValueVector(1+i+sum(reactionType(addReactions(1:i-1))))=((priorRange(2,1+addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,1+addReactions(i)+sum(reactionType(1:addReactions(i)-1))))/1)^(-2);
                
            else
                
                eigenValueVector(i+sum(reactionType(addReactions(1:i-1))))=((priorRange(2,addReactions(i)+sum(reactionType(1:addReactions(i)-1)))-priorRange(1,addReactions(i)+sum(reactionType(1:addReactions(i)-1))))/1)^(-2);
                
            end
            
        end

        eigenVectorMatrix=eye(numberOfVariables);   
        centeringLocation=x0;
        centeringLocationOnBoundary=zeros(1,numberOfVariables);
        gradientVector=ones(numberOfVariables,1);
                
    end

% If reactions are being deleted    
elseif isempty(addReactions) && numberOfVariables~=0
    
    variableNumber=0;
    
    if shouldDisplay==1
      disp('initializing parameters')    
    end
    
    for i=1:length(deleteReactions)
        
        if reactionType(deleteReactions(i))==1
            
            % Initial point is the mid-point of the prior range
            a1=rand(1);
            a2=rand(1);
            
            lb(variableNumber+1)=priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))));
            lb(variableNumber+2)=priorRange(1,1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))));
            
            ub(variableNumber+1)=priorRange(2,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))));
            ub(variableNumber+2)=priorRange(2,1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))));
            
            %x0(variableNumber+1)=priorRange(1,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))+a1*(-priorRange(1,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))+priorRange(2,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1))));
            %x0(variableNumber+2)=priorRange(1,1+deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))+a2*(-priorRange(1,1+deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))+priorRange(2,1+deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1))));

            x0(variableNumber+1)=meanVector(deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))+sqrt(valVector(deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))))*randn(1);
            x0(variableNumber+2)=meanVector(1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))+sqrt(valVector(1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))))*randn(1);
                        
            variableNumber=variableNumber+2;
            
        else
            
            % Initial point is the mid-point of the prior range
            a1=rand(1);
            
            lb(variableNumber+1)=priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))));
            
            ub(variableNumber+1)=priorRange(2,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))));

            x0(variableNumber+1)=meanVector(deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))+sqrt(valVector(deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))))*randn(1);
            
            %x0(variableNumber+1)=priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))+a1*(-priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))+priorRange(2,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
            
            variableNumber=variableNumber+1;
            
        end
        
    end    

    if shouldDisplay==1
   
       disp('deleting reactions') 
        
       disp('optimizing parameters')
    
       %save x0.mat x0
       
       %x0
       
       %load x0.mat
       
    end    
    
   
    if optimizedParameters==1
        
        %if modelProblem==8 || modelProblem==9 || modelProblem==10
        if proposalType==1
            
                warning('off','all')
                
                opt = optimset('Display','off','TolFun',1E-6,'TolX',1E-6);
                
                [centeringLocation,~,~,~,~,HessianMatrix]=fminunc(@(centeringLocation)optimizationObjective(centeringLocation,deleteReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector),x0,opt);

                
                
                centeringLocationOnBoundary=0;                
                
        else
            
            if derivativeBased==1
               
                if shouldDisplay==1

                  opt = optimset('Algorithm','active-set','FunValCheck','on','Display','iter-detailed','TolFun',1E-6,'TolX',1E-6);
                                       
                else
                    
                 opt = optimset('Algorithm','active-set','Display','off','TolFun',1E-6,'TolX',1E-6);
                 
                end
                            
                
                [centeringLocation,fval,~,~,lambda,grad,H]=fmincon(@(centeringLocation)optimizationObjective(centeringLocation,deleteReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector),x0,[],[],[],[],lb,ub,[],opt);

                
                
            else
                
                
                opt.algorithm = NLOPT_LN_COBYLA;
                opt.lower_bounds = lb;
                opt.upper_bounds = ub;
                opt.min_objective = (@(centeringLocation)optimizationObjective(centeringLocation,deleteReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector));
                opt.fc_tol = 1e-5*ones(1,length(lb));
                opt.xtol_rel = 1e-6;
                opt.ftol_rel = 1e-6;
                %opt.verbose=1;
                
                [centeringLocation, fval, retcode] = nlopt_optimize(opt, x0);
                
                                
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check if the centering location is the boundary

        %if modelProblem~=8 && modelProblem~=9 && modelProblem~=10
        if proposalType==3
            
            l=length(deleteReactions);
            
            variableNumber=0;
            
            for i=1:l
                
                if reactionType(deleteReactions(i))==1
                    
                    if centeringLocation(variableNumber+1)<=(priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocation(variableNumber+1)=priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))+...
                            offSetRatio*(priorRange(2,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+1)=1;
                    elseif centeringLocation(variableNumber+1)>=(priorRange(2,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocation(variableNumber+1)=priorRange(2,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))-...
                            offSetRatio*(priorRange(2,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+1)=2;
                    else
                        centeringLocationOnBoundary(variableNumber+1)=0;
                    end                    
 
                    if centeringLocation(variableNumber+2)<=(priorRange(1,1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocation(variableNumber+2)=priorRange(1,1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))+...
                            offSetRatio*(priorRange(2,1+deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+2)=1;
                    elseif centeringLocation(variableNumber+2)>=(priorRange(2,1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocation(variableNumber+2)=priorRange(2,1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))-...
                            offSetRatio*(priorRange(2,1+deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,1+deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+2)=2;
                    else
                        centeringLocationOnBoundary(variableNumber+2)=0;
                    end                     
                    
                    variableNumber=variableNumber+2;
                else

                    if centeringLocation(variableNumber+1)<=(priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocation(variableNumber+1)=priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))+...
                            offSetRatio*(priorRange(2,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+1)=1;
                    elseif centeringLocation(variableNumber+1)>=(priorRange(2,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocation(variableNumber+1)=priorRange(2,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1))))-...
                            offSetRatio*(priorRange(2,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,deleteReactions(i)+sum(reactionType(1:(deleteReactions(i)-1)))));
                        centeringLocationOnBoundary(variableNumber+1)=2;
                    else
                        centeringLocationOnBoundary(variableNumber+1)=0;
                    end                    
                                        
                    variableNumber=variableNumber+1;
                end
                
            end

        end   
        
         if proposalType==1
             
             if shouldDisplay==1
                 
                 disp('determining finite differences')
                 
             end
             
             
             if selfCode==1
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % my finite difference code for gradient and Hessian
                 
                 if shouldDisplay==1
                     
                     disp('determining the gradients')
                     
                 end
                 
                 
                 %[gradientVector]=gradientCalculator(centeringLocation,deleteReactions,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType);
                 
                 if shouldDisplay==1
                     
                     disp('determining the hessian')
                     
                 end
                 
                 
                 [HessianMatrix]=hessianCalculator(centeringLocation,deleteReactions,RJMCMCType,key,keyCount,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,N,valArray,val,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,initialIncludedSpeciesPool,targetNode,reactionType,speciesReference,stoichiometricMatrix,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector);
                 
                 gradientVector=-gradientVector;
                 
             else
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % DERIVEST code for gradient and Hessian
                 
                 [hessianDERIVEST,err] = hessian(@(x)optimizationObjective(x,deleteReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector),centeringLocation);
                 
                 [gradientDERIVEST,err] = gradest(@(x)optimizationObjective(x,deleteReactions,key,keyCount,intermediateProposedInfluentialReactions,N,valArray,leadLogLikelihood,logPos,pos,priorRange,input,data,delta,dataVariance,nReactions,nSpecies,modelProblem,reactionType,typeSampling,numForwardModelSolves,shouldDisplay,proposalType,meanVector,valVector),centeringLocation);
                 
                 
                 % The true gradient is the negative of the computed gradient
                 gradientVector=-gradientDERIVEST;
                 HessianMatrix=hessianDERIVEST;
                 
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             end
         end
       
        if shouldDisplay==1
           disp('Initial Hessian') 
           HessianMatrix 
        end    
         
        [eigenVectorMatrix,eigenValueMatrix]=eig(HessianMatrix);
        
        % Return the main diagonal
        eigenValueVector=diag(eigenValueMatrix);         
        
                
        %if modelProblem==8 || modelProblem==9 || modelProblem==10
        if proposalType==1
          
                        
            indexLocations=find(eigenValueVector<=hessianThreshold);
            
            % Truncate eigenvalues below zero
            if ~isempty(indexLocations)
                eigenValueVector(indexLocations)=hessianThreshold;
            end
            
            indexLocations=find(eigenValueVector>=hessianUpperBound);
            
            % Truncate eigenvalues below zero
            if ~isempty(indexLocations)
                eigenValueVector(indexLocations)=hessianUpperBound;
            end            
        
        else    
                          
               for variableNumber=1:numberOfVariables
                  
                   if centeringLocationOnBoundary(variableNumber)==1
                       if gradientVector(variableNumber)>0
                          gradientVector(variableNumber)=-1E-6; 
                       end
                       if HessianMatrix(variableNumber,variableNumber)<1E-6
                          HessianMatrix(variableNumber,variableNumber)=1E-6; 
                       end        
                       if HessianMatrix(variableNumber,variableNumber)>1E6
                          HessianMatrix(variableNumber,variableNumber)=1E6; 
                       end                                                                              
                   elseif centeringLocationOnBoundary(variableNumber)==2
                       if gradientVector(variableNumber)<0
                           gradientVector(variableNumber)=1E-6;
                       end
                       if HessianMatrix(variableNumber,variableNumber)<1E-6
                          HessianMatrix(variableNumber,variableNumber)=1E-6; 
                       end    
                       if HessianMatrix(variableNumber,variableNumber)>1E6
                          HessianMatrix(variableNumber,variableNumber)=1E6; 
                       end                                
                   elseif centeringLocationOnBoundary(variableNumber)==0
                       if HessianMatrix(variableNumber,variableNumber)<hessianThreshold
                          HessianMatrix(variableNumber,variableNumber)=hessianThreshold;
                       end   
                        if HessianMatrix(variableNumber,variableNumber)>1E6
                          HessianMatrix(variableNumber,variableNumber)=1E6; 
                       end                               
                   end    
                   
               end             
                                
        end    
            
    else
        
        eigenValueVector=zeros(1,numberOfVariables);
        
        for i=1:length(deleteReactions)
            
            if reactionType(deleteReactions(i))==1
                                
                eigenValueVector(i+sum(reactionType(deleteReactions(1:i-1))))=((priorRange(2,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1))))/1)^(-2);
                eigenValueVector(1+i+sum(reactionType(deleteReactions(1:i-1))))=((priorRange(2,1+deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,1+deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1))))/1)^(-2);
                
                
            else
                                
                eigenValueVector(i+sum(reactionType(deleteReactions(1:i-1))))=((priorRange(2,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1)))-priorRange(1,deleteReactions(i)+sum(reactionType(1:deleteReactions(i)-1))))/1)^(-2);
                
            end
            
        end        
    
        eigenVectorMatrix=eye(numberOfVariables);   
        centeringLocation=x0;
        centeringLocationOnBoundary=zeros(1,numberOfVariables);
        gradientVector=ones(numberOfVariables,1);
                
    end    
        
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if shouldDisplay==1
 disp('end of proposal parameters')
end    

if optimizedParameters==1 
    
    %if modelProblem==8 || modelProblem==9 || modelProblem==10
    if proposalType==1
        
        % Diagonalize the eigenvalues
        eigenValueMatrix=diag(eigenValueVector);
        HessianMatrix=eigenVectorMatrix*eigenValueMatrix*eigenVectorMatrix';
        
    end
    
else
    
        eigenValueMatrix=diag(eigenValueVector);    
        HessianMatrix=eigenVectorMatrix*eigenValueMatrix*eigenVectorMatrix';
        
end



end
