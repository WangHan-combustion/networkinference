function [formatSamples,leadLogLikelihood,leadLogPosterior,logPosTrue,numForwardModelSolves,logPosterior] = initializeChain(typeSampling,inputX,data,modelProblem,N,numIteration,pos,priorRange,numForwardModelSolves,dataVariance,nReactions,reactionType,nSpecies,currentInfluentialReactions,intermediateCurrentInfluentialReactions,delta,shouldDisplay,key,keyCount,meanVector,valVector)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remember the parameter values here are after taking logarithm, so you
% don't have to initialize with log(parameter value)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logPosterior=zeros(1,numIteration);

    %% 13 reaction (model problem 1)
    if modelProblem==1
        formatSamples=delta*ones(N,numIteration);
        
        if any(1==intermediateCurrentInfluentialReactions)
         formatSamples(1,numIteration)=(priorRange(2,1)+priorRange(1,1))/2+(-1+2*rand(1))*(priorRange(2,1)-priorRange(1,1))/4;
        else
         formatSamples(1,numIteration)=delta;            
        end

        if any(2==intermediateCurrentInfluentialReactions)
         formatSamples(2,numIteration)=(priorRange(2,2)+priorRange(1,2))/2+(-1+2*rand(1))*(priorRange(2,2)-priorRange(1,2))/4;
         formatSamples(3,numIteration)=(priorRange(2,3)+priorRange(1,3))/2+(-1+2*rand(1))*(priorRange(2,3)-priorRange(1,3))/4;
        else
         formatSamples(2,numIteration)=delta;
         formatSamples(3,numIteration)=delta;            
        end
        
        if any(3==intermediateCurrentInfluentialReactions)
         formatSamples(4,numIteration)=(priorRange(2,4)+priorRange(1,4))/2+(-1+2*rand(1))*(priorRange(2,4)-priorRange(1,4))/4;
        else
         formatSamples(4,numIteration)=delta;              
        end

        if any(4==intermediateCurrentInfluentialReactions)
         formatSamples(5,numIteration)=(priorRange(2,5)+priorRange(1,5))/2+(-1+2*rand(1))*(priorRange(2,5)-priorRange(1,5))/4;
        else
         formatSamples(5,numIteration)=delta;               
        end
        
        if any(5==intermediateCurrentInfluentialReactions)  
         formatSamples(6,numIteration)=(priorRange(2,6)+priorRange(1,6))/2+(-1+2*rand(1))*(priorRange(2,6)-priorRange(1,6))/4;
        else
         formatSamples(6,numIteration)=delta;                 
        end

        if any(6==intermediateCurrentInfluentialReactions)    
         formatSamples(7,numIteration)=(priorRange(2,7)+priorRange(1,7))/2+(-1+2*rand(1))*(priorRange(2,7)-priorRange(1,7))/4;
        else
         formatSamples(7,numIteration)=delta;               
        end        

        if any(7==intermediateCurrentInfluentialReactions)
         formatSamples(8,numIteration)=(priorRange(2,8)+priorRange(1,8))/2+(-1+2*rand(1))*(priorRange(2,8)-priorRange(1,8))/4;
        else
         formatSamples(8,numIteration)=delta;              
        end
        
        if any(8==intermediateCurrentInfluentialReactions)  
         formatSamples(9,numIteration)=(priorRange(2,9)+priorRange(1,9))/2+(-1+2*rand(1))*(priorRange(2,9)-priorRange(1,9))/4;
        else
         formatSamples(9,numIteration)=delta;             
        end
        
        if any(9==intermediateCurrentInfluentialReactions)   
         %formatSamples(10,numIteration)=(priorRange(2,10)+priorRange(1,10))/2+(-1+2*rand(1))*(priorRange(2,10)-priorRange(1,10))/4;
         formatSamples(10,numIteration)=1.0;
        else
         formatSamples(10,numIteration)=delta;             
        end
        
        if any(10==intermediateCurrentInfluentialReactions)   
         formatSamples(11,numIteration)=(priorRange(2,11)+priorRange(1,11))/2+(-1+2*rand(1))*(priorRange(2,11)-priorRange(1,11))/4;
        else
         formatSamples(11,numIteration)=delta;              
        end
        
        if any(11==intermediateCurrentInfluentialReactions)   
         formatSamples(12,numIteration)=(priorRange(2,12)+priorRange(1,12))/2+(-1+2*rand(1))*(priorRange(2,12)-priorRange(1,12))/4;
        else
         formatSamples(12,numIteration)=delta;            
        end
        
        if any(12==intermediateCurrentInfluentialReactions)    
         formatSamples(13,numIteration)=(priorRange(2,13)+priorRange(1,13))/2+(-1+2*rand(1))*(priorRange(2,13)-priorRange(1,13))/4;
        else
         formatSamples(13,numIteration)=delta;            
        end
        
        if any(13==intermediateCurrentInfluentialReactions) 
         formatSamples(14,numIteration)=(priorRange(2,14)+priorRange(1,14))/2+(-1+2*rand(1))*(priorRange(2,14)-priorRange(1,14))/4;
        else
         formatSamples(14,numIteration)=delta;             
        end
        
        %{
        formatSamples(1:14,numIteration)=[-0.0948
                                           1.4144
                                           0.3948
                                           0.9309
                                        -100.0000
                                           2.3588
                                           0.7419
                                        -100.0000
                                           0.0964
                                        -100.0000
                                           0.1620
                                           0.0033
                                           4.0987
                                           2.0295];   
        %}
        
        %{
        formatSamples(1:14,numIteration)=  [ -0.3059
                                              1.5146
                                              0.1880
                                           -100.0000
                                           -100.0000
                                           -100.0000
                                           -100.0000
                                           -100.0000
                                              0.4785
                                           -100.0000
                                           -100.0000
                                              0.5510
                                           -100.0000
                                              2.4303];
        
        
        load ff.mat
        
        formatSamples(1:14,numIteration)=ff;
        %}
        
        
        [leadLogLikelihood,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,inputX,data,formatSamples(1:N,numIteration),pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);

        lprior=logprior(N,formatSamples(1:N,numIteration),priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateCurrentInfluentialReactions,meanVector,valVector);        
        
        leadLogPosterior=leadLogLikelihood+lprior;  
        
        if typeSampling==1
        
          if ~all(sort(currentInfluentialReactions)==intermediateCurrentInfluentialReactions)
           
             diffReactions=setdiff(currentInfluentialReactions,intermediateCurrentInfluentialReactions); 
              
             logPosTrue=leadLogPosterior+logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions,meanVector,valVector);

          else
            
             logPosTrue=leadLogPosterior;
            
          end    
        
        else
          
            logPosTrue=leadLogPosterior;
            
        end    
        
        
    end
    %%

     %% 13 reaction (model problem 11)
    if modelProblem==11
        formatSamples=delta*ones(N,numIteration);                     
                
        if any(1==intermediateCurrentInfluentialReactions)
         %formatSamples(1,numIteration)=meanVector(1)+sqrt(valVector(1))*randn(1);
         formatSamples(1,numIteration)=0.0;
        else
         formatSamples(1,numIteration)=delta;            
        end

        if any(2==intermediateCurrentInfluentialReactions)  
         %formatSamples(2,numIteration)=meanVector(2)+sqrt(valVector(2))*randn(1);
         %formatSamples(3,numIteration)=meanVector(3)+sqrt(valVector(3))*randn(1);
         formatSamples(2,numIteration)=1.5;
         formatSamples(3,numIteration)=0.0;
        else
         formatSamples(2,numIteration)=delta;
         formatSamples(3,numIteration)=delta;            
        end
        
        if any(3==intermediateCurrentInfluentialReactions)
         formatSamples(4,numIteration)=meanVector(4)+sqrt(valVector(4))*randn(1);
        else
         formatSamples(4,numIteration)=delta;              
        end

        if any(4==intermediateCurrentInfluentialReactions)
          formatSamples(5,numIteration)=2.0;
          %formatSamples(5,numIteration)=meanVector(5)+sqrt(valVector(5))*randn(1);
        else
         formatSamples(5,numIteration)=delta;               
        end
        
        if any(5==intermediateCurrentInfluentialReactions)  
         formatSamples(6,numIteration)=meanVector(6)+sqrt(valVector(6))*randn(1);
        else
         formatSamples(6,numIteration)=delta;                 
        end

        if any(6==intermediateCurrentInfluentialReactions)    
         formatSamples(7,numIteration)=meanVector(7)+sqrt(valVector(7))*randn(1);
        else
         formatSamples(7,numIteration)=delta;               
        end        

        if any(7==intermediateCurrentInfluentialReactions)
         formatSamples(8,numIteration)=meanVector(8)+sqrt(valVector(8))*randn(1);
        else
         formatSamples(8,numIteration)=delta;              
        end
        
        if any(8==intermediateCurrentInfluentialReactions)  
         %formatSamples(9,numIteration)=0.4953;
         formatSamples(9,numIteration)=meanVector(9)+sqrt(valVector(9))*randn(1);
        else
         formatSamples(9,numIteration)=delta;             
        end
        
        if any(9==intermediateCurrentInfluentialReactions)   
         %formatSamples(10,numIteration)=meanVector(10)+sqrt(valVector(10))*randn(1);
         formatSamples(10,numIteration)=1.0;
        else
         formatSamples(10,numIteration)=delta;             
        end
        
        if any(10==intermediateCurrentInfluentialReactions)   
         %formatSamples(11,numIteration)=meanVector(11)+sqrt(valVector(11))*randn(1);
         formatSamples(11,numIteration)=0.0;        
        else
         formatSamples(11,numIteration)=delta;              
        end
        
        if any(11==intermediateCurrentInfluentialReactions)   
         formatSamples(12,numIteration)=meanVector(12)+sqrt(valVector(12))*randn(1);
         %formatSamples(12,numIteration)=0.5175;         
         
        else
         formatSamples(12,numIteration)=delta;            
        end
        
        if any(12==intermediateCurrentInfluentialReactions)    
         %formatSamples(13,numIteration)=meanVector(13)+sqrt(valVector(13))*randn(1);
         formatSamples(13,numIteration)=4.0; 
        else
         formatSamples(13,numIteration)=delta;            
        end
        
        if any(13==intermediateCurrentInfluentialReactions) 
         formatSamples(14,numIteration)=meanVector(14)+sqrt(valVector(14))*randn(1);
         %formatSamples(14,numIteration)=2.4422;        
        else
         formatSamples(14,numIteration)=delta;             
        end
                
        [leadLogLikelihood,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,inputX,data,formatSamples(1:N,numIteration),pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);

        lprior=logprior(N,formatSamples(1:N,numIteration),priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateCurrentInfluentialReactions,meanVector,valVector);        
        
        leadLogPosterior=leadLogLikelihood+lprior;               
       
        if typeSampling==1
            
            if ~all(sort(currentInfluentialReactions)==intermediateCurrentInfluentialReactions)
                
                diffReactions=setdiff(currentInfluentialReactions,intermediateCurrentInfluentialReactions); 
                
                logPosTrue=leadLogPosterior+logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions,meanVector,valVector);
                
            else
                
                logPosTrue=leadLogPosterior;
                
            end
            
        else
            
            logPosTrue=leadLogPosterior;
            
        end
        
        
    end
    %%   
    
 
    %%    
       
    %% 10 reaction (model problem 4)
    if modelProblem==4
        
        formatSamples=delta*ones(N,numIteration);
        
        if any(1==intermediateCurrentInfluentialReactions)
         formatSamples(1,numIteration)=priorRange(1,1)+rand(1)*(priorRange(2,1)-priorRange(1,1));
        else
         formatSamples(1,numIteration)=delta;            
        end

        if any(2==intermediateCurrentInfluentialReactions)
         formatSamples(2,numIteration)=priorRange(1,2)+rand(1)*(priorRange(2,2)-priorRange(1,2));
        else
         formatSamples(2,numIteration)=delta;     
        end
        
        if any(3==intermediateCurrentInfluentialReactions)  
         formatSamples(3,numIteration)=priorRange(1,3)+rand(1)*(priorRange(2,3)-priorRange(1,3));
        else
         formatSamples(3,numIteration)=delta;                
        end

        if any(4==intermediateCurrentInfluentialReactions)   
         formatSamples(4,numIteration)=priorRange(1,4)+rand(1)*(priorRange(2,4)-priorRange(1,4));
        else
         formatSamples(4,numIteration)=delta;              
        end
        
        if any(5==intermediateCurrentInfluentialReactions)   
         formatSamples(5,numIteration)=priorRange(1,5)+rand(1)*(priorRange(2,5)-priorRange(1,5));
        else
         formatSamples(5,numIteration)=delta;               
        end

        if any(6==intermediateCurrentInfluentialReactions)    
         formatSamples(6,numIteration)=priorRange(1,6)+rand(1)*(priorRange(2,6)-priorRange(1,6));
        else
         formatSamples(6,numIteration)=delta;               
        end        

        if any(7==intermediateCurrentInfluentialReactions)
         formatSamples(7,numIteration)=priorRange(1,7)+rand(1)*(priorRange(2,7)-priorRange(1,7));
        else
         formatSamples(7,numIteration)=delta;               
        end
        
        if any(8==intermediateCurrentInfluentialReactions)   
         formatSamples(8,numIteration)=priorRange(1,8)+rand(1)*(priorRange(2,8)-priorRange(1,8));
        else
         formatSamples(8,numIteration)=delta;               
        end
        
        if any(9==intermediateCurrentInfluentialReactions)   
         formatSamples(9,numIteration)=priorRange(1,9)+rand(1)*(priorRange(2,9)-priorRange(1,9));
        else
         formatSamples(9,numIteration)=delta;             
        end
        
        if any(10==intermediateCurrentInfluentialReactions)   
         formatSamples(10,numIteration)=priorRange(1,10)+rand(1)*(priorRange(2,10)-priorRange(1,10));
        else
         formatSamples(10,numIteration)=delta;               
        end        
                                                      
        [leadLogLikelihood,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,inputX,data,formatSamples(1:N,numIteration),pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);
       
        lprior=logprior(N,formatSamples(1:N,numIteration),priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateCurrentInfluentialReactions);        
        
        leadLogPosterior=leadLogLikelihood+lprior;         
        
        if typeSampling==1
            
            if ~all(sort(currentInfluentialReactions)==intermediateCurrentInfluentialReactions)
                
                diffReactions=setdiff(currentInfluentialReactions,intermediateCurrentInfluentialReactions); 
                
                logPosTrue=leadLogPosterior+logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions);
                
            else
                
                logPosTrue=leadLogPosterior;
                
            end
            
        else
            
            logPosTrue=leadLogPosterior;
            
        end
                
    end
    %%    
    
    %% 10 reaction (model problem 14)
    if modelProblem==14
        
        formatSamples=delta*ones(N,numIteration);
        
        if any(1==intermediateCurrentInfluentialReactions)
         formatSamples(1,numIteration)=priorRange(1,1)+rand(1)*(priorRange(2,1)-priorRange(1,1));
        else
         formatSamples(1,numIteration)=delta;            
        end

        if any(2==intermediateCurrentInfluentialReactions)
         formatSamples(2,numIteration)=priorRange(1,2)+rand(1)*(priorRange(2,2)-priorRange(1,2));
        else
         formatSamples(2,numIteration)=delta;     
        end
        
        if any(3==intermediateCurrentInfluentialReactions)  
         formatSamples(3,numIteration)=priorRange(1,3)+rand(1)*(priorRange(2,3)-priorRange(1,3));
        else
         formatSamples(3,numIteration)=delta;                
        end

        if any(4==intermediateCurrentInfluentialReactions)   
         formatSamples(4,numIteration)=priorRange(1,4)+rand(1)*(priorRange(2,4)-priorRange(1,4));
        else
         formatSamples(4,numIteration)=delta;              
        end
        
        if any(5==intermediateCurrentInfluentialReactions)   
         formatSamples(5,numIteration)=priorRange(1,5)+rand(1)*(priorRange(2,5)-priorRange(1,5));
        else
         formatSamples(5,numIteration)=delta;               
        end

        if any(6==intermediateCurrentInfluentialReactions)    
         formatSamples(6,numIteration)=priorRange(1,6)+rand(1)*(priorRange(2,6)-priorRange(1,6));
        else
         formatSamples(6,numIteration)=delta;               
        end        

        if any(7==intermediateCurrentInfluentialReactions)
         formatSamples(7,numIteration)=priorRange(1,7)+rand(1)*(priorRange(2,7)-priorRange(1,7));
        else
         formatSamples(7,numIteration)=delta;               
        end
        
        if any(8==intermediateCurrentInfluentialReactions)   
         formatSamples(8,numIteration)=priorRange(1,8)+rand(1)*(priorRange(2,8)-priorRange(1,8));
        else
         formatSamples(8,numIteration)=delta;               
        end
        
        if any(9==intermediateCurrentInfluentialReactions)   
         formatSamples(9,numIteration)=priorRange(1,9)+rand(1)*(priorRange(2,9)-priorRange(1,9));
        else
         formatSamples(9,numIteration)=delta;             
        end
        
        if any(10==intermediateCurrentInfluentialReactions)   
         formatSamples(10,numIteration)=priorRange(1,10)+rand(1)*(priorRange(2,10)-priorRange(1,10));
        else
         formatSamples(10,numIteration)=delta;               
        end              
      
        [leadLogLikelihood,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,inputX,data,formatSamples(1:N,numIteration),pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);
       
        lprior=logprior(N,formatSamples(1:N,numIteration),priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateCurrentInfluentialReactions);        
        
        leadLogPosterior=leadLogLikelihood+lprior;               
        
        if typeSampling==1
            
            if ~all(sort(currentInfluentialReactions)==intermediateCurrentInfluentialReactions)
                
                diffReactions=setdiff(currentInfluentialReactions,intermediateCurrentInfluentialReactions); 
                
                logPosTrue=leadLogPosterior+logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions);
                
            else
                
                logPosTrue=leadLogPosterior;
                
            end
            
        else
            
            logPosTrue=leadLogPosterior;
            
        end
        
    end
    %%     
    
    %% 2 variable linear Gaussian problem (model problem 8)
    if modelProblem==8
        
        formatSamples=delta*ones(N,numIteration);
        
        if any(1==intermediateCurrentInfluentialReactions)
         formatSamples(1,numIteration)=0.1;
        else
         formatSamples(1,numIteration)=delta;            
        end
        
        if any(2==intermediateCurrentInfluentialReactions)  
         formatSamples(2,numIteration)=0.1;
        else
         formatSamples(2,numIteration)=delta;                     
        end

         [leadLogLikelihood,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,inputX,data,formatSamples(1:N,numIteration),pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);

         lprior=logprior(N,formatSamples(1:N,numIteration),priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateCurrentInfluentialReactions);        
        
         leadLogPosterior=leadLogLikelihood+lprior;                
         
         if typeSampling==1
             
             if ~all(sort(currentInfluentialReactions)==intermediateCurrentInfluentialReactions)
                 
                 diffReactions=setdiff(currentInfluentialReactions,intermediateCurrentInfluentialReactions); 
                 
                 logPosTrue=leadLogPosterior+logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions);
                 
             else
                 
                 logPosTrue=leadLogPosterior;
                 
             end
             
         else
             
             logPosTrue=leadLogPosterior;
             
         end
         
                  
    end 
    %%    

    %% 6 dimensional linear Gaussian reaction (model problem 9)
    if modelProblem==9

        formatSamples=delta*ones(N,numIteration);

        if any(1==intermediateCurrentInfluentialReactions)   
         formatSamples(1,numIteration)=0.1;
        else
         formatSamples(1,numIteration)=delta;               
        end

        if any(2==intermediateCurrentInfluentialReactions)   
         formatSamples(2,numIteration)=0.1;
        else
         formatSamples(2,numIteration)=delta;             
        end
        
        if any(3==intermediateCurrentInfluentialReactions)    
         formatSamples(3,numIteration)=0.1;
        else
         formatSamples(3,numIteration)=delta; 
        end 
        
        if any(4==intermediateCurrentInfluentialReactions)   
         formatSamples(4,numIteration)=0.1;
        else
         formatSamples(4,numIteration)=delta;               
        end

        if any(5==intermediateCurrentInfluentialReactions)   
         formatSamples(5,numIteration)=0.1;
        else
         formatSamples(5,numIteration)=delta;             
        end
        
        if any(6==intermediateCurrentInfluentialReactions)    
         formatSamples(6,numIteration)=0.1;
        else
         formatSamples(6,numIteration)=delta; 
        end 
            
        [leadLogLikelihood,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,inputX,data,formatSamples(1:N,numIteration),pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);        

        lprior=logprior(N,formatSamples(1:N,numIteration),priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateCurrentInfluentialReactions,meanVector,valVector);        
        
        leadLogPosterior=leadLogLikelihood+lprior;              
        
        if typeSampling==1
            
            if ~all(sort(currentInfluentialReactions)==intermediateCurrentInfluentialReactions)
                
                diffReactions=setdiff(currentInfluentialReactions,intermediateCurrentInfluentialReactions); 
                
                logPosTrue=leadLogPosterior+logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions,meanVector,valVector);
                
            else
                
                logPosTrue=leadLogPosterior;
                
            end
            
        else
            
            logPosTrue=leadLogPosterior;
            
        end
        
        
        
    end 
    
    %% 3 dimensional linear Gaussian reaction (model problem 10)
    if modelProblem==10

        formatSamples=delta*ones(N,numIteration);

        if any(1==intermediateCurrentInfluentialReactions)   
         formatSamples(1,numIteration)=0.1;
        else
         formatSamples(1,numIteration)=delta;               
        end

        if any(2==intermediateCurrentInfluentialReactions)   
         formatSamples(2,numIteration)=0.1;
        else
         formatSamples(2,numIteration)=delta;             
        end
        
        if any(3==intermediateCurrentInfluentialReactions)    
         formatSamples(3,numIteration)=0.1;
        else
         formatSamples(3,numIteration)=delta; 
        end        

        [leadLogLikelihood,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,inputX,data,formatSamples(1:N,numIteration),pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay);        
      
        lprior=logprior(N,formatSamples(1:N,numIteration),priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,intermediateCurrentInfluentialReactions);        
        
        leadLogPosterior=leadLogLikelihood+lprior;          
                
        if typeSampling==1
            
            if ~all(sort(currentInfluentialReactions)==intermediateCurrentInfluentialReactions)
                
                diffReactions=setdiff(currentInfluentialReactions,intermediateCurrentInfluentialReactions); 
                
                logPosTrue=leadLogPosterior+logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions);
                
            else
                
                logPosTrue=leadLogPosterior;
                
            end
            
        else
            
            logPosTrue=leadLogPosterior;
            
        end
        
    end
