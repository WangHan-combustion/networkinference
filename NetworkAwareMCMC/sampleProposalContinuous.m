function [y]=sampleProposalContinuous(typeSampling,currentModelIndex,proposedModelIndex,modelProblem,algorithm,nReactions,reactionType,N,priorRange,d1,intermediateCurrentInfluentialReactions,intermediateProposedInfluentialReactions,proposedWithinModelMoveNumber,delta,withinModel,proposalType,centeringLocation,centeringLocationOnBoundary,eigenVectorMatrix,eigenValueMatrix,gradientVector,HessianMatrix,numberOfVariables,AreReactionAdded,K,K1,optimizedParameters)

y=zeros(N,1);

% Default alpha, beta, and gamma values

alphaParameter=0.0;
betaParameter=0.0;
gammaParameter=0.0;
normalMean=0.0;
normalSigma=0.0;


if proposalType==1

   E=diag(diag(eigenValueMatrix).^(-1/2));     
    
   % Construct the A matrix
   A=eigenVectorMatrix*E;

end

% Z: vector of standard normal variables
Z=zeros(numberOfVariables,1);

if AreReactionAdded==1 && (typeSampling~=3 || (typeSampling==3 && currentModelIndex~=proposedModelIndex))
    
    for i=1:numberOfVariables
        Z(i)=randn(1);
    end
    
end

variableNumber=0;

%% Standard algorithm (pilot parameters and conditional maximization)
if algorithm==2
            
        % normal proposal
        
        if withinModel==0
            % outside model move
            for n=1:nReactions
                
                if (any(n==intermediateProposedInfluentialReactions) && any(n==intermediateCurrentInfluentialReactions)) 
                    
                    % if n is not the proposed reaction and n present in both
                    % currentInfluentialReactions and
                    % proposedInfluentialReactions
                    
                    if reactionType(n)==1
                        % if n is a reversible reaction
                        y(n+sum(reactionType(1:n-1)))=d1(n+sum(reactionType(1:n-1)));
                        y(1+n+sum(reactionType(1:n-1)))=d1(1+n+sum(reactionType(1:n-1)));
                    else
                        % if n is an irreversible reaction
                        y(n+sum(reactionType(1:n-1)))=d1(n+sum(reactionType(1:n-1)));
                    end
                    
                else
                                       
                    % either n is present only in proposedInfluentialReactions
                    % or only in currentInfluentialReactions or not present in
                    % both
                    
                    if reactionType(n)==1
                        
                        % if n is a reversible reaction
                        if any(n==intermediateProposedInfluentialReactions)
                            
                            % Add reaction
                                                       
                            % Regular proposal
                            
                            m1=centeringLocation(variableNumber+1);
                            m2=centeringLocation(variableNumber+2);
                            
                            if proposalType==1
                                

                                   if n~=2 && n~=1 && n~=10 && n~=12 && n~=7 && n~=8 && n~=11 && n~=13
                                      y(n+sum(reactionType(1:n-1)))=m1+A(variableNumber+1,:)*Z;
                                      y(1+n+sum(reactionType(1:n-1)))=m2+A(variableNumber+2,:)*Z;                                        
                                      
                                   else
                                     if n==1
                                      y(n+sum(reactionType(1:n-1)))=0.0;
                                      y(1+n+sum(reactionType(1:n-1)))=0.0;                                          
                                     elseif n==2
                                      y(n+sum(reactionType(1:n-1)))=1.5;
                                      y(1+n+sum(reactionType(1:n-1)))=0.0; 
                                     elseif n==4
                                      y(n+sum(reactionType(1:n-1)))=2.0;
                                      y(1+n+sum(reactionType(1:n-1)))=2.0;                                          
                                     elseif n==8  
                                      y(n+sum(reactionType(1:n-1)))=0.5;
                                      y(1+n+sum(reactionType(1:n-1)))=0.5; 
                                     elseif n==9  
                                      y(n+sum(reactionType(1:n-1)))=1.0;
                                      y(1+n+sum(reactionType(1:n-1)))=1.0;                                       
                                     elseif n==10 
                                      y(n+sum(reactionType(1:n-1)))=0.0;
                                      y(1+n+sum(reactionType(1:n-1)))=0.0;                                       
                                     elseif n==11
                                      y(n+sum(reactionType(1:n-1)))=0.5;
                                      y(1+n+sum(reactionType(1:n-1)))=0.5;
                                     elseif n==12 
                                       y(n+sum(reactionType(1:n-1)))=4.0;
                                      y(1+n+sum(reactionType(1:n-1)))=4.0;                                      
                                     elseif n==13
                                      y(n+sum(reactionType(1:n-1)))=2.5;
                                      y(1+n+sum(reactionType(1:n-1)))=2.5;                                        
                                     end    
                                   end 
                                   
                                
                                variableNumber=variableNumber+2;
                                                                                                                                
                            elseif proposalType==3
                                
                               if optimizedParameters==0
                                  
                                  normalMean(variableNumber+1)=centeringLocation(variableNumber+1);
                                  normalMean(variableNumber+2)=centeringLocation(variableNumber+2);
                                  
                                  normalSigma(variableNumber+1)=sqrt(1/HessianMatrix(variableNumber+1,variableNumber+1));
                                  normalSigma(variableNumber+2)=sqrt(1/HessianMatrix(variableNumber+2,variableNumber+2));
                                  
                               else
                                  
                                  if centeringLocationOnBoundary(variableNumber+1)==1 || centeringLocationOnBoundary(variableNumber+1)==2
                                     
                                     % Hessian is of -loglikelihood
                                      
                                     normalMean(variableNumber+1)=m1+gradientVector(variableNumber+1)/HessianMatrix(variableNumber+1,variableNumber+1);
                                     normalSigma(variableNumber+1)=sqrt(1/HessianMatrix(variableNumber+1,variableNumber+1));
                                      
                                  else
                                     
                                     normalMean(variableNumber+1)=m1;
                                     normalSigma(variableNumber+1)=sqrt(1/HessianMatrix(variableNumber+1,variableNumber+1));
                                      
                                  end
                                  
                                  if centeringLocationOnBoundary(variableNumber+2)==1 || centeringLocationOnBoundary(variableNumber+2)==2
                                     
                                     normalMean(variableNumber+2)=m2+gradientVector(variableNumber+2)/HessianMatrix(variableNumber+2,variableNumber+2);
                                     normalSigma(variableNumber+2)=sqrt(1/HessianMatrix(variableNumber+2,varialbeNumber+2));
                                      
                                  else
                                     
                                     normalMean(variableNumber+2)=m2;
                                     normalSigma(variableNumber+2)=sqrt(1/HessianMatrix(variableNumber+2,variableNumber+2));
                                      
                                  end                                  
                                   
                               end     
                                   
                                                                

                                   
                                   p1=normcdf(priorRange(1,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1))+rand(1)*(normcdf(priorRange(2,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1))-normcdf(priorRange(1,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1)));
                                   
                                   if (centeringLocationOnBoundary(variableNumber+1)==1) && (normcdf(priorRange(2,n+sum(reactionType(1:n-1))))==1 || (normcdf(priorRange(2,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1))-normcdf(priorRange(1,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1)))<=1E-10)
                                     % Uniform proposal
                                     y(n+sum(reactionType(1:n-1)))=priorRange(1,n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))));    
                                   elseif (centeringLocationOnBoundary(variableNumber+1)==2) && (normcdf(priorRange(1,n+sum(reactionType(1:n-1))))==0 || (normcdf(priorRange(2,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1))-normcdf(priorRange(1,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1)))<=1E-10)
                                     % Uniform proposal
                                     y(n+sum(reactionType(1:n-1)))=priorRange(1,n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))));                                          
                                   else  
                                     % Truncated normal proposal  
                                     y(n+sum(reactionType(1:n-1)))=norminv(p1,normalMean(variableNumber+1),normalSigma(variableNumber+1));                                       
                                   end    
                                   
                                   p2=normcdf(priorRange(1,1+n+sum(reactionType(1:n-1))),normalMean(variableNumber+2),normalSigma(variableNumber+2))+rand(1)*(normcdf(priorRange(2,1+n+sum(reactionType(1:n-1))),normalMean(variableNumber+2),normalSigma(variableNumber+2))-normcdf(priorRange(1,1+n+sum(reactionType(1:n-1))),normalMean(variableNumber+2),normalSigma(variableNumber+2)));
                                         
                                   if (centeringLocationOnBoundary(variableNumber+2)==1) && (normcdf(priorRange(2,1+n+sum(reactionType(1:n-1))))==1 || (normcdf(priorRange(2,1+n+sum(reactionType(1:n-1))),normalMean(variableNumber+2),normalSigma(variableNumber+2))-normcdf(priorRange(1,1+n+sum(reactionType(1:n-1))),normalMean(variableNumber+2),normalSigma(variableNumber+2)))<=1E-10)                                   
                                     % Uniform proposal
                                     y(1+n+sum(reactionType(1:n-1)))=priorRange(1,1+n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,1+n+sum(reactionType(1:n-1)))-priorRange(1,1+n+sum(reactionType(1:n-1))));                                   
                                   elseif (centeringLocationOnBoundary(variableNumber+2)==2) && (normcdf(priorRange(1,1+n+sum(reactionType(1:n-1))))==0 || (normcdf(priorRange(2,1+n+sum(reactionType(1:n-1))),normalMean(variableNumber+2),normalSigma(variableNumber+2))-normcdf(priorRange(1,1+n+sum(reactionType(1:n-1))),normalMean(variableNumber+2),normalSigma(variableNumber+2)))<=1E-10)
                                     % Uniform proposal
                                     y(1+n+sum(reactionType(1:n-1)))=priorRange(1,1+n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,1+n+sum(reactionType(1:n-1)))-priorRange(1,1+n+sum(reactionType(1:n-1))));                                                                                                   
                                   else  
                                     % Truncated normal proposal  
                                     y(1+n+sum(reactionType(1:n-1)))=norminv(p2,normalMean(variableNumber+2),normalSigma(variableNumber+2));                                                                      
                                   end    
   
                                
                               variableNumber=variableNumber+2;
                               
                            end
                                                                                                                                                                                                                                   
                        else
                            
                            % delete reaction
                            
                            y(n+sum(reactionType(1:n-1)))=delta;
                            y(n+sum(reactionType(1:n-1))+1)=delta;
                            
                            
                        end
                    else
                        
                        % if n is irreversible
                        
                        % either n is present only in proposedInfluentialReactions
                        % or only in currentInfluentialReactions or not present in
                        % both
                        
                        if any(n==intermediateProposedInfluentialReactions)
                            
                            % Add reaction
                                                                                    
                            m1=centeringLocation(variableNumber+1);
                            
                            if proposalType==1
                            
                                if n~=1 && n~=2 && n~=10 && n~=12 && n~=8 && n~=11 && n~=13 && n~=7
                                    y(n+sum(reactionType(1:n-1)))=m1+A(variableNumber+1,:)*Z;
                                    
                                else
                                    if n==1
                                        y(n+sum(reactionType(1:n-1)))=0.0;                                          
                                    elseif n==2
                                        y(n+sum(reactionType(1:n-1)))=1.5;    
                                    elseif n==4
                                        y(n+sum(reactionType(1:n-1)))=2.0;                                         
                                    elseif n==8
                                        y(n+sum(reactionType(1:n-1)))=0.5;
                                    elseif n==9
                                        y(n+sum(reactionType(1:n-1)))=1.0;                                        
                                    elseif n==10    
                                        y(n+sum(reactionType(1:n-1)))=0.0;                                          
                                    elseif n==11
                                        y(n+sum(reactionType(1:n-1)))=0.5;
                                    elseif n==12    
                                        y(n+sum(reactionType(1:n-1)))=4.0;                                          
                                    elseif n==13
                                        y(n+sum(reactionType(1:n-1)))=2.5;
                                    end
                                    
                                end
                                
                            
                                 variableNumber=variableNumber+1;                            
                               
                            elseif proposalType==3
                                
                               % Truncated normal proposal 
                                
                               if optimizedParameters==0
                                  
                                  normalMean(variableNumber+1)=centeringLocation(variableNumber+1);
                                  normalSigma(variableNumber+1)=sqrt(1/HessianMatrix(variableNumber+1,variableNumber+1));
                                  
                               else

                                  if centeringLocationOnBoundary(variableNumber+1)==1 || centeringLocationOnBoundary(variableNumber+1)==2
                                     
                                     % Hessian is postive definite 
                                      
                                     normalMean(variableNumber+1)=m1+gradientVector(variableNumber+1)/HessianMatrix(variableNumber+1,variableNumber+1);
                                     normalSigma(variableNumber+1)=sqrt(1/HessianMatrix(variableNumber+1,variableNumber+1));
                                      
                                  else
                                     
                                     % Hessian is positive definite 
                                      
                                     normalMean(variableNumber+1)=m1;
                                     normalSigma(variableNumber+1)=sqrt(1/HessianMatrix(variableNumber+1,variableNumber+1));
                                      
                                  end     
                                                                     
                               end                                 
                                
                                

                                                                                                                                                                                                                                                               
                                   p1=normcdf(priorRange(1,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1))+rand(1)*(normcdf(priorRange(2,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1))-normcdf(priorRange(1,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1)));
                                 
                                   if (centeringLocationOnBoundary(variableNumber+1)==1) && (normcdf(priorRange(2,n+sum(reactionType(1:n-1))))==1 || (normcdf(priorRange(2,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1))-normcdf(priorRange(1,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1)))<=1E-10)
                                     % Uniform proposal
                                     y(n+sum(reactionType(1:n-1)))=priorRange(1,n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))));    
                                   elseif (centeringLocationOnBoundary(variableNumber+1)==2) && (normcdf(priorRange(1,n+sum(reactionType(1:n-1))))==0 || (normcdf(priorRange(2,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1))-normcdf(priorRange(1,n+sum(reactionType(1:n-1))),normalMean(variableNumber+1),normalSigma(variableNumber+1)))<=1E-10)
                                     % Uniform proposal  
                                     y(n+sum(reactionType(1:n-1)))=priorRange(1,n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))));                                          
                                   else  
                                     % Truncated normal proposal  
                                     y(n+sum(reactionType(1:n-1)))=norminv(p1,normalMean(variableNumber+1),normalSigma(variableNumber+1));                                       
                                   end                                   
                             
                                
                                variableNumber=variableNumber+1;
                               
                            end    
                                                                                                               
                        else
                            
                            % Delete reaction
                            y(n+sum(reactionType(1:n-1)))=delta;
                            
                        end
                    end   
                    
                end
            end
                   
        elseif withinModel==1
            
           
            % Gaussian random walk
            % within model move
            if proposedWithinModelMoveNumber==0
                y=d1;
            else
                for n=1:nReactions
                    % if n is the proposed reaction
                    %if n==proposedWithinModelMoveNumber
                    if any(intermediateProposedInfluentialReactions==n)
                        if reactionType(n)==1

                           if n~=1 && n~=2 %&& n~=4 && n~=9 && n~=10 && n~=12 && n~=7                            
                            y(n+sum(reactionType(1:n-1)))=d1(n+sum(reactionType(1:n-1)))+(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))))/20*randn(1);
                            y(n+sum(reactionType(1:n-1))+1)=d1(1+n+sum(reactionType(1:n-1)))+(priorRange(2,1+n+sum(reactionType(1:n-1)))-priorRange(1,1+n+sum(reactionType(1:n-1))))/20*randn(1);                           
                           
                           else
                                     if n==1
                                      y(n+sum(reactionType(1:n-1)))=0.0;
                                      y(1+n+sum(reactionType(1:n-1)))=0.0;                                          
                                     elseif n==2
                                      y(n+sum(reactionType(1:n-1)))=1.5;
                                      y(1+n+sum(reactionType(1:n-1)))=0.0;             
                                     elseif n==4
                                      y(n+sum(reactionType(1:n-1)))=2.0;
                                      y(1+n+sum(reactionType(1:n-1)))=2.0;                                          
                                     elseif n==8  
                                      y(n+sum(reactionType(1:n-1)))=0.5;
                                      y(1+n+sum(reactionType(1:n-1)))=0.5; 
                                     elseif n==9  
                                      y(n+sum(reactionType(1:n-1)))=1.0;
                                      y(1+n+sum(reactionType(1:n-1)))=1.0;                                       
                                     elseif n==10 
                                      y(n+sum(reactionType(1:n-1)))=0.0;
                                      y(1+n+sum(reactionType(1:n-1)))=0.0;                                       
                                     elseif n==11
                                      y(n+sum(reactionType(1:n-1)))=0.5;
                                      y(1+n+sum(reactionType(1:n-1)))=0.5;
                                     elseif n==12 
                                       y(n+sum(reactionType(1:n-1)))=4.0;
                                      y(1+n+sum(reactionType(1:n-1)))=4.0;                                      
                                     elseif n==13
                                      y(n+sum(reactionType(1:n-1)))=2.5;
                                      y(1+n+sum(reactionType(1:n-1)))=2.5;                                        
                                     end                                   
                           end    
                            
                        else
                            if n~=1 && n~=2 %&& n~=4 && n~=9 && n~=10 && n~=12 && n~=7 
                             y(n+sum(reactionType(1:n-1)))=d1(n+sum(reactionType(1:n-1)))+(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))))/20*randn(1);
                             
                            else
                                    if n==1
                                        y(n+sum(reactionType(1:n-1)))=0.0;                                          
                                    elseif n==2
                                        y(n+sum(reactionType(1:n-1)))=1.5;   
                                    elseif n==4
                                        y(n+sum(reactionType(1:n-1)))=2.0;                                         
                                    elseif n==8
                                        y(n+sum(reactionType(1:n-1)))=0.5;
                                    elseif n==9
                                        y(n+sum(reactionType(1:n-1)))=1.0;                                        
                                    elseif n==10    
                                        y(n+sum(reactionType(1:n-1)))=0.0;                                          
                                    elseif n==11
                                        y(n+sum(reactionType(1:n-1)))=0.5;
                                    elseif n==12    
                                        y(n+sum(reactionType(1:n-1)))=4.0;                                          
                                    elseif n==13
                                        y(n+sum(reactionType(1:n-1)))=2.5;
                                    end                               
                            end    
                            
                        end
                        
                        
                    else
                        % if n is not the proposed reaction
                        if reactionType(n)==1
                            % if n is a reversible reaction
                            y(n+sum(reactionType(1:n-1)))=d1(n+sum(reactionType(1:n-1)));
                            y(n+sum(reactionType(1:n-1))+1)=d1(n+sum(reactionType(1:n-1))+1);
                        else
                            % if n is an irreversible reaction
                            y(n+sum(reactionType(1:n-1)))=d1(n+sum(reactionType(1:n-1)));
                        end
                    end
                    
                end
            end
                       
        end

end    
