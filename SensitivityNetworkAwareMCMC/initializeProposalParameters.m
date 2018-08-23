function [mVector,varVector,wVector]=initializeProposalParameters(priorRange,delta,N,nReactions,reactionType,K,K1)

% Parameters being updated
mVector=zeros(K,N);
varVector=zeros(K,N);

wVector=zeros(K,nReactions);

%===================================
% Initialization of delta components
 mVector(1:K1,1:N)=delta;
 varVector(1:K1,1:N)=0.0;
%===================================

%==========================================================================
% Initialization of the parameters

% Initialization of mean
for n=1:nReactions

    if reactionType(n)==1
        for k=K1+1:K
            mVector(k,n+sum(reactionType(1:n-1)))=priorRange(1,n+sum(reactionType(1:n-1)))+(k-K1)*(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))))/K;
            mVector(k,1+n+sum(reactionType(1:n-1)))=priorRange(1,1+n+sum(reactionType(1:n-1)))+(k-K1)*(priorRange(2,1+n+sum(reactionType(1:n-1)))-priorRange(1,1+n+sum(reactionType(1:n-1))))/K;     
            
            varVector(k,n+sum(reactionType(1:n-1)))=((priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))))/2)^2;
            varVector(k,1+n+sum(reactionType(1:n-1)))=((priorRange(2,1+n+sum(reactionType(1:n-1)))-priorRange(1,1+n+sum(reactionType(1:n-1))))/2)^2;
            
        end
    else
        for k=K1+1:K
            mVector(k,n+sum(reactionType(1:n-1)))=priorRange(1,n+sum(reactionType(1:n-1)))+(k-K1)*(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))))/K;
            varVector(k,n+sum(reactionType(1:n-1)))=((priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))))/2)^2;            
        end        
    end
        
end
    
%==========================================================================


for n=1:nReactions
   wVector(1:K1,n)=0.5*drchrnd(K1);
   wVector(K1+1:K,n)=0.5*drchrnd(K-K1);
end
