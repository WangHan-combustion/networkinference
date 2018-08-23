function [lnprior]=logAdditionalPrior(nReactions,reactionType,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,diffReactions,meanVector,valVector)

% signalling pathway example


if modelProblem==8
    
    lnprior=0.0;
    
    for n=1:length(diffReactions)
       
          
           if reactionType(diffReactions(n))==1 
           
             val1=sqrt(1/0.1)*randn(1);  
             lnprior=lnprior+lognormal(val1,0,1/0.1);        
             
             val2=sqrt(1/0.1)*randn(1);
             lnprior=lnprior+lognormal(val2,0,1/0.1);
                                      
           else
             
             val1=sqrt(1/0.1)*randn(1);
             lnprior=lnprior+lognormal(val1,0,1/0.1);
               
           end    

        
    end
    
elseif modelProblem==9
    
    lnprior=0.0;
    
    for n=1:length(diffReactions)
        

          
           if reactionType(diffReactions(n))==1 
           
             val1=sqrt(1/0.1)*randn(1);  
             lnprior=lnprior+lognormal(val1,0,1/0.1);        
             
             val2=sqrt(1/0.1)*randn(1);
             lnprior=lnprior+lognormal(val2,0,1/0.1);
                                      
           else
             
             val1=sqrt(1/0.1)*randn(1);
             lnprior=lnprior+lognormal(val1,0,1/0.1);
               
           end    

        
    end
    
        
    
elseif modelProblem==10
    
    lnprior=0.0;
    
    for n=1:length(diffReactions)
        

          
           if reactionType(diffReactions(n))==1 
           
             val1=sqrt(1/0.1)*randn(1);  
             lnprior=lnprior+lognormal(val1,0,1/0.1);        
             
             val2=sqrt(1/0.1)*randn(1);
             lnprior=lnprior+lognormal(val2,0,1/0.1);
                                      
           else
             
             val1=sqrt(1/0.1)*randn(1);
             lnprior=lnprior+lognormal(val1,0,1/0.1);
               
           end    
            

        
    end
    
            
elseif modelProblem==11
    
    
    lnprior=0.0;
    
    for m=1:length(diffReactions)
        
           n=diffReactions(m);
          
           if reactionType(n)==1 
           
             val1=meanVector(n+sum(reactionType(1:n-1)))+sqrt(valVector(n+sum(reactionType(1:n-1))))*randn(1); 
             
             if n==1
                 val1=0.0;    
             elseif n==2
                 val1=1.5; 
             %elseif n==8
             %    val1=0.5;
             %elseif n==10
             %    val1=0.0;
             %elseif n==11
             %    val1=0.5;
             %elseif n==12
             %    val1=4.0;
             %elseif n==13
             %    val1=2.5;
             %elseif n==4
             %    val1=2.0;
             %elseif n==9
             %    val1=1.0;
             end
                 
             lnprior=lnprior+lognormal(val1,meanVector(n+sum(reactionType(1:n-1))),valVector(n+sum(reactionType(1:n-1))));        
             
             val2=meanVector(n+sum(reactionType(1:n-1))+1)+sqrt(valVector(n+sum(reactionType(1:n-1))+1))*randn(1);
             
             if n==1
                 val2=0.0;
             elseif n==2
                 val2=0.0;
             %elseif n==8
             %    val2=0.5;
             %elseif n==10
             %    val2=0.0;
             %elseif n==11
             %    val2=0.5;
             %elseif n==12
             %    val2=4.0;
             %elseif n==13
             %    val2=2.5;
             %elseif n==4
             %    val2=2.0;
             %elseif n==9
             %    val2=1.0;
             end
                          
             
             lnprior=lnprior+lognormal(val2,meanVector(n+sum(reactionType(1:n-1))+1),valVector(n+sum(reactionType(1:n-1))+1));
                                      
           else
             
             val1=meanVector(n+sum(reactionType(1:n-1)))+sqrt(valVector(n+sum(reactionType(1:n-1))))*randn(1); 
             
             if n==1
                 val1=0.0;    
             elseif n==2
                 val1=1.5; 
             %elseif n==8
             %    val1=0.5;
             %elseif n==10
             %    val1=0.0;
             %elseif n==11
             %    val1=0.5;
             %elseif n==12
             %    val1=4.0;
             %elseif n==13
             %    val1=2.5;
             %elseif n==4
             %    val1=2.0;
             %elseif n==9
             %    val1=1.0;
             end        
                                       
             lnprior=lnprior+lognormal(val1,meanVector(n+sum(reactionType(1:n-1))),valVector(n+sum(reactionType(1:n-1)))); 
               
           end    
                  
    end
    
elseif modelProblem==14
    
    lnprior=0.0;
    
    for m=1:length(diffReactions)       
          
           n=diffReactions(m);
        
           if reactionType(n)==1 
           
             val1=sqrt(1/0.5)*randn(1);  
             lnprior=lnprior+lognormal(val1,0,1/0.5);        
             
             val2=sqrt(1/0.5)*randn(1);
             lnprior=lnprior+lognormal(val2,0,1/0.5);
                                      
           else
             
             val1=sqrt(1/0.5)*randn(1);
             lnprior=lnprior+lognormal(val1,0,1/0.5);
               
           end    

        
    end
    
    
    
else
    
    lnprior=0.0;
    
    %{
    for n=1:N
        if val(n)==delta
            lnprior=lnprior+log(0.5);
        else
            lnprior=lnprior+log(0.5)-log(priorRange(2,n)-priorRange(1,n));
        end
    end
    %}
    
    for m=1:length(diffReactions)
        
           n=diffReactions(m);
          
           if reactionType(n)==1 
           
             %val1=priorRange(1,n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))));
             lnprior=lnprior-log(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))));        
             
             %val2=priorRange(1,1+n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,1+n+sum(reactionType(1:n-1)))-priorRange(1,1+n+sum(reactionType(1:n-1))));
             lnprior=lnprior-log(priorRange(2,1+n+sum(reactionType(1:n-1)))-priorRange(1,1+n+sum(reactionType(1:n-1))));
                                      
           else
             
             %val1=priorRange(1,n+sum(reactionType(1:n-1)))+rand(1)*(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))));
             lnprior=lnprior-log(priorRange(2,n+sum(reactionType(1:n-1)))-priorRange(1,n+sum(reactionType(1:n-1))));    
               
           end    
            
        
    end    
    
    
    
end

if lnprior==-Inf
    lnprior=-1E300;
end

if lnprior==Inf
    lnprior=1E300;
end

