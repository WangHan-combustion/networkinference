function [priorRange]=priorSetup(pos,modelProblem,N,data,valVector)
    
    %% 13 reaction (model problem 1)
    if modelProblem==1
         
        % Explicit prior specification

            %{
        
            % Reaction 1
            priorRange(1,1)=-0.5;
            priorRange(2,1)=0.5;
            
            % Reaction 2
            priorRange(1,2)=1.0;
            priorRange(2,2)=2.0;
            
            priorRange(1,3)=-0.5;
            priorRange(2,3)=0.5;
            
            % Reaction 3
            priorRange(1,4)=0.0;
            priorRange(2,4)=1.0;
            
            % Reaction 4
            priorRange(1,5)=1.5;
            priorRange(2,5)=2.5;
            
            % Reaction 5
            priorRange(1,6)=2.0;
            priorRange(2,6)=3.0;
            
            % Reaction 6
            priorRange(1,7)=-0.5;
            priorRange(2,7)=0.5;
            
            % Reaction 7
            priorRange(1,8)=2.5;
            priorRange(2,8)=3.5;
            
            % Reaction 8          %%
            priorRange(1,9)=0.0;
            priorRange(2,9)=1.0;
            
            % Reaction 9
            priorRange(1,10)=0.5;
            priorRange(2,10)=1.5;
            
            % Reaction 10         %%
            priorRange(1,11)=-0.5;
            priorRange(2,11)=0.5;
            
            % Reaction 11         %%
            priorRange(1,12)=0.0;
            priorRange(2,12)=1.0;
            
            % Reaction 12         %%
            priorRange(1,13)=3.5;
            priorRange(2,13)=4.5;
            
            % Reaction 13         %%
            priorRange(1,14)=2.0;
            priorRange(2,14)=3.0;
  
            %}
            
            % Reaction 1
            priorRange(1,1)=0.0;
            priorRange(2,1)=10.0;
            
            % Reaction 2
            priorRange(1,2)=0.0;
            priorRange(2,2)=10.0;
            
            priorRange(1,3)=0.0;
            priorRange(2,3)=10.0;
            
            % Reaction 3
            priorRange(1,4)=0.0;
            priorRange(2,4)=10.0;
            
            % Reaction 4
            priorRange(1,5)=0.0;
            priorRange(2,5)=10.0;
            
            % Reaction 5
            priorRange(1,6)=0.0;
            priorRange(2,6)=10.0;
            
            % Reaction 6
            priorRange(1,7)=0.0;
            priorRange(2,7)=10.0;
            
            % Reaction 7
            priorRange(1,8)=0.0;
            priorRange(2,8)=10.0;
            
            % Reaction 8          %%
            priorRange(1,9)=0.0;
            priorRange(2,9)=10.0;
            
            % Reaction 9
            priorRange(1,10)=0.0;
            priorRange(2,10)=10.0;
            
            % Reaction 10         %%
            priorRange(1,11)=0.0;
            priorRange(2,11)=10.0;
            
            % Reaction 11         %%
            priorRange(1,12)=0.0;
            priorRange(2,12)=10.0;
            
            % Reaction 12         %%
            priorRange(1,13)=0.0;
            priorRange(2,13)=10.0;
            
            % Reaction 13         %%
            priorRange(1,14)=0.0;
            priorRange(2,14)=10.0;            
            
            
            
            
    %% 13 reaction (model problem 11)   
    elseif modelProblem==11

            m=1;
                
            % Reaction 1
            priorRange(1,1)=0.0-sqrt(valVector(1))*m;
            priorRange(2,1)=0.0+sqrt(valVector(1))*m;
            
            % Reaction 2
            priorRange(1,2)=1.5-sqrt(valVector(2))*m;
            priorRange(2,2)=1.5+sqrt(valVector(2))*m;
            
            priorRange(1,3)=0.0-sqrt(valVector(3))*m;
            priorRange(2,3)=0.0+sqrt(valVector(3))*m;
            
            % Reaction 3
            priorRange(1,4)=0.5-sqrt(valVector(4))*m;
            priorRange(2,4)=0.5+sqrt(valVector(4))*m;
            
            % Reaction 4
            priorRange(1,5)=2.0-sqrt(valVector(5))*m;
            priorRange(2,5)=2.0+sqrt(valVector(5))*m;
            
            % Reaction 5
            priorRange(1,6)=2.0-sqrt(valVector(6))*m;
            priorRange(2,6)=2.0+sqrt(valVector(6))*m;
            
            % Reaction 6
            priorRange(1,7)=0.4-sqrt(valVector(7))*m;
            priorRange(2,7)=0.4+sqrt(valVector(7))*m;
            
            % Reaction 7
            priorRange(1,8)=3.0-sqrt(valVector(8))*m;
            priorRange(2,8)=3.0+sqrt(valVector(8))*m;
            
            % Reaction 8          %%
            priorRange(1,9)=0.5-sqrt(valVector(9))*m;
            priorRange(2,9)=0.5+sqrt(valVector(9))*m;
            
            % Reaction 9
            priorRange(1,10)=1.0-sqrt(valVector(10))*m;
            priorRange(2,10)=1.0+sqrt(valVector(10))*m;
            
            % Reaction 10         %%
            priorRange(1,11)=0.0-sqrt(valVector(11))*m;
            priorRange(2,11)=0.0+sqrt(valVector(11))*m;
            
            % Reaction 11         %%
            priorRange(1,12)=0.5-sqrt(valVector(12))*m;
            priorRange(2,12)=0.5+sqrt(valVector(12))*m;
            
            % Reaction 12         %%
            priorRange(1,13)=4.0-sqrt(valVector(13))*m;
            priorRange(2,13)=4.0+sqrt(valVector(13))*m;
            
            % Reaction 13         %%
            priorRange(1,14)=2.5-sqrt(valVector(14))*m;
            priorRange(2,14)=2.5+sqrt(valVector(14))*m;                              
       
       
    %% 30 reaction (model problem 2)
    elseif modelProblem==2
        
        
        % Explicit prior specification
        priorRange(1,1)=0;
        priorRange(2,1)=2;
        
        priorRange(1,2)=0;
        priorRange(2,2)=2;
        
        priorRange(1,3)=0;
        priorRange(2,3)=2;
        
        priorRange(1,4)=0;
        priorRange(2,4)=2;
        
        priorRange(1,5)=0;
        priorRange(2,5)=2;
        
        priorRange(1,6)=0;
        priorRange(2,6)=2;
        
        priorRange(1,7)=0;
        priorRange(2,7)=2;
        
        priorRange(1,8)=0;
        priorRange(2,8)=2;
        
        priorRange(1,9)=0;
        priorRange(2,9)=2;
        
        priorRange(1,10)=0;
        priorRange(2,10)=2;
        
        priorRange(1,11)=0;
        priorRange(2,11)=2;
        
        priorRange(1,12)=0;
        priorRange(2,12)=2;
        
        priorRange(1,13)=0;
        priorRange(2,13)=2;
        
        priorRange(1,14)=0;
        priorRange(2,14)=2;

        priorRange(1,15)=0;
        priorRange(2,15)=2;
        
        priorRange(1,16)=0;
        priorRange(2,16)=2;
        
        priorRange(1,17)=0;
        priorRange(2,17)=2;
        
        priorRange(1,18)=0;
        priorRange(2,18)=2;
        
        priorRange(1,19)=0;
        priorRange(2,19)=2;
        
        priorRange(1,20)=0;
        priorRange(2,20)=2;
        
        priorRange(1,21)=0;
        priorRange(2,21)=2;
        
        priorRange(1,22)=0;
        priorRange(2,22)=2;
        
        priorRange(1,23)=0;
        priorRange(2,23)=2;
        
        priorRange(1,24)=0;
        priorRange(2,24)=2;
        
        priorRange(1,25)=0;
        priorRange(2,25)=2;
        
        priorRange(1,26)=0;
        priorRange(2,26)=2;
        
        priorRange(1,27)=0;
        priorRange(2,27)=2;
        
        priorRange(1,28)=0;
        priorRange(2,28)=2;
        
        priorRange(1,29)=0;
        priorRange(2,29)=2;
        
        priorRange(1,30)=0;
        priorRange(2,30)=2;
        
        priorRange(1,31)=0;
        priorRange(2,31)=2;
        
        
    %%          

    %% 20 reaction (model problem 3)
    elseif modelProblem==3
        
        
        % Explicit prior specification
        priorRange(1,1)=0.9;
        priorRange(2,1)=1.1;
        
        priorRange(1,2)=0.9;
        priorRange(2,2)=1.1;
        
        priorRange(1,3)=0.9;
        priorRange(2,3)=1.1;
        
        priorRange(1,4)=0.9;
        priorRange(2,4)=1.1;
        
        priorRange(1,5)=0.3;
        priorRange(2,5)=0.5;
        
        priorRange(1,6)=0.1;
        priorRange(2,6)=0.3;
        
        priorRange(1,7)=0.4;
        priorRange(2,7)=0.6;
        
        priorRange(1,8)=0;
        priorRange(2,8)=0.2;
        
        priorRange(1,9)=0.4;
        priorRange(2,9)=0.6;
        
        priorRange(1,10)=0.5;
        priorRange(2,10)=0.7;
        
        priorRange(1,11)=0.2;
        priorRange(2,11)=0.4;
        
        priorRange(1,12)=1.9;
        priorRange(2,12)=2.1;
        
        priorRange(1,13)=0.9;
        priorRange(2,13)=1.1;
        
        priorRange(1,14)=0.6;
        priorRange(2,14)=0.8;

        priorRange(1,15)=0.2;
        priorRange(2,15)=0.4;
        
        priorRange(1,16)=0.9;
        priorRange(2,16)=1.1;
        
        priorRange(1,17)=0.4;
        priorRange(2,17)=0.6;
        
        priorRange(1,18)=1.1;
        priorRange(2,18)=1.3;
        
                      
    %%              
    
    
    
    %% 10 reaction (model problem 4)
    elseif modelProblem==4
        
        priorRange(1,1)=0;
        priorRange(2,1)=2;
    
        priorRange(1,2)=1;
        priorRange(2,2)=3;
    
        priorRange(1,3)=2;
        priorRange(2,3)=4;
    
        priorRange(1,4)=1;
        priorRange(2,4)=3;
    
        priorRange(1,5)=3;
        priorRange(2,5)=5;
    
        priorRange(1,6)=0.5;
        priorRange(2,6)=2.5;
    
        priorRange(1,7)=1.0;
        priorRange(2,7)=3.0;
    
        priorRange(1,8)=1.5;
        priorRange(2,8)=3.5;
    
        priorRange(1,9)=0;
        priorRange(2,9)=2;
    
        priorRange(1,10)=2.0;
        priorRange(2,10)=4.0;
        
    %% 10 reaction (model problem 14)    
    elseif modelProblem==14
        
        priorRange(1,1)=0;
        priorRange(2,1)=2;
    
        priorRange(1,2)=1;
        priorRange(2,2)=3;
    
        priorRange(1,3)=2;
        priorRange(2,3)=4;
    
        priorRange(1,4)=1;
        priorRange(2,4)=3;
    
        priorRange(1,5)=3;
        priorRange(2,5)=5;
    
        priorRange(1,6)=0.5;
        priorRange(2,6)=2.5;
    
        priorRange(1,7)=1.0;
        priorRange(2,7)=3.0;
    
        priorRange(1,8)=1.5;
        priorRange(2,8)=3.5;
    
        priorRange(1,9)=0;
        priorRange(2,9)=2;
    
        priorRange(1,10)=2.0;
        priorRange(2,10)=4.0;        
        
            
    %% 6 reaction (model problem 5)
    elseif modelProblem==5
        
        % Explicit prior specification

            % Reaction 1
            priorRange(1,1)=-4.5;
            priorRange(2,1)=-3.5;
            
            % Reaction 2
            priorRange(1,2)=0.5;
            priorRange(2,2)=1.5;
            
            priorRange(1,3)=0.0;
            priorRange(2,3)=1.0;
            
            % Reaction 3
            priorRange(1,4)=0.0;
            priorRange(2,4)=1.0;
            
            % Reaction 4
            priorRange(1,5)=-1.0;
            priorRange(2,5)=-0.0;
            
            % Reaction 5
            priorRange(1,6)=-0.1;
            priorRange(2,6)=0.9;
            
            % Reaction 6
            priorRange(1,7)=-0.1;
            priorRange(2,7)=0.9;
            
            % Reaction 7
            priorRange(1,8)=1.0;
            priorRange(2,8)=2.0;
            
            % Reaction 8          %%
            priorRange(1,9)=0.0;
            priorRange(2,9)=1.0;
            
            % Reaction 9
            priorRange(1,10)=-4.0;
            priorRange(2,10)=-3.0;
            
            % Reaction 10         %%
            priorRange(1,11)=-1.0;
            priorRange(2,11)=-0.0;
            
            % Reaction 11         %%
            priorRange(1,12)=1.0;
            priorRange(2,12)=2.0;
            
            % Reaction 12         %%
            priorRange(1,13)=-4.5;
            priorRange(2,13)=-3.5;
            
            % Reaction 13         %%
            priorRange(1,14)=2.0;
            priorRange(2,14)=3.0;
  

    %% 5 reaction (model problem 6)
    elseif modelProblem==6
        
        priorRange(1,1)=0;
        priorRange(2,1)=2;
    
        priorRange(1,2)=1;
        priorRange(2,2)=3;
    
        priorRange(1,3)=2;
        priorRange(2,3)=4;
    
        priorRange(1,4)=1;
        priorRange(2,4)=3;
    
        priorRange(1,5)=3;
        priorRange(2,5)=5;
    
        priorRange(1,6)=0.5;
        priorRange(2,6)=2.5;
    
        priorRange(1,7)=1.0;
        priorRange(2,7)=3.0;
    
        priorRange(1,8)=1.5;
        priorRange(2,8)=3.5;
    
        priorRange(1,9)=0;
        priorRange(2,9)=2;
    
        priorRange(1,10)=2.0;
        priorRange(2,10)=4.0;        

    
    %% 2 dimensional variable selection (model problem 7)
    elseif modelProblem==7    
        
        
        priorRange(1,1)=-3;
        priorRange(2,1)=3;
        
        priorRange(1,2)=-3;
        priorRange(2,2)=3;
    
    %% 2 linear Gaussian variable selection (model problem 8)     
    elseif modelProblem==8 
                             
        priorRange(1,1)=0;
        priorRange(2,1)=2;
    
        priorRange(1,2)=0;
        priorRange(2,2)=1;      
    
    %% 6 dimensional linear Gaussian variable selection (model problem 9)    
    elseif modelProblem==9    
                              
        priorRange(1,1)=-5;
        priorRange(2,1)=5;
    
        priorRange(1,2)=-5;
        priorRange(2,2)=5;
    
        priorRange(1,3)=-5;
        priorRange(2,3)=5;
    
        priorRange(1,4)=-5;
        priorRange(2,4)=5;
    
        priorRange(1,5)=-5;
        priorRange(2,5)=5;
    
        priorRange(1,6)=-5;
        priorRange(2,6)=5;        
    
    elseif modelProblem==10
                
        
        priorRange(1,1)=-5;
        priorRange(2,1)=5;
    
        priorRange(1,2)=-5;
        priorRange(2,2)=5;
    
        priorRange(1,3)=-5;
        priorRange(2,3)=5;        
        
    end
    
