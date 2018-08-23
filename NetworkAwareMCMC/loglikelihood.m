function [l,numForwardModelSolves]=loglikelihood(numForwardModelSolves,typeSampling,input,data,valVector,pos,modelProblem,dataVariance,nReactions,nSpecies,shouldDisplay)

% This function calculates the loglikelihood given experimental data 
% and a forward model
    
    % Directing the loglikelihood evaluation for model problem 11 and 14 to
    % 1 and 4 respectively.
    if modelProblem==11
       modelProblem=1; 
    end
    
    if modelProblem==14
       modelProblem=4;  
    end
    
    
    %% (model problem 1)
    if modelProblem==1
        
        %{
        p56=valVector(1);
        p5=valVector(2);
        p6=valVector(3);
        p47=valVector(4);
        p49=valVector(5);
        p50=valVector(6);
        p41=valVector(7);
        p44=valVector(8);
        p52=valVector(9);
        p39=valVector(10);
        p11=valVector(11);
        p9=valVector(12);
        p8=valVector(13);
        p3=valVector(14);
        %}
        
        % (Reaction 7)
        p13=log10(2.24);
        % (Reaction 8)
        p16=log10(1.299);
        % (Reaction 9)
        p17=log10(4239.6527);
        % (Reaction 10)
        p20=log10(4.443);
        % (Reaction 11)
        p21=log10(10188.05);
        % (Reaction 12 and Reaction 13)
        p1=-100;
        % (Reaction 14)
        p23=-100;
        % (Reaction 15)
        p25=log10(44.58);
        % (Reaction 16)
        p27=log10(6064);
        % (Reaction 17)
        p30=log10(0.995);
        % (Reaction 18)
        p31=log10(810.869);
        % (Reaction 19)
        p33=log10(1.21);
        % (Reaction 20)
        p36=log10(4.23);
        % (Reaction 21)
        p37=log10(13413.194);
        % (Reaction 25)
        p45=log10(549.7896);
        % (Reaction 30)
        p55=log10(28.28);
                        
       
        % rate Parameters
        rateParameters=[valVector(14) valVector(2) valVector(3) valVector(13) valVector(12) valVector(11) p13 p16 p17 p20 p21 p1 p1 p23 p25 p27 p30 p31 p33 p36 p37 valVector(10) valVector(7) valVector(8) p45 valVector(4) valVector(5) valVector(6) valVector(9) p55 valVector(1)]';
        
        %rateParameters=[p3 p5 p6 p8 p9 p11 p13 p16 p17 p20 p21 p1 p1 p23 p25 p27 p30 p31 p33 p36 p37 p39 p41 p44 p45 p47 p49 p50 p52 p55 p56]';
        
        l=0.0;
        
               
        nData=size(data,2);
        
        numForwardModelSolves=numForwardModelSolves+1;
         
        
        % This is the call to the main reaction rate computation routine
                
        Rates=main(modelProblem,nReactions,nSpecies,nData,input,rateParameters);
               
        m=[];
                        
        for t=1:nData
            
             %modelPrediction=[1:nSpecies;Rates((t-1)*(nReactions+nSpecies)+nReactions+1:(t-1)*(nReactions+nSpecies)+nReactions+nSpecies)];
                             
             %modelPrediction=speciesConcns(2,8);
             
             modelPrediction=Rates((t-1)*(nReactions+nSpecies)+nReactions+8);             
                                       
             l=l-((data(t)-modelPrediction)^2)/2/dataVariance; 
             
             % This is a check to see if the ODE integrator was
             % successful, else we set loglikelihood to -infinity
             if modelPrediction==-1
                l=-1E100;
                break;
             end 
            
        end
             
        
 
    %%

    %% (model problem 2)
    elseif modelProblem==2
        
        k1=valVector(1);
        k2=valVector(2);
        k2r=valVector(3);
        k3=valVector(4);
        k4=valVector(5);
        k5=valVector(6);
        k6=valVector(7);
        k7=valVector(8);
        k8=valVector(9);
        k9=valVector(10);
        k10=valVector(11);
        k11=valVector(12);
        k12=valVector(13);
        k13=valVector(14);
        k14=valVector(15);
        k15=valVector(16);
        k16=valVector(17);
        k17=valVector(18);
        k18=valVector(19);
        k19=valVector(20);
        k20=valVector(21);
        k21=valVector(22);
        k22=valVector(23);
        k23=valVector(24);
        k24=valVector(25);
        k25=valVector(26);
        k26=valVector(27);
        k27=valVector(28);
        k28=valVector(29);
        k29=valVector(30);
        k30=valVector(31);        
                
        % rate Parameters
        rateParameters=[k1 k2 k2r k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 k19 k20 k21 k22 k23 k24 k25 k26 k27 k28 k29 k30]';        
        
        l=0.0;
        
        nData=size(data,2);
           
        %for t=1:s
            
            numForwardModelSolves=numForwardModelSolves+1;
            
            time=input;
            
            % This is the call to the main reaction rate computation routine
            
            Rates=main(modelProblem,nReactions,nSpecies,nData,time,rateParameters);
                
            for t=1:nData   
                
               % compute reaction rates for DRG reactions
               speciesConcns=[1:nSpecies;Rates((t-1)*(nReactions+nSpecies)+nReactions+1:(t-1)*(nReactions+nSpecies)+nReactions+nSpecies)];
        
               modelPrediction=speciesConcns(2,19);
                
               l=l-((data(t)-modelPrediction)^2)/2/dataVariance;
               
            end            
            
            
            
            
        %end
    %%    

    %% (model problem 3)
    elseif modelProblem==3
        
        k1=valVector(1);
        k2=valVector(2);
        k2r=valVector(3);
        k3=valVector(4);
        k4=valVector(5);
        k5=valVector(6);
        k6=valVector(7);
        k7=valVector(8);
        k8=valVector(9);
        k9=valVector(10);
        k10=valVector(11);
        k11=valVector(12);
        k12=valVector(13);
        k13=valVector(14);
        k14=valVector(15);
        k15=valVector(16);
        k16=valVector(17);
        k17=valVector(18);
                        
        % rate Parameters
        rateParameters=[k1 k2 k2r k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17]';        
        
        l=0.0;
        
        %{
        s=size(data,2);
           
        for t=1:s
            
            numForwardModelSolves=numForwardModelSolves+1;
            
            time=input(t);
            
            % This is the call to the main reaction rate computation routine
            
            % The variable Rates actually also has species concns in it
            Rates=main(modelProblem,nReactions,nSpecies,time,rateParameters);
           
            % compute reaction rates for DRG reactions
            speciesConcns=[1:nSpecies;Rates(nReactions+1:nReactions+nSpecies)];
            
            modelPrediction=speciesConcns(2,17);

            l=l-((data(t)-modelPrediction)^2)/2/dataVariance;
            
        end
        %}
        
        
         nData=size(data,2);
           
         numForwardModelSolves=numForwardModelSolves+1;
            
         time=input;
            
         % This is the call to the main reaction rate computation routine
            
         Rates=main(modelProblem,nReactions,nSpecies,nData,time,rateParameters);
                
         m=[];
         
         for t=1:nData   
                
             % compute reaction rates for DRG reactions
             speciesConcns=[1:nSpecies;Rates((t-1)*(nReactions+nSpecies)+nReactions+1:(t-1)*(nReactions+nSpecies)+nReactions+nSpecies)];
                             
             modelPrediction=speciesConcns(2,16);
             
             m=[m modelPrediction];
             
             l=l-((data(t)-modelPrediction)^2)/2/dataVariance;
                          
         end        
         
    %%      
    
    
    
    %% (model problem 4)
    elseif modelProblem==4
        
        % Older working version of loglikelihood   
      
        % rate Parameters
        rateParameters=[valVector(1) valVector(2) valVector(3) valVector(4) valVector(5) valVector(6) valVector(7) valVector(8) valVector(9) valVector(10)]';
        
        l=0.0;
                
        nData=size(data,2);
           
        numForwardModelSolves=numForwardModelSolves+1;
            
        % This is the call to the main reaction rate computation routine
                
        Rates=main(modelProblem,nReactions,nSpecies,nData,input,rateParameters);
                        
        for t=1:nData   
                
             % compute reaction rates for DRG reactions
             %speciesConcns=[1:nSpecies;Rates((t-1)*(nReactions+nSpecies)+nReactions+1:(t-1)*(nReactions+nSpecies)+nReactions+nSpecies)];
                                 
             %modelPrediction=speciesConcns(2,11);             
             
             modelPrediction=Rates((t-1)*(nReactions+nSpecies)+nReactions+11);
                
             l=l-((data(t)-modelPrediction)^2)/2/dataVariance;        
             
             % This is a check to see if the ODE integrator was
             % successfull, else we set loglikelihood to -infinity
             
             if modelPrediction==-1 
                l=-1E100;
                break;
             end   
                 
             
        end         
          

       
    
    %% (model problem 8)
    elseif modelProblem==8
       
        
       for i=1:2
           if valVector(i)==-100
               valVector(i)=0.0;
           end
       end
        
       numForwardModelSolves=numForwardModelSolves+1;
       
       
       modelPrediction=input*valVector;
       
       l=0.0;
       
       s=size(data,2);
       
       for t=1:s
            l=l-((data(t)-modelPrediction(t))^2)/2/dataVariance;  
       end    
          
    
    %% (model problem 9)
    elseif modelProblem==9
    
       % if standard reversible jump 
       if typeSampling==0 || typeSampling==2 || typeSampling==3 || typeSampling==4
          
          if valVector(1)==-100
              valVector(2)=-100;
              valVector(3)=-100;
          elseif valVector(2)==-100
              valVector(3)=-100;
          end
          
          if valVector(6)==-100
              valVector(5)=-100;
              valVector(4)=-100;
          elseif valVector(5)==-100
              valVector(4)=-100;
          end    
                    
       end    
        
       if valVector(3)==-100 && valVector(4)==-100
        for i=1:6
          valVector(i)=0.0; 
        end   
       else    
        for i=1:6
           if valVector(i)==-100
               valVector(i)=0.0;
           end
        end
       end 
       
       modelPrediction=input*valVector;
       
       l=0.0;

       
       numForwardModelSolves=numForwardModelSolves+1;
              
       s=size(data,2);
       
       for t=1:s
            l=l-((data(t)-modelPrediction(t))^2)/2/dataVariance;  
       end        
    
       
    
    %% (model problem 10)
    elseif modelProblem==10
         
        % if standard reversible jump
        if typeSampling==0 || typeSampling==2 || typeSampling==3 || typeSampling==4
            if valVector(1)==-100
                valVector(2)=-100;
            end
        end
                
        for i=1:3
            if valVector(i)==-100
                valVector(i)=0.0;
            end
        end
       
       modelPrediction=input*valVector;

       numForwardModelSolves=numForwardModelSolves+1;      
       
       l=0.0;
       
       s=size(data,2);
       
       for t=1:s                      
            l=l-((data(t)-modelPrediction(t))^2)/2/dataVariance;  
       end       
                
    end    

