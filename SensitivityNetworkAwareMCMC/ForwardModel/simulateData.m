function [inputX,data] = simulateData(modelProblem,nReactions,nSpecies)

rng(1,'twister');
rng(1,'twister');

tic

if modelProblem==2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model Problem 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dataStdDev=2.0;
    
    nData=30;
    time=linspace(0.1,5.0,nData);
    
    % Set parameter values
    
    k1=1.0;
    
    k2=1.0;
    k2r=1.0;
    
    k3=1.0;
    
    k4=1.0;
    
    k5=1.0;
    
    k6=1.0;
    
    k7=1.0;
    
    k8=1.0;
    
    k9=1.0;
    
    k10=1.0;
    
    k11=1.0;
    
    k12=1.0;
    
    k13=1.0;
    
    k14=1.0;
    
    k15=1.0;
    
    k16=1.0;
    
    k17=1.0;
    
    k18=1.0;
    
    k19=1.0;
    
    k20=1.0;
    
    k21=1.0;
    
    k22=1.0;
    
    k23=1.0;
    
    k24=1.0;
    
    k25=1.0;
    
    k26=1.0;
    
    k27=1.0;
    
    k28=1.0;
    
    k29=1.0;
    
    k30=1.0;
    
    % rate Parameters
    rateParameters=[k1 k2 k2r k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 k19 k20 k21 k22 k23 k24 k25 k26 k27 k28 k29 k30]';
    
    s=[];
        
    % The variable Rates actually also has species concns in it
    Rates=main(modelProblem,nReactions,nSpecies,nData,time,rateParameters);
    
    dataSpecies=[4 19];
    
     for t=1:nData   


        inputX(t)=time(t);             
                  
        % compute reaction rates for DRG reactions
        speciesConcns=[1:nSpecies;Rates((t-1)*(nReactions+nSpecies)+nReactions+1:(t-1)*(nReactions+nSpecies)+nReactions+nSpecies)];
        
        for j=1:length(dataSpecies)
            
            modelPrediction=speciesConcns(2,dataSpecies(j));
            
            s=[s speciesConcns(2,dataSpecies(j))'];
            
            data(j,t)=modelPrediction+randn(1)*dataStdDev;
            
        end
        
     end
     
     
     s 
     save inputXMP2nData30STD2MultipleOutputs.mat inputX    
     save dataMP2nData30STD2MultipleOutputs.mat data
   
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Model Problem 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Num reactions 13
 
 % Num species 15

 
 
 
    
%% (model problem 1)
if modelProblem==1
    
    dataStdDev=0.2;
 
    nData=15361;    
    
    time1=linspace(0.0,3,nData);
           
    %{
    time2=[time1 linspace(2.1,3.0,770)];
    
    time3=[time2 linspace(3.1,4,1000)];
    
    time4=[time3 linspace(4.1,5,2000)];
    
    time5=[time4 linspace(5.1,7,3000)];
    
    time6=[time5 linspace(7.1,9,3000)];    
    %}
    
    time=time1;
    
    r1=1;
    r2=1;
    r3=1;
    r4=1;
    r5=1;
    r6=1;
    r7=0;
    r8=1;
    r9=1;
    r10=1;
    r11=1;
    r12=1;
    r13=1;
    
    
    if r1==1
        p56=0.0;
    else
        p56=-100;
    end
    
    if r2==1
        p5=1.5;
        p6=0.0;
    else
        p5=-100;
        p6=-100;
    end
    
    if r3==1
       p47=0.5;
    else
        p47=-100;
    end
    
    if r4==1
        p49=2.0;
    else
        p49=-100;
    end
    
    if r5==1
        p50=2.0;
    else
        p50=-100;
    end
    
    if r6==1
        p41=0.4;
    else
        p41=-100;
    end
    
    if r7==1
        p44=3.0; 
    else
        p44=-100;
    end
    
    if r8==1
        p52=0.5;
    else
        p52=-100;
    end
    
    if r9==1
        p39=1.0;
    else
        p39=-100;
    end
    
    if r10==1
        p11=0.0;
    else
        p11=-100;
    end
    
    if r11==1
        p9=0.5;
    else
        p9=-100;
    end
    
    if r12==1
        p8=4.0;
    else
        p8=-100;
    end
    
    if r13==1
        p3=2.5;
    else
        p3=-100;
    end
                    
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
        rateParameters=[p3 p5 p6 p8 p9 p11 p13 p16 p17 p20 p21 p1 p1 p23 p25 p27 p30 p31 p33 p36 p37 p39 p41 p44 p45 p47 p49 p50 p52 p55 p56]';
        
        s=[];
                                
        % This is the call to the main reaction rate computation routine
        
        Rates=main(modelProblem,nReactions,nSpecies,nData,time,rateParameters);
               
        for t=1:nData
            
             % compute reaction rates for DRG reactions
             speciesConcns=[1:nSpecies;Rates((t-1)*(nReactions+nSpecies)+nReactions+1:(t-1)*(nReactions+nSpecies)+nReactions+nSpecies)];
                             
             modelPrediction=speciesConcns(2,8);
             
             s=[s speciesConcns(2,8)']; 
             
             inputX(t)=time(t);
             
             data(t)=modelPrediction+randn(1)*dataStdDev;
                         
        end
        
        s
                              
        save inputXMP1nData15360STD0andTwentyExcl7.mat inputX
        save dataMP1nData15360STD0andTwentyExcl7.mat data
        
end
    %%




if modelProblem==3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model Problem 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Data generated excluding reaction 13
    
    dataStdDev=2.0;
    
    nData=20;
    time=linspace(1.0,2.0,nData);
    
    % Set parameter values
    
    r1=1;
    r2=1;
    r3=1;
    r4=1;
    r5=1;
    r6=1;
    r7=1;
    r8=1;
    r9=1;
    r10=1;
    r11=1;
    r12=1;
    r13=1;
    r14=1;
    r15=1;
    r16=1;
    r17=1;
    
    
    if r1==1
     k1=1.0;
    else
     k1=-100;
    end 
    
    if r2==1
     k2=1.0;
     k2r=1.0;
    else
     k2=-100;
     k2r=-100;
    end  
    
    if r3==1
     k3=1.0;
    else
     k3=-100;
    end 
     
    if r4==1
     k4=0.4;
    else
     k4=-100;
    end 
    
    if r5==1
     k5=0.2;
    else
     k5=-100;
    end 

    if r6==1 
     k6=0.5;
    else
     k6=-100;
    end 

    if r7==1
     k7=0.1;
    else
     k7=-100;
    end
    
    if r8==1
     k8=0.5;
    else
     k8=-100;
    end 
     
    if r9==1
     k9=0.6;
    else
     k9=-100;
    end 

    if r10==1
     k10=0.3;
    else
     k10=-100;
    end 
     
    if r11==1
     k11=2.0;
    else
     k11=-100;   
    end
    
    if r12==1
     k12=1.0;
    else
     k12=-100;
    end 
     
    if r13==1
     k13=0.7;
    else
     k13=-100;
    end 

    if r14==1
     k14=0.3;
    else
     k14=-100;
    end 

    if r15==1
     k15=1.0;
    else
     k15=-100;
    end
     
    if r16==1
     k16=0.5;
    else
     k16=-100;
    end
    
    if r17==1
     k17=1.2;
    else
     k17=-100;
    end 
    
    % rate Parameters
    rateParameters=[k1 k2 k2r k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17]';
    
    %{
    rateParameters=[0.4174
    0.0058
    1.9991
    1.9991
    0.9945
    1.0021
    0.9942
    1.0109
    0.9281
    0.0520
    0.8262
    2.0000
 -100.0000
    0.8337
    0.8164
    0.7316
    1.8802
    1.8301];
    %}
    
    s=[];
        
    %for t=1:10
        
    % The variable Rates actually also has species concns in it
    Rates=main(modelProblem,nReactions,nSpecies,nData,time,rateParameters);
        
        
     for t=1:nData   
        % compute reaction rates for DRG reactions
        speciesConcns=[1:nSpecies;Rates((t-1)*(nReactions+nSpecies)+nReactions+1:(t-1)*(nReactions+nSpecies)+nReactions+nSpecies)];
        
        modelPrediction=speciesConcns(2,16);
        
        s=[s speciesConcns(2,16)'];
        
        inputX(t)=time(t);
        %data(t)=modelPrediction;
        data(t)=modelPrediction+randn(1)*dataStdDev;
     end
    
     s
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modelProblem 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if modelProblem==4
    
    % true model generated by removing reaction 5
    
    dataStdDev=0.5;
        
    nData=10;
    time=linspace(0.1,0.5,nData);
    
    % set true rate parameters
    
    k1=1.0;
    
    k2=2.0;
    
    k3=1.4;
    k3=-100;
    
    k4=1.0;
    
    k5=1.5;
    %k5=-100;
    
    k6=0.8;
    
    k7=1.0;
    
    k8=1.0;
    
    k9=0.5;
    
    k10=1.0;
        
    s=[(11)'];    
        
    % rate Parameters
    rateParameters=[k1 k2 k3 k4 k5 k6 k7 k8 k9 k10]';

    % The variable Rates actually also has species concns in it
    Rates=main(modelProblem,nReactions,nSpecies,nData,time,rateParameters);   
        
    for t=1:nData
        
       speciesConcns=[1:nSpecies;Rates((t-1)*(nReactions+nSpecies)+nReactions+1:(t-1)*(nReactions+nSpecies)+nReactions+nSpecies)];     

       modelPrediction=speciesConcns(2,11);
       
       s=[s speciesConcns(2,11)'];       
               
       inputX(t)=time(t);
       %data(t)=modelPrediction;
       data(t)=modelPrediction+randn(1)*dataStdDev;
        
    end
    
    s
    
    save s.mat s
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model problem 9

if modelProblem==9
    
%alpha=[0.0; 0.0; 0.0; 2.0; 1.0; 2.0];

alpha=[0.0; 0.0; 0.0; 0.5; 0.5; 0.5];

beta=0.5;

inputX =[-0.6490   -0.8519   -0.2752    1.2708    0.6033   -0.5581
          1.1812    0.8003    0.6037    0.0660   -0.5352   -0.0285
         -0.7585   -1.5094    1.7813    0.4513   -0.1551   -1.4763
         -1.1096    0.8759    1.7737   -0.3222    0.6121    0.2589
         -0.8456   -0.2428   -1.8651    0.7884   -1.0443   -2.0187
         -0.5727    0.1668   -1.0511    0.9287   -0.3456    0.1997
         -0.5587   -1.9654   -0.4174   -0.4908   -1.1714    0.4259
          0.1784   -1.2701    1.4022    1.7972   -0.6856   -1.2700
         -0.1969    1.1752   -1.3677    0.5907    0.9262   -0.4852
          0.5864    2.0292   -0.2925   -0.6358   -1.4817    0.5943];

data=inputX*alpha+1/sqrt(beta)*randn(10,1);  

data=data'    
    
    
end    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model problem 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if modelProblem==10

alpha=[2.0; 1.0; 2.0];

beta=0.5;

inputX=[1.0000    0.0685   -0.8299
        1.0000    0.9512   -0.1262
        1.0000   -0.3448    0.3087
        1.0000    0.0359   -0.5701
        1.0000    1.1221    0.5388
        1.0000    0.3237    0.4669
        1.0000   -0.4112   -1.9741
        1.0000   -0.6418    0.3630
        1.0000   -0.6615   -0.6617
        1.0000   -0.2253   -0.2770];

data=inputX*alpha+1/sqrt(beta)*randn(10,1);  

data=data'
    
end

toc
 
