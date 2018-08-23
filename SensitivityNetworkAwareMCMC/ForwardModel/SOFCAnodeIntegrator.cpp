#include "SOFCAnodeIntegrator.h"
#include "SOFCAnode.h"
#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>    // defines class IdealGasMix
#include <cantera/equilibrium.h>    // chemical equilibrium
#include <cantera/transport.h>      // transport properties
#include <cantera/numerics.h>
#include <cantera/integrators.h>
#include <cantera/kernel/Array.h>
#include <cantera/kernel/FuncEval.h>
//#include <ida/ida.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>


#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cantera/Edge.h>
#include <cantera/importPhase.h>
#include <cantera/Interface.h>

#include <ctime>
const double F=96485.0;
const double pi=22.0/7.0;

using namespace Cantera;
using namespace std;
using namespace Cantera_CXX;

SOFCAnodeIntegrator::SOFCAnodeIntegrator(int mP,int nR,int nS,double* rateParameters) // constructor
{
    
    m_time=0.0;
    m_init=false;
    m_nv=0;
    m_rtol=1.0e-11;
    
    m_rtolsens=1.0e-11;
    m_atols=1.0e-13;
    m_atolsens=1.0e-13;
    m_maxstep=-1.0;
    m_verbose=false;
    m_integ= newIntegrator("CVODE");
    
    m_integ->setMethod(BDF_Method);
    m_integ->setProblemType(DENSE+NOJAC);
    m_integ->setIterator(Newton_Iter);
    m_integ->setMaxSteps(1000000);
    
    modelProblem=mP;
    nReactions=nR;
    neqvar=nS;
    
    delta=-100;
    
    reactionYDT.resize(nReactions);
    speciesYDT.resize(neqvar);

    if (modelProblem==3)
    {    
        // Reaction 1:
        if (rateParameters[0]==delta)
         k1=0.0;
        else
         k1=pow(10,rateParameters[0]);
                
        // Reaction 2:
        if (rateParameters[1]==delta)
         k2=0.0;
        else
         k2=pow(10,rateParameters[1]);   
        
        // Reaction 2 reverse
        if (rateParameters[2]==delta)
         k2r=0.0;
        else
         k2r=pow(10,rateParameters[2]);   
                
        // Reaction 3:
        if (rateParameters[3]==delta)
         k3=0.0;
        else
         k3=pow(10,rateParameters[3]); 
        
        // Reaction 4:        
        if (rateParameters[4]==delta)
         k4=0.0;
        else
         k4=pow(10,rateParameters[4]);
        
        // Reaction 5:        
        if (rateParameters[5]==delta)
         k5=0.0;
        else
         k5=pow(10,rateParameters[5]);
        
        // Reaction 6:
        
        if (rateParameters[6]==delta)
         k6=0.0;
        else
         k6=pow(10,rateParameters[6]);
        
        // Reaction 7:
        
        if (rateParameters[7]==delta) 
         k7=0.0;
        else
         k7=pow(10,rateParameters[7]);
        
        // Reaction 8:
        
        if (rateParameters[8]==delta)
         k8=0.0;
        else
         k8=pow(10,rateParameters[8]);   
        
        // Reaction 9:
        
        if (rateParameters[9]==delta) 
         k9=0.0;
        else
         k9=pow(10,rateParameters[9]);
        
        // Reaction 10:
        
        if (rateParameters[10]==delta)
         k10=0.0;
        else
         k10=pow(10,rateParameters[10]);
        
        // Reaction 11:
        
        if (rateParameters[11]==delta) 
         k11=0.0;
        else
         k11=pow(10,rateParameters[11]);
        
        // Reaction 12:
        
        if (rateParameters[12]==delta)
         k12=0.0;
        else
         k12=pow(10,rateParameters[12]);
        
        // Reaction 13:
        if (rateParameters[13]==delta) 
         k13=0.0;
        else
         k13=pow(10,rateParameters[13]);   
        
        // Reaction 14:
        
        if (rateParameters[14]==delta)
         k14=0.0;
        else
         k14=pow(10,rateParameters[14]); 
            
        // Reaction 15:
        
        if (rateParameters[15]==delta) 
         k15=0.0;
        else
         k15=pow(10,rateParameters[15]);
        
        // Reaction 16:
        
        if (rateParameters[16]==delta)
         k16=0.0;
        else
         k16=pow(10,rateParameters[16]);
        
        // Reaction 17:
        
        if (rateParameters[17]==delta)
         k17=0.0;
        else
         k17=pow(10,rateParameters[17]);   
               
        k1p=1.0;
        k3p=1.0;   
        k4p=1.0;          
        k5p=1.0;        
        k6p=1.0;        
        k7p=1.0;   
        k8p=1.0; 
        k9p=1.0;     
        k10p=1.0;        
        
        k11p=1.0;
        k12p=1.0;   
        k13p=1.0;          
        k14p=1.0;        
        k15p=1.0;        
        k16p=1.0;   
        k17p=1.0; 
        k18p=1.0;     
        k19p=1.0;  
        k20p=1.0;
        
        k21p=1.0;
        k22p=1.0;   
        k23p=1.0;          
        k24p=1.0;        
        k25p=1.0;        
        k27p=1.0;   
        k28p=1.0; 
        k29p=1.0;    
                    
    }       
    
   
    if (modelProblem==2)
    {    
        // Reaction 1:
        if (rateParameters[0]==delta)
         k1=0.0;
        else
         k1=pow(10,rateParameters[0]);
                
        // Reaction 2:
        if (rateParameters[1]==delta)
         k2=0.0;
        else
         k2=pow(10,rateParameters[1]);   
        
        // Reaction 2 reverse
        if (rateParameters[2]==delta)
         k2r=0.0;
        else
         k2r=pow(10,rateParameters[2]);   
                
        // Reaction 3:
        if (rateParameters[3]==delta)
         k3=0.0;
        else
         k3=pow(10,rateParameters[3]); 
        
        // Reaction 4:        
        if (rateParameters[4]==delta)
         k4=0.0;
        else
         k4=pow(10,rateParameters[4]);
        
        // Reaction 5:        
        if (rateParameters[5]==delta)
         k5=0.0;
        else
         k5=pow(10,rateParameters[5]);
        
        // Reaction 6:
        
        if (rateParameters[6]==delta)
         k6=0.0;
        else
         k6=pow(10,rateParameters[6]);
        
        // Reaction 7:
        
        if (rateParameters[7]==delta) 
         k7=0.0;
        else
         k7=pow(10,rateParameters[7]);
        
        // Reaction 8:
        
        if (rateParameters[8]==delta)
         k8=0.0;
        else
         k8=pow(10,rateParameters[8]);   
        
        // Reaction 9:
        
        if (rateParameters[9]==delta) 
         k9=0.0;
        else
         k9=pow(10,rateParameters[9]);
        
        // Reaction 10:
        
        if (rateParameters[10]==delta)
         k10=0.0;
        else
         k10=pow(10,rateParameters[10]);
        
        // Reaction 11:
        
        if (rateParameters[11]==delta) 
         k11=0.0;
        else
         k11=pow(10,rateParameters[11]);
        
        // Reaction 12:
        
        if (rateParameters[12]==delta)
         k12=0.0;
        else
         k12=pow(10,rateParameters[12]);
        
        // Reaction 13:
        if (rateParameters[13]==delta) 
         k13=0.0;
        else
         k13=pow(10,rateParameters[13]);   
        
        // Reaction 14:
        
        if (rateParameters[14]==delta)
         k14=0.0;
        else
         k14=pow(10,rateParameters[14]); 
            
        // Reaction 15:
        
        if (rateParameters[15]==delta) 
         k15=0.0;
        else
         k15=pow(10,rateParameters[15]);
        
        // Reaction 16:
        
        if (rateParameters[16]==delta)
         k16=0.0;
        else
         k16=pow(10,rateParameters[16]);
        
        // Reaction 17:
        
        if (rateParameters[17]==delta)
         k17=0.0;
        else
         k17=pow(10,rateParameters[17]);   
        
        // Reaction 18:
        
        if (rateParameters[18]==delta) 
         k18=0.0;
        else
         k18=pow(10,rateParameters[18]);   
        
        // Reaction 19:
        
        if (rateParameters[19]==delta) 
         k19=0.0;
        else
         k19=pow(10,rateParameters[19]);   
        
        // Reaction 20:
        
        if (rateParameters[20]==delta)
         k20=0.0;
        else
         k20=pow(10,rateParameters[20]);   
        
        // Reaction 21:
        
        if (rateParameters[21]==delta)
         k21=0.0;
        else
         k21=pow(10,rateParameters[21]);
        
        // Reaction 22:
        
        if (rateParameters[22]==delta) 
         k22=0.0;
        else
         k22=pow(10,rateParameters[22]);   
        
        // Reaction 23:
        
        if (rateParameters[23]==delta)
         k23=0.0;
        else
         k23=pow(10,rateParameters[23]);   
        
        // Reaction 24:
        
        if (rateParameters[24]==delta) 
         k24=0.0;
        else
         k24=pow(10,rateParameters[24]);   
        
        // Reaction 25:
        if (rateParameters[25]==delta) 
         k25=0.0;
        else
         k25=pow(10,rateParameters[25]);   
        
        // Reaction 26:
        if (rateParameters[26]==delta)
         k26=0.0;
        else
         k26=pow(10,rateParameters[26]);   
        
        // Reaction 27:
        if (rateParameters[27]==delta)
         k27=0.0;
        else
         k27=pow(10,rateParameters[27]);   
        
        // Reaction 28:
        
        if (rateParameters[28]==delta)
         k28=0.0;
        else
         k28=pow(10,rateParameters[28]);   
        
        // Reaction 29:
        if (rateParameters[29]==delta)
         k29=0.0;
        else
         k29=pow(10,rateParameters[29]);   
        
        // Reaction 30:
        
        if (rateParameters[30]==delta)
         k30=0.0;
        else
         k30=pow(10,rateParameters[30]);   
        

        k1p=1.0;
        k3p=1.0;   
        k4p=1.0;          
        k5p=1.0;        
        k6p=1.0;        
        k7p=1.0;   
        k8p=1.0; 
        k9p=1.0;     
        k10p=1.0;        
        
        k11p=1.0;
        k12p=1.0;   
        k13p=1.0;          
        k14p=1.0;        
        k15p=1.0;        
        k16p=1.0;   
        k17p=1.0; 
        k18p=1.0;     
        k19p=1.0;  
        k20p=1.0;
        
        k21p=1.0;
        k22p=1.0;   
        k23p=1.0;          
        k24p=1.0;        
        k25p=1.0;        
        k27p=1.0;   
        k28p=1.0; 
        k29p=1.0;    
                    
    }    
    if (modelProblem==1)
    {

        // Reaction 1:
        if (rateParameters[0]==delta)
         p3=0.0;
        else
         p3=pow(10,rateParameters[0]);
                
        // Reaction 2:
        if (rateParameters[1]==delta)
         p5=0.0;
        else
         p5=pow(10,rateParameters[1]);   
        
        // Reaction 3:
        if (rateParameters[2]==delta)
         p6=0.0;
        else
         p6=pow(10,rateParameters[2]); 
        
        // Reaction 4:
        
        if (rateParameters[3]==delta)
         p8=0.0;
        else
         p8=pow(10,rateParameters[3]);
        
        // Reaction 5:
        
        if (rateParameters[4]==delta)
         p9=0.0;
        else
         p9=pow(10,rateParameters[4]);
        
        // Reaction 6:
        
        if (rateParameters[5]==delta)
         p11=0.0;
        else
         p11=pow(10,rateParameters[5]);
        
        // Reaction 7:
        
        if (rateParameters[6]==delta) 
         p13=0.0;
        else
         p13=pow(10,rateParameters[6]);
        
        // Reaction 8:
        
        if (rateParameters[7]==delta)
         p16=0.0;
        else
         p16=pow(10,rateParameters[7]);   
        
        // Reaction 9:
        
        if (rateParameters[8]==delta) 
         p17=0.0;
        else
         p17=pow(10,rateParameters[8]);
        
        // Reaction 10:
        
        if (rateParameters[9]==delta)
         p20=0.0;
        else
         p20=pow(10,rateParameters[9]);
        
        // Reaction 11:
        
        if (rateParameters[10]==delta) 
         p21=0.0;
        else
         p21=pow(10,rateParameters[10]);
        
        // Reaction 12:
        
        if (rateParameters[11]==delta)
         p1=0.0;
        else
         p1=pow(10,rateParameters[11]);
        
        // Reaction 13:
        if (rateParameters[12]==delta) 
         p1=0.0;
        else
         p1=pow(10,rateParameters[12]);   
        
        // Reaction 14:
        
        if (rateParameters[13]==delta)
         p23=0.0;
        else
         p23=pow(10,rateParameters[13]); 
            
        // Reaction 15:
        
        if (rateParameters[14]==delta) 
         p25=0.0;
        else
         p25=pow(10,rateParameters[14]);
        
        // Reaction 16:
        
        if (rateParameters[15]==delta)
         p27=0.0;
        else
         p27=pow(10,rateParameters[15]);
        
        // Reaction 17:
        
        if (rateParameters[16]==delta)
         p30=0.0;
        else
         p30=pow(10,rateParameters[16]);   
        
        // Reaction 18:
        
        if (rateParameters[17]==delta) 
         p31=0.0;
        else
         p31=pow(10,rateParameters[17]);   
        
        // Reaction 19:
        
        if (rateParameters[18]==delta) 
         p33=0.0;
        else
         p33=pow(10,rateParameters[18]);   
        
        // Reaction 20:
        
        if (rateParameters[19]==delta)
         p36=0.0;
        else
         p36=pow(10,rateParameters[19]);   
        
        // Reaction 21:
        
        if (rateParameters[20]==delta)
         p37=0.0;
        else
         p37=pow(10,rateParameters[20]);
        
        // Reaction 22:
        
        if (rateParameters[21]==delta) 
         p39=0.0;
        else
         p39=pow(10,rateParameters[21]);   
        
        // Reaction 23:
        
        if (rateParameters[22]==delta)
         p41=0.0;
        else
         p41=pow(10,rateParameters[22]);   
        
        // Reaction 24:
        
        if (rateParameters[23]==delta) 
         p44=0.0;
        else
         p44=pow(10,rateParameters[23]);   
        
        // Reaction 25:
        if (rateParameters[24]==delta) 
         p45=0.0;
        else
         p45=pow(10,rateParameters[24]);   
        
        // Reaction 26:
        if (rateParameters[25]==delta)
         p47=0.0;
        else
         p47=pow(10,rateParameters[25]);   
        
        // Reaction 27:
        if (rateParameters[26]==delta)
         p49=0.0;
        else
         p49=pow(10,rateParameters[26]);   
        
        // Reaction 28:
        
        if (rateParameters[27]==delta)
         p50=0.0;
        else
         p50=pow(10,rateParameters[27]);   
        
        // Reaction 29:
        if (rateParameters[28]==delta)
         p52=0.0;
        else
         p52=pow(10,rateParameters[28]);   
        
        // Reaction 30:
        
        if (rateParameters[29]==delta)
         p55=0.0;
        else
         p55=pow(10,rateParameters[29]);   
        
        // Reaction 31:
        
        if (rateParameters[30]==delta)
         p56=0.0;
        else
         p56=pow(10,rateParameters[30]);   
        
        
        p2=1.0;
        p4=8176.5664;
        p7=9834.13;
        p10=13.73;
        p12=12457.816;
        p14=1961.79;
        p15=5930.78;
        p18=9815.0778;
        p19=712.53;
        p22=378.406;
        p24=0;
        p26=3767.4921;
        p28=4395.21;
        p29=15124.941;
        p32=6000.7723;
        p34=1130.76;
        p35=2251.46;
        p38=994.69;
        p40=6808.32;
        p42=17991.179;
        p43=5507.89;
        p46=4.02;
        p48=3386.3875;
        p51=3566;
        p53=7631.63;
        p54=5249;
         
         
         
        /*
        p48=100;
        p40=100;
        p51=100;
        p42=100;
        p43=100;
        p53=100;
        p10=100;
        p12=100;
        p4=100;
        p7=100;
        */
        
        
    }
    if (modelProblem==4)
    {
        
        // Reaction 1
        if (rateParameters[0]==delta)
            k1=0.0;
        else
            k1=pow(10,rateParameters[0]);
                
        // Reaction 2
        if (rateParameters[1]==delta)
            k2=0.0;
        else
            k2=pow(10,rateParameters[1]);
                
        // Reaction 3
        if (rateParameters[2]==delta) 
            k3=0.0;
        else
            k3=pow(10,rateParameters[2]);
        
        
        //Reaction 4
        if (rateParameters[3]==delta)
            k4=0.0;
        else
            k4=pow(10,rateParameters[3]);
        
        if (rateParameters[4]==delta)
            k5=0.0;
        else
            k5=pow(10,rateParameters[4]);
        
        if (rateParameters[5]==delta)
            k6=0.0;
        else
            k6=pow(10,rateParameters[5]);
        
        if (rateParameters[6]==delta)
            k7=0.0;
        else
            k7=pow(10,rateParameters[6]);

        if (rateParameters[7]==delta)
            k8=0.0;
        else
            k8=pow(10,rateParameters[7]);
        
        if (rateParameters[8]==delta)
            k9=0.0;
        else
            k9=pow(10,rateParameters[8]);
        
        if (rateParameters[9]==delta)
            k10=0.0;
        else
            k10=pow(10,rateParameters[9]);        
    }
    if (modelProblem==5)
    {
        
        // Reaction 1 forward
        if (rateParameters[0]==delta)
            modelProblem5k1f=0.0;
        else
            modelProblem5k1f=pow(10,rateParameters[0]);
                
        // Reaction 1 reverse
        if (rateParameters[1]==delta)
            modelProblem5k1r=0.0;
        else
            modelProblem5k1r=pow(10,rateParameters[1]);
        
        // Reaction 2
        if (rateParameters[2]==delta)
            modelProblem5k2=0.0;
        else
            modelProblem5k2=pow(10,rateParameters[2]);
        
        // Reaction 3
        if (rateParameters[3]==delta) 
            modelProblem5k3=0.0;
        else
            modelProblem5k3=pow(10,rateParameters[3]);
        
        // Reaction 4
        if (rateParameters[4]==delta)
            modelProblem5k4=0.0;
        else
            modelProblem5k4=pow(10,rateParameters[4]);
        
        // Reaction 5
        if (rateParameters[5]==delta)
            modelProblem5k5=0.0;
        else
            modelProblem5k5=pow(10,rateParameters[5]);
        
        // Reaction 6
        if (rateParameters[6]==delta)
            modelProblem5k6=0.0;
        else
            modelProblem5k6=pow(10,rateParameters[6]);
    }    
    if (modelProblem==6)
    {
        
        // Reaction 1 forward
        if (rateParameters[0]==delta)
            modelProblem6k1f=0.0;
        else
            modelProblem6k1f=pow(10,rateParameters[0]);
                
        // Reaction 1 reverse
        if (rateParameters[1]==delta)
            modelProblem6k1r=0.0;
        else
            modelProblem6k1r=pow(10,rateParameters[1]);
        
        // Reaction 2
        if (rateParameters[2]==delta)
            modelProblem6k2=0.0;
        else
            modelProblem6k2=pow(10,rateParameters[2]);
        
        // Reaction 3
        if (rateParameters[3]==delta) 
            modelProblem6k3=0.0;
        else
            modelProblem6k3=pow(10,rateParameters[3]);
        
        // Reaction 4
        if (rateParameters[4]==delta)
            modelProblem6k4=0.0;
        else
            modelProblem6k4=pow(10,rateParameters[4]);
        
        // Reaction 5
        if (rateParameters[5]==delta)
            modelProblem6k5=0.0;
        else
            modelProblem6k5=pow(10,rateParameters[5]);
        
    }      
    
}

SOFCAnodeIntegrator::~SOFCAnodeIntegrator()  // destructor
{
    deleteIntegrator(m_integ);
}

int SOFCAnodeIntegrator::neq()
{
    return neqvar;
}

void SOFCAnodeIntegrator::eval(doublereal t,doublereal* y,doublereal* ydot,doublereal* p)
{

    if (modelProblem==3)
    {
        double uBEGFR=y[0];
        double remcRaf=y[1];
        double remSOS=y[2];
        double inactiveSOS=y[3];
        double inactiveRas=y[4];
        double inactivePKA=y[5];
        double cRafPP=y[6];
        double cRaf=y[7];
        double bEGFR=y[8];
        double activeSOS=y[9];
        double activeRas=y[10];
        double PKA=y[11];
        double MEKPP=y[12];
        double MEK=y[13];
        double ERKPP=y[14];
        double ERK=y[15];
        double EGF=y[16];
        double Gap=y[17];
        double PKAA=y[18];
        double Cilostamide=y[19];
        
        ydot[0]=-k2*EGF*uBEGFR+k2r*bEGFR;
                
        ydot[1]=k13*PKA*cRaf/(k13p+cRaf);
                
        ydot[2]=k11*ERKPP*inactiveSOS/(k11p+inactiveSOS)+k12*ERKPP*activeSOS/(k12p+activeSOS);        
        
        ydot[3]=-k1*bEGFR*inactiveSOS/(k1p+inactiveSOS)+k3*activeSOS/(k3p+activeSOS)-k11*ERKPP*inactiveSOS/(k11p+inactiveSOS);
                
        ydot[4]=-k4*activeSOS*inactiveRas/(k4p+inactiveRas)+k5*Gap*activeRas/(k5p+activeRas);                
                
        ydot[5]=-k14*PKAA*inactivePKA/(k14p+inactivePKA)-k15*Cilostamide*inactivePKA/(k15p+inactivePKA)+k16*PKA/(k16p+PKA);                
                
        ydot[6]=k6*activeRas*cRaf/(k6p+cRaf)-k7*cRafPP/(k7p+cRafPP);
                
        ydot[7]=-k6*activeRas*cRaf/(k6p+cRaf)+k7*cRafPP/(k7p+cRafPP)-k13*PKA*cRaf/(k13p+cRaf);
                
        ydot[8]=k2*EGF*uBEGFR-k2r*bEGFR;     
                
        ydot[9]=k1*bEGFR*inactiveSOS/(k1p+inactiveSOS)-k3*activeSOS/(k3p+activeSOS)-k12*ERKPP*activeSOS/(k12p+activeSOS);
                
        ydot[10]=k4*activeSOS*inactiveRas/(k4p+inactiveRas)-k5*Gap*activeRas/(k5p+activeRas);        
                
        ydot[11]=k14*PKAA*inactivePKA/(k14p+inactivePKA)+k15*Cilostamide*inactivePKA/(k15p+inactivePKA)-k16*PKA/(k16p+PKA);
                
        ydot[12]=k8*cRafPP*MEK/(k8p+MEK)-k9*MEKPP/(k9p+MEKPP);
                
        ydot[13]=-k8*cRafPP*MEK/(k8p+MEK)+k9*MEKPP/(k9p+MEKPP);
                
        ydot[14]=k10*MEKPP*ERK/(k10p+ERK)-k17*ERKPP/(k17p+ERKPP);
                
        ydot[15]=-k10*MEKPP*ERK/(k10p+ERK)+k17*ERKPP/(k17p+ERKPP);
                                
        ydot[16]=-k2*EGF*uBEGFR+k2r*bEGFR;  
                                
        ydot[17]=0;
                
        ydot[18]=0;
                        
        ydot[19]=0;                
                             
    }    
    
    
    if (modelProblem==2)
    {
        double uBEGFR=y[0];
        double remcRaf=y[1];
        double remSOS=y[2];
        double inactiveSOS=y[3];
        double inactiveRas=y[4];
        double inactiveRap1=y[5];
        double inactivePKA=y[6];
        double inactiveEPAC=y[7];
        double cRafPP=y[8];
        double cRaf=y[9];
        double bEGFR=y[10];
        double activeSOS=y[11];
        double activeRas=y[12];
        double activeRap1=y[13];
        double PKA=y[14];
        double MEKPP=y[15];
        double MEK=y[16];
        double ERKPP=y[17];
        double ERK=y[18];
        double EPAC=y[19];
        double EGF=y[20];
        double BRafPP=y[21];
        double BRaf=y[22];
        double activeC3G=y[23];
        double inactiveC3G=y[24];
        double dEGFR=y[25];
        double EPACA=y[26];
        double Gap=y[27];
        double PKAA=y[28];
        double Cilostamide=y[29];
        
        ydot[0]=-k2*EGF*uBEGFR+k2r*bEGFR;
                
        ydot[1]=k13*PKA*cRaf/(k13p+cRaf);
                
        ydot[2]=k11*ERKPP*inactiveSOS/(k11p+inactiveSOS)+k12*ERKPP*activeSOS/(k12p+activeSOS);        
        
        ydot[3]=-k1*bEGFR*inactiveSOS/(k1p+inactiveSOS)+k3*activeSOS/(k3p+activeSOS)-k11*ERKPP*inactiveSOS/(k11p+inactiveSOS);
                
        ydot[4]=-k4*activeSOS*inactiveRas/(k4p+inactiveRas)+k5*Gap*activeRas/(k5p+activeRas);
                
        ydot[5]=-k20*EPAC*inactiveRap1/(k20p+inactiveRap1)+k21*Gap*activeRap1/(k21p+activeRap1)-k27*activeC3G*inactiveRap1/(k27p+inactiveRap1);
                
        ydot[6]=-k14*PKAA*inactivePKA/(k14p+inactivePKA)-k15*Cilostamide*inactivePKA/(k15p+inactivePKA)+k16*PKA/(k16p+PKA);
                
        ydot[7]=-k17*EPACA*inactiveEPAC/(k17p+inactiveEPAC)-k18*Cilostamide*inactiveEPAC/(k18p+inactiveEPAC)+k19*EPAC/(k19p+EPAC);
                
        ydot[8]=k6*activeRas*cRaf/(k6p+cRaf)-k7*cRafPP/(k7p+cRafPP);
                
        ydot[9]=-k6*activeRas*cRaf/(k6p+cRaf)+k7*cRafPP/(k7p+cRafPP)-k13*PKA*cRaf/(k13p+cRaf);
                
        ydot[10]=k2*EGF*uBEGFR-k2r*bEGFR-k30*bEGFR;     
                
        ydot[11]=k1*bEGFR*inactiveSOS/(k1p+inactiveSOS)-k3*activeSOS/(k3p+activeSOS)-k12*ERKPP*activeSOS/(k12p+activeSOS);
                
        ydot[12]=k4*activeSOS*inactiveRas/(k4p+inactiveRas)-k5*Gap*activeRas/(k5p+activeRas);        
        
        ydot[13]=k20*EPAC*inactiveRap1/(k20p+inactiveRap1)-k21*Gap*activeRap1/(k21p+activeRap1)+k27*activeC3G*inactiveRap1/(k27p+inactiveRap1);
                
        ydot[14]=k14*PKAA*inactivePKA/(k14p+inactivePKA)+k15*Cilostamide*inactivePKA/(k15p+inactivePKA)-k16*PKA/(k16p+PKA);
                
        ydot[15]=k8*cRafPP*MEK/(k8p+MEK)-k9*MEKPP/(k9p+MEKPP)+k24*BRafPP*MEK/(k24p+MEK);
                
        ydot[16]=-k8*cRafPP*MEK/(k8p+MEK)+k9*MEKPP/(k9p+MEKPP)-k24*BRafPP*MEK/(k24p+MEK);
                
        ydot[17]=k10*MEKPP*ERK/(k10p+ERK)-k29*ERKPP/(k29p+ERKPP);
                
        ydot[18]=-k10*MEKPP*ERK/(k10p+ERK)+k29*ERKPP/(k29p+ERKPP);
                
        ydot[19]=k17*EPACA*inactiveEPAC/(k17p+inactiveEPAC)+k18*Cilostamide*inactiveEPAC/(k18p+inactiveEPAC)-k19*EPAC/(k19p+EPAC);
                
        ydot[20]=-k2*EGF*uBEGFR+k2r*bEGFR;  
                
        ydot[21]=k22*activeRap1*BRaf/(k22p+BRaf)-k23*BRafPP/(k23p+BRafPP)+k28*activeRas*BRaf/(k28p+BRaf);
                
        ydot[22]=-k22*activeRap1*BRaf/(k22p+BRaf)+k23*BRafPP/(k23p+BRafPP)-k28*activeRas*BRaf/(k28p+BRaf);
                
        ydot[23]=k25*bEGFR*inactiveC3G/(k25p+inactiveC3G)-k26*activeC3G;        
        
        ydot[24]=-k25*bEGFR*inactiveC3G/(k25p+inactiveC3G)+k26*activeC3G;
                
        ydot[25]=k30*bEGFR;
                
        ydot[26]=0;
                
        ydot[27]=0;
                        
        ydot[28]=0;
                
        ydot[29]=0;
                             
    }
    if (modelProblem==1)
    {       

        double BEGFR=y[0];
        double DEGFR=y[1];
        double EGF=y[2];
        double UEGFR=y[3];
        double inactiveC3G=y[4];
        double activeC3G=y[5];
        double activeRap1=y[6];
        double BRaf=y[7];
        double BRafPP=y[8];
        double activeRas=y[9];
        double Gap=y[10];
        double inactiveRap1=y[11];
        double inactiveRas=y[12];
        double activeSOS=y[13];
        double inactiveSOS=y[14];
        
        ydot[0]=p5*EGF*UEGFR-p6*BEGFR-p56*BEGFR;
        
        ydot[1]=p56*BEGFR;
        
        ydot[2]=-p5*EGF*UEGFR+p6*BEGFR;
        
        ydot[3]=-p5*EGF*UEGFR+p6*BEGFR;
        
        ydot[4]=-p47*BEGFR*inactiveC3G/(p48+inactiveC3G)+p49*activeC3G;
        
        ydot[5]=p47*BEGFR*inactiveC3G/(p48+inactiveC3G)-p49*activeC3G;
        
        ydot[6]=-p39*Gap*activeRap1/(p40+activeRap1)+p50*activeC3G*inactiveRap1/(p51+inactiveRap1);
        
        ydot[7]=-p41*activeRap1*BRaf/(p42+BRaf)+p44*BRafPP/(p43+BRafPP)-p52*activeRas*BRaf/(p53+BRaf);
        
        ydot[8]=p41*activeRap1*BRaf/(p42+BRaf)-p44*BRafPP/(p43+BRafPP)+p52*activeRas*BRaf/(p53+BRaf);
        
        ydot[9]=p9*activeSOS*inactiveRas/(p10+inactiveRas)-p11*Gap*activeRas/(p12+activeRas);
        
        ydot[10]=0;
        
        ydot[11]=-p50*activeC3G*inactiveRap1/(p51+inactiveRap1)+p39*Gap*activeRap1/(p40+activeRap1);
        
        ydot[12]=-p9*activeSOS*inactiveRas/(p10+inactiveRas)+p11*Gap*activeRas/(p12+activeRas);
        
        ydot[13]=p3*BEGFR*inactiveSOS/(p4+inactiveSOS)-p8*activeSOS/(p7+activeSOS);
        
        ydot[14]=-p3*BEGFR*inactiveSOS/(p4+inactiveSOS)+p8*activeSOS/(p7+activeSOS);


    }
    if (modelProblem==4)
    {
        
        
        //ofstream myfile;
        //myfile.open ("outputs.txt",ios::app);
        //myfile<<scientific<<setprecision(16)<<y[0]<<" "<<y[1]<<" "<<y[2]<<" "<<y[3]<<" "<<y[4]<<" "<<y[5]<<" "<<y[6]<<" "<<y[7]<<" "<<y[8]<<" "<<y[9]<<" "<<y[10]<<" "<<"\n\n";
        //myfile.close();
        
        
        double A=y[0];
        double B=y[1];
        double C=y[2];
        double D=y[3];
        double E=y[4];
        double F=y[5];
        double G=y[6];
        double H=y[7];
        double I=y[8];
        double J=y[9];
        double K=y[10];
        
        ydot[0]=-k1*A-k3*A-k5*A+k2*C*D;
        ydot[1]=k1*A;
        ydot[2]=-k2*C*D;
        ydot[3]=-k2*C*D;
        ydot[4]=k3*A-k4*E*F;
        ydot[5]=-k4*E*F-k8*F*J;
        ydot[6]=k4*E*F-k10*G;
        ydot[7]=k5*A-k6*H-k7*H;
        ydot[8]=k6*H-k9*I;
        ydot[9]=k7*H-k8*F*J;
        ydot[10]=k8*F*J+k9*I+k10*G;
            
    }
    if (modelProblem==5)
    {
       double A=y[0];
       double B=y[1];
       double C=y[2];
       double D=y[3];
       double E=y[4];
       double F=y[5];
       double G=y[6];
       double H=y[7];
       double I=y[8];
       double J=y[9];
       
       ydot[0]=-modelProblem5k1r*A+modelProblem5k1f*B*C-modelProblem5k2*A-modelProblem5k5*A;
       
       ydot[1]=modelProblem5k1r*A-modelProblem5k1f*B*C;
       
       ydot[2]=modelProblem5k1r*A-modelProblem5k1f*B*C;
       
       ydot[3]=modelProblem5k2*A-modelProblem5k3*D*F;
       
       ydot[4]=modelProblem5k5*A-modelProblem5k6*E*G;
       
       ydot[5]=-modelProblem5k3*D*F-modelProblem5k4*F*G;
       
       ydot[6]=-modelProblem5k4*F*G-modelProblem5k6*E*G;
       
       ydot[7]=modelProblem5k6*E*G;
       
       ydot[8]=modelProblem5k4*F*G;
       
       ydot[9]=modelProblem5k3*D*F;
         
    }
    if (modelProblem==6)
    {
      double A=y[0];
      double B=y[1];
      double C=y[2];
      double D=y[3];
      double E=y[4];
      double F=y[5];
      
      ydot[0]=modelProblem6k1f*B*C-modelProblem6k1r*A-modelProblem6k2*A-modelProblem6k3*A;
      
      ydot[1]=-modelProblem6k1f*B*C+modelProblem6k1r*A;
      
      ydot[2]=-modelProblem6k1f*B*C+modelProblem6k1r*A;
      
      ydot[3]=modelProblem6k2*A-modelProblem6k4*D;
      
      ydot[4]=modelProblem6k3*A-modelProblem6k5*E;
      
      ydot[5]=modelProblem6k4*D+modelProblem6k5*E;
    }
    
    
    
}

void SOFCAnodeIntegrator::getInitialConditions(doublereal t0,size_t leny,doublereal* y)
{
    m_time=t0;

     if (modelProblem==3)
    {
        y[0]=500;
        y[1]=0;
        y[2]=0;
        y[3]=1200;
        y[4]=1200;
        y[5]=1000;
        y[6]=0;
        y[7]=1500;
        y[8]=0;
        y[9]=0;
        y[10]=0;
        y[11]=0;
        y[12]=0;
        y[13]=3000;
        y[14]=0;
        y[15]=10000;
        y[16]=1000;
        y[17]=2400;
        y[18]=10;
        y[19]=10;
    }   
    
    
    
    if (modelProblem==2)
    {
        y[0]=500;
        y[1]=0;
        y[2]=0;
        y[3]=1200;
        y[4]=1200;
        y[5]=1200;
        y[6]=1000;
        y[7]=1000;
        y[8]=0;
        y[9]=1500;
        y[10]=0;
        y[11]=0;
        y[12]=0;
        y[13]=0;
        y[14]=0;
        y[15]=0;
        y[16]=3000;
        y[17]=0;
        y[18]=10000;
        y[19]=0;
        y[20]=1000;
        y[21]=0;
        y[22]=1500;
        y[23]=0;
        y[24]=1200;
        y[25]=0;
        y[26]=1000;
        y[27]=2400;
        y[28]=1000;
        y[29]=1000;
    }
    if (modelProblem==1)
    {        
        
        double BEGFR=0;
        double DEGFR=0;
        double EGF=1000;
        double UEGFR=500;
        double inactiveC3G=1200;
        double activeC3G=0;
        double activeRap1=0;
        double BRaf=1500;
        double BRafPP=0;
        double activeRas=0;
        double Gap=2400;
        double inactiveRap1=1200;
        double inactiveRas=1200;
        double activeSOS=0;
        double inactiveSOS=1200;
        
        
        y[0]=BEGFR;
        y[1]=DEGFR;
        y[2]=EGF;
        y[3]=UEGFR;
        y[4]=inactiveC3G;
        y[5]=activeC3G;
        y[6]=activeRap1;
        y[7]=BRaf;
        y[8]=BRafPP;
        y[9]=activeRas;
        y[10]=Gap;
        y[11]=inactiveRap1;
        y[12]=inactiveRas;
        y[13]=activeSOS;
        y[14]=inactiveSOS;
         
        
    }
    if (modelProblem==4)
    {
        double A=0;
        double B=0;
        double C=10;
        double D=10;
        double E=0;
        double F=10;
        double G=0;
        double H=0;
        double I=0;
        double J=0;
        double K=0;
        
        y[0]=A;
        y[1]=B;
        y[2]=C;
        y[3]=D;
        y[4]=E;
        y[5]=F;
        y[6]=G;
        y[7]=H;
        y[8]=I;
        y[9]=J;
        y[10]=K;
             
    }
    if (modelProblem==5)
    {
       double A=0;
       double B=10;
       double C=10;
       double D=0;
       double E=0;
       double F=10;
       double G=10;
       double H=0;
       double I=0;
       double J=0;
       
        y[0]=A;
        y[1]=B;
        y[2]=C;
        y[3]=D;
        y[4]=E;
        y[5]=F;
        y[6]=G;
        y[7]=H;
        y[8]=I;
        y[9]=J;
               
    }
    if (modelProblem==6)
    {
       double A=0;
       double B=10;
       double C=10;
       double D=0;
       double E=0;
       double F=0;
       
       y[0]=A;
       y[1]=B;
       y[2]=C;
       y[3]=D;
       y[4]=E;
       y[5]=F;
        
    }
    
    
}

void SOFCAnodeIntegrator::getInitialConditionsIDA(doublereal t0,size_t leny,doublereal* y, doublereal* ydot)
{
}

void SOFCAnodeIntegrator::setID(doublereal* Id)
{
    
}

void SOFCAnodeIntegrator::printInitialConditions(double* a, double* b)
{
}


/*
 * int SOFCAnodeIntegrator::nparams()
 * {
 * return m_ntotpar;
 * }
 */
void SOFCAnodeIntegrator::initialize(double t0)
{
    m_size.push_back(neqvar);
    m_ydot.resize(neqvar,0.0);
    m_atol.resize(neq());
    fill(m_atol.begin(), m_atol.end(), m_atols);
    m_integ->setTolerances(m_rtol, neq(), DATA_PTR(m_atol));
    m_integ->setSensitivityTolerances(m_rtolsens, m_atolsens);
    m_integ->setMaxStepSize(m_maxstep);
    m_integ->initialize(t0, *this);
    m_init = true;
}

double* SOFCAnodeIntegrator::advance(int nData, double* time)
{
    
    double* modelOutput;
    
    modelOutput= new double [nData*(nReactions+neqvar)];
    
    if (!m_init)
    {
        if (m_maxstep < 0.0)
            m_maxstep = time[0] - m_time;
        initialize(0.0);
    }
    
    
    // Call to the implicit solver
    int l=1;
    double sum=1E10;
    int x;
    
    
    for (int t=0;t<nData;t++)
    {
                
        x=m_integ->integrate(time[t]);
        sol=m_integ->solution();
        ++l;
        
        // Get reactions rates
        reactionRates(sol);
        
        // Get species rates
        speciesRates(sol);
                
        // Reaction rates
         
        
        if (x==-1)
        { 
            
          // Note concentrations cannot be negative; this is used just as a reference  
          for (int i=0;i<nReactions;i++)
          {
              modelOutput[t*(nReactions+neqvar)+i]=x;
        
          }    
          // Species concentrations
        
          for (int i=nReactions;i<nReactions+neqvar;i++)
              modelOutput[t*(nReactions+neqvar)+i]=x;  

        }                
        else
        {
          for (int i=0;i<nReactions;i++)
              modelOutput[t*(nReactions+neqvar)+i]=reactionYDT[i];
        
          // Species concentrations
        
          for (int i=nReactions;i<nReactions+neqvar;i++)
              modelOutput[t*(nReactions+neqvar)+i]=sol[i-nReactions];//speciesYDT[i-nReactions];
        }
        
    }
    
    return modelOutput;
}

void SOFCAnodeIntegrator::reactionRates(double* sol)
{
    if (modelProblem==2)
    {
        double uBEGFR=sol[0];
        double remcRaf=sol[1];
        double remSOS=sol[2];
        double inactiveSOS=sol[3];
        double inactiveRas=sol[4];
        double inactiveRap1=sol[5];
        double inactivePKA=sol[6];
        double inactiveEPAC=sol[7];
        double cRafPP=sol[8];
        double cRaf=sol[9];
        double bEGFR=sol[10];
        double activeSOS=sol[11];
        double activeRas=sol[12];
        double activeRap1=sol[13];
        double PKA=sol[14];
        double MEKPP=sol[15];
        double MEK=sol[16];
        double ERKPP=sol[17];
        double ERK=sol[18];
        double EPAC=sol[19];
        double EGF=sol[20];
        double BRafPP=sol[21];
        double BRaf=sol[22];
        double activeC3G=sol[23];
        double inactiveC3G=sol[24];
        double dEGFR=sol[25];
        double EPACA=sol[26];
        double Gap=sol[27];
        double PKAA=sol[28];
        double Cilostamide=sol[29];
        
        reactionYDT[0]=k1*bEGFR*inactiveSOS/(k1p+inactiveSOS);
        reactionYDT[1]=k2*EGF*uBEGFR-k2r*bEGFR;
        reactionYDT[2]=k3*activeSOS/(k3p+activeSOS);
        reactionYDT[3]=k4*activeSOS*inactiveRas/(k4p+inactiveRas);
        reactionYDT[4]=k5*Gap*activeRas/(k5p+activeRas);
        reactionYDT[5]=k6*activeRas*cRaf/(k6p+cRaf);
        reactionYDT[6]=k7*cRafPP/(k7p+cRafPP);
        reactionYDT[7]=k8*cRafPP*MEK/(k8p+MEK);
        reactionYDT[8]=k9*MEKPP/(k9p+MEKPP);
        reactionYDT[9]=k10*MEKPP*ERK/(k10p+ERK);
        reactionYDT[10]=k11*ERKPP*inactiveSOS/(k11p+inactiveSOS);
        reactionYDT[11]=k12*ERKPP*activeSOS/(k12p+activeSOS);
        reactionYDT[12]=k13*PKA*cRaf/(k13p+cRaf);
        reactionYDT[13]=k14*PKAA*inactivePKA/(k14p+inactivePKA);
        reactionYDT[14]=k15*Cilostamide*inactivePKA/(k15p+inactivePKA);
        reactionYDT[15]=k16*PKA/(k16p+PKA);
        reactionYDT[16]=k17*EPACA*inactiveEPAC/(k17p+inactiveEPAC);
        reactionYDT[17]=k18*Cilostamide*inactiveEPAC/(k18p+inactiveEPAC);
        reactionYDT[18]=k19*EPAC/(k19p+EPAC);
        reactionYDT[19]=k20*EPAC*inactiveRap1/(k20p+inactiveRap1);
        reactionYDT[20]=k21*Gap*activeRap1/(k21p+activeRap1);
        reactionYDT[21]=k22*activeRap1*BRaf/(k22p+BRaf);
        reactionYDT[22]=k23*BRafPP/(k23p+BRafPP);
        reactionYDT[23]=k24*BRafPP*MEK/(k24p+MEK);
        reactionYDT[24]=k25*bEGFR*inactiveC3G/(k25p+inactiveC3G);
        reactionYDT[25]=k26*activeC3G;
        reactionYDT[26]=k27*activeC3G*inactiveRap1/(k27p+inactiveRap1);
        reactionYDT[27]=k28*activeRas*BRaf/(k28p+BRaf);
        reactionYDT[28]=k29*ERKPP/(k29p+ERKPP);
        reactionYDT[29]=k30*bEGFR;      
                       
    }
    else if (modelProblem==1)
    {
        double BEGFR=sol[0];
        double DEGFR=sol[1];
        double EGF=sol[2];
        double UEGFR=sol[3];
        double inactiveC3G=sol[4];
        double activeC3G=sol[5];
        double activeRap1=sol[6];
        double BRaf=sol[7];
        double BRafPP=sol[8];
        double activeRas=sol[9];
        double Gap=sol[10];
        double inactiveRap1=sol[11];
        double inactiveRas=sol[12];
        double activeSOS=sol[13];
        double inactiveSOS=sol[14];
        
        reactionYDT[0]=p56*DEGFR;
        reactionYDT[1]=p5*EGF*UEGFR-p6*BEGFR;
        reactionYDT[2]=p47*BEGFR*inactiveC3G/(p48+inactiveC3G);
        reactionYDT[3]=p49*activeC3G;
        reactionYDT[4]=p50*activeC3G*inactiveRap1/(p51+inactiveRap1);
        reactionYDT[5]=p41*activeRap1*BRaf/(p42+BRaf);
        reactionYDT[6]=p44*BRafPP/(p43+BRafPP);
        reactionYDT[7]=p52*activeRas*BRaf/(p53+BRaf);
        reactionYDT[8]=p39*Gap*activeRap1/(p40+activeRap1);
        reactionYDT[9]=p11*Gap*activeRas/(p12+activeRas);
        reactionYDT[10]=p9*activeSOS*inactiveRas/(p10+inactiveRas);
        reactionYDT[11]=p8*activeSOS/(p7+activeSOS);
        reactionYDT[12]=p3*BEGFR*inactiveSOS/(p4+inactiveSOS);
    }
    if (modelProblem==3)
    {
        double uBEGFR=sol[0];
        double remcRaf=sol[1];
        double remSOS=sol[2];
        double inactiveSOS=sol[3];
        double inactiveRas=sol[4];
        double inactivePKA=sol[5];
        double cRafPP=sol[6];
        double cRaf=sol[7];
        double bEGFR=sol[8];
        double activeSOS=sol[9];
        double activeRas=sol[10];
        double PKA=sol[11];
        double MEKPP=sol[12];
        double MEK=sol[13];
        double ERKPP=sol[14];
        double ERK=sol[15];
        double EGF=sol[16];
        double Gap=sol[17];
        double PKAA=sol[18];
        double Cilostamide=sol[19];
        
        reactionYDT[0]=k1*bEGFR*inactiveSOS/(k1p+inactiveSOS);
        reactionYDT[1]=k2*EGF*uBEGFR-k2r*bEGFR;
        reactionYDT[2]=k3*activeSOS/(k3p+activeSOS);
        reactionYDT[3]=k4*activeSOS*inactiveRas/(k4p+inactiveRas);
        reactionYDT[4]=k5*Gap*activeRas/(k5p+activeRas);
        reactionYDT[5]=k6*activeRas*cRaf/(k6p+cRaf);
        reactionYDT[6]=k7*cRafPP/(k7p+cRafPP);
        reactionYDT[7]=k8*cRafPP*MEK/(k8p+MEK);
        reactionYDT[8]=k9*MEKPP/(k9p+MEKPP);
        reactionYDT[9]=k10*MEKPP*ERK/(k10p+ERK);
        reactionYDT[10]=k11*ERKPP*inactiveSOS/(k11p+inactiveSOS);
        reactionYDT[11]=k12*ERKPP*activeSOS/(k12p+activeSOS);
        reactionYDT[12]=k13*PKA*cRaf/(k13p+cRaf);
        reactionYDT[13]=k14*PKAA*inactivePKA/(k14p+inactivePKA);
        reactionYDT[14]=k15*Cilostamide*inactivePKA/(k15p+inactivePKA);
        reactionYDT[15]=k16*PKA/(k16p+PKA);
        reactionYDT[16]=k17*ERKPP/(k17p+ERKPP);
     
    }
    if (modelProblem==4)
    {
       double A=sol[0];
       double B=sol[1];
       double C=sol[2];
       double D=sol[3];
       double E=sol[4];
       double F=sol[5];
       double G=sol[6];
       double H=sol[7];
       double I=sol[8];
       double J=sol[9];
       double K=sol[10];
       
       reactionYDT[0]=k1*A;
       reactionYDT[1]=k2*C*D;
       reactionYDT[2]=k3*A;
       reactionYDT[3]=k4*E*F;
       reactionYDT[4]=k5*A;
       reactionYDT[5]=k6*H;
       reactionYDT[6]=k7*H;
       reactionYDT[7]=k8*F*J;
       reactionYDT[8]=k9*I;
       reactionYDT[9]=k10*G;
    }
    if (modelProblem==5)
    {
       double A=sol[0];
       double B=sol[1];
       double C=sol[2];
       double D=sol[3];
       double E=sol[4];
       double F=sol[5];
       double G=sol[6];
       double H=sol[7];
       double I=sol[8];
       
       reactionYDT[0]=modelProblem5k1f*B*C-modelProblem5k1r*A;
       reactionYDT[1]=modelProblem5k2*A;
       reactionYDT[2]=modelProblem5k3*D*F;
       reactionYDT[3]=modelProblem5k4*F*G;
       reactionYDT[4]=modelProblem5k5*A;
       reactionYDT[5]=modelProblem5k6*E*G;    
    }
    if (modelProblem==6)
    {
       double A=sol[0];
       double B=sol[1];
       double C=sol[2];
       double D=sol[3];
       double E=sol[4];
       double F=sol[5];
       
       reactionYDT[0]=modelProblem6k1f*B*C-modelProblem6k1r*A;
       reactionYDT[1]=modelProblem6k2*A;
       reactionYDT[2]=modelProblem6k3*A;
       reactionYDT[3]=modelProblem6k4*D;
       reactionYDT[4]=modelProblem6k5*E;
       
    }
    
}

void SOFCAnodeIntegrator::speciesRates(double* sol)
{
    if (modelProblem==2)
    {
        
        double uBEGFR=sol[0];
        double remcRaf=sol[1];
        double remSOS=sol[2];
        double inactiveSOS=sol[3];
        double inactiveRas=sol[4];
        double inactiveRap1=sol[5];
        double inactivePKA=sol[6];
        double inactiveEPAC=sol[7];
        double cRafPP=sol[8];
        double cRaf=sol[9];
        double bEGFR=sol[10];
        double activeSOS=sol[11];
        double activeRas=sol[12];
        double activeRap1=sol[13];
        double PKA=sol[14];
        double MEKPP=sol[15];
        double MEK=sol[16];
        double ERKPP=sol[17];
        double ERK=sol[18];
        double EPAC=sol[19];
        double EGF=sol[20];
        double BRafPP=sol[21];
        double BRaf=sol[22];
        double activeC3G=sol[23];
        double inactiveC3G=sol[24];
        double dEGFR=sol[25];
        double EPACA=sol[26];
        double Gap=sol[27];
        double PKAA=sol[28];
        double Cilostamide=sol[29];        

        speciesYDT[0]=-k2*EGF*uBEGFR+k2r*bEGFR;
                
        speciesYDT[1]=k13*PKA*cRaf/(k13p+cRaf);
                
        speciesYDT[2]=k11*ERKPP*inactiveSOS/(k11p+inactiveSOS)+k12*ERKPP*activeSOS/(k12p+activeSOS);        
        
        speciesYDT[3]=-k1*bEGFR*inactiveSOS/(k1p+inactiveSOS)+k3*activeSOS/(k3p+activeSOS)-k11*ERKPP*inactiveSOS/(k11p+inactiveSOS);
                
        speciesYDT[4]=-k4*activeSOS*inactiveRas/(k4p+inactiveRas)+k5*Gap*activeRas/(k5p+activeRas);
                
        speciesYDT[5]=-k20*EPAC*inactiveRap1/(k20p+inactiveRap1)+k21*Gap*activeRap1/(k21p+activeRap1)-k27*activeC3G*inactiveRap1/(k27p+inactiveRap1);
                
        speciesYDT[6]=-k14*PKAA*inactivePKA/(k14p+inactivePKA)-k15*Cilostamide*inactivePKA/(k15p+inactivePKA)+k16*PKA/(k16p+PKA);
                
        speciesYDT[7]=-k17*EPACA*inactiveEPAC/(k17p+inactiveEPAC)-k18*Cilostamide*inactiveEPAC/(k18p+inactiveEPAC)+k19*EPAC/(k19p+EPAC);
                
        speciesYDT[8]=k6*activeRas*cRaf/(k6p+cRaf)-k7*cRafPP/(k7p+cRafPP);
                
        speciesYDT[9]=-k6*activeRas*cRaf/(k6p+cRaf)+k7*cRafPP/(k7p+cRafPP)-k13*PKA*cRaf/(k13p+cRaf);
                
        speciesYDT[10]=k2*EGF*uBEGFR-k2r*bEGFR-k30*bEGFR;     
                
        speciesYDT[11]=k1*bEGFR*inactiveSOS/(k1p+inactiveSOS)-k3*activeSOS/(k3p+activeSOS)-k12*ERKPP*activeSOS/(k12p+activeSOS);
                
        speciesYDT[12]=k4*activeSOS*inactiveRas/(k4p+inactiveRas)-k5*Gap*activeRas/(k5p+activeRas);        
        
        speciesYDT[13]=k20*EPAC*inactiveRap1/(k20p+inactiveRap1)-k21*Gap*activeRap1/(k21p+activeRap1)+k27*activeC3G*inactiveRap1/(k27p+inactiveRap1);
                
        speciesYDT[14]=k14*PKAA*inactivePKA/(k14p+inactivePKA)+k15*Cilostamide*inactivePKA/(k15p+inactivePKA)-k16*PKA/(k16p+PKA);
                
        speciesYDT[15]=k8*cRafPP*MEK/(k8p+MEK)-k9*MEKPP/(k9p+MEKPP)+k24*BRafPP*MEK/(k24p+MEK);
                
        speciesYDT[16]=-k8*cRafPP*MEK/(k8p+MEK)+k9*MEKPP/(k9p+MEKPP)-k24*BRafPP*MEK/(k24p+MEK);
                
        speciesYDT[17]=k10*MEKPP*ERK/(k10p+ERK)-k29*ERKPP/(k29p+ERKPP);
                
        speciesYDT[18]=-k10*MEKPP*ERK/(k10p+ERK)+k29*ERKPP/(k29p+ERKPP);
                
        speciesYDT[19]=k17*EPACA*inactiveEPAC/(k17p+inactiveEPAC)+k18*Cilostamide*inactiveEPAC/(k18p+inactiveEPAC)-k19*EPAC/(k19p+EPAC);
                
        speciesYDT[20]=-k2*EGF*uBEGFR+k2r*bEGFR;  
                
        speciesYDT[21]=k22*activeRap1*BRaf/(k22p+BRaf)-k23*BRafPP/(k23p+BRafPP)+k28*activeRas*BRaf/(k28p+BRaf);
                
        speciesYDT[22]=-k22*activeRap1*BRaf/(k22p+BRaf)+k23*BRafPP/(k23p+BRafPP)-k28*activeRas*BRaf/(k28p+BRaf);
                
        speciesYDT[23]=k25*bEGFR*inactiveC3G/(k25p+inactiveC3G)-k26*activeC3G;        
        
        speciesYDT[24]=-k25*bEGFR*inactiveC3G/(k25p+inactiveC3G)+k26*activeC3G;
                
        speciesYDT[25]=k30*bEGFR;
                
        speciesYDT[26]=0;
                
        speciesYDT[27]=0;
                        
        speciesYDT[28]=0;
                
        speciesYDT[29]=0;        
        

    }
    if (modelProblem==1)
    {
        double BEGFR=sol[0];
        double DEGFR=sol[1];
        double EGF=sol[2];
        double UEGFR=sol[3];
        double inactiveC3G=sol[4];
        double activeC3G=sol[5];
        double activeRap1=sol[6];
        double BRaf=sol[7];
        double BRafPP=sol[8];
        double activeRas=sol[9];
        double Gap=sol[10];
        double inactiveRap1=sol[11];
        double inactiveRas=sol[12];
        double activeSOS=sol[13];
        double inactiveSOS=sol[14];
        
        speciesYDT[0]=p5*EGF*UEGFR-p6*BEGFR-p56*BEGFR;
        
        speciesYDT[1]=p56*BEGFR;
        
        speciesYDT[2]=-p5*EGF*UEGFR+p6*BEGFR;
        
        speciesYDT[3]=-p5*EGF*UEGFR+p6*BEGFR;
        
        speciesYDT[4]=-p47*BEGFR*inactiveC3G/(p48+inactiveC3G)+p49*activeC3G;
        
        speciesYDT[5]=p47*BEGFR*inactiveC3G/(p48+inactiveC3G)-p49*activeC3G;
        
        speciesYDT[6]=-p39*Gap*activeRap1/(p40+activeRap1)+p50*activeC3G*inactiveRap1/(p51+inactiveRap1);
        
        speciesYDT[7]=-p41*activeRap1*BRaf/(p42+BRaf)+p44*BRafPP/(p43+BRafPP)-p52*activeRas*BRaf/(p53+BRaf);
        
        speciesYDT[8]=p41*activeRap1*BRaf/(p42+BRaf)-p44*BRafPP/(p43+BRafPP)+p52*activeRas*BRaf/(p53+BRaf);
        
        speciesYDT[9]=p9*activeSOS*inactiveRas/(p10+inactiveRas)-p11*Gap*activeRas/(p12+activeRas);
        
        speciesYDT[10]=0;
        
        speciesYDT[11]=-p50*activeC3G*inactiveRap1/(p51+inactiveRap1)+p39*Gap*activeRap1/(p40+activeRap1);
        
        speciesYDT[12]=-p9*activeSOS*inactiveRas/(p10+inactiveRas)+p11*Gap*activeRas/(p12+activeRas);
        
        speciesYDT[13]=p3*BEGFR*inactiveSOS/(p4+inactiveSOS)-p8*activeSOS/(p7+activeSOS);
        
        speciesYDT[14]=-p3*BEGFR*inactiveSOS/(p4+inactiveSOS)+p8*activeSOS/(p7+activeSOS);
    }
    if (modelProblem==3)
    {
        double uBEGFR=sol[0];
        double remcRaf=sol[1];
        double remSOS=sol[2];
        double inactiveSOS=sol[3];
        double inactiveRas=sol[4];
        double inactivePKA=sol[5];
        double cRafPP=sol[6];
        double cRaf=sol[7];
        double bEGFR=sol[8];
        double activeSOS=sol[9];
        double activeRas=sol[10];
        double PKA=sol[11];
        double MEKPP=sol[12];
        double MEK=sol[13];
        double ERKPP=sol[14];
        double ERK=sol[15];
        double EGF=sol[16];
        double Gap=sol[17];
        double PKAA=sol[18];
        double Cilostamide=sol[19];
        
        speciesYDT[0]=-k2*EGF*uBEGFR+k2r*bEGFR;
                
        speciesYDT[1]=k13*PKA*cRaf/(k13p+cRaf);
                
        speciesYDT[2]=k11*ERKPP*inactiveSOS/(k11p+inactiveSOS)+k12*ERKPP*activeSOS/(k12p+activeSOS);        
        
        speciesYDT[3]=-k1*bEGFR*inactiveSOS/(k1p+inactiveSOS)+k3*activeSOS/(k3p+activeSOS)-k11*ERKPP*inactiveSOS/(k11p+inactiveSOS);
                
        speciesYDT[4]=-k4*activeSOS*inactiveRas/(k4p+inactiveRas)+k5*Gap*activeRas/(k5p+activeRas);                
                
        speciesYDT[5]=-k14*PKAA*inactivePKA/(k14p+inactivePKA)-k15*Cilostamide*inactivePKA/(k15p+inactivePKA)+k16*PKA/(k16p+PKA);                
                
        speciesYDT[6]=k6*activeRas*cRaf/(k6p+cRaf)-k7*cRafPP/(k7p+cRafPP);
                
        speciesYDT[7]=-k6*activeRas*cRaf/(k6p+cRaf)+k7*cRafPP/(k7p+cRafPP)-k13*PKA*cRaf/(k13p+cRaf);
                
        speciesYDT[8]=k2*EGF*uBEGFR-k2r*bEGFR;     
                
        speciesYDT[9]=k1*bEGFR*inactiveSOS/(k1p+inactiveSOS)-k3*activeSOS/(k3p+activeSOS)-k12*ERKPP*activeSOS/(k12p+activeSOS);
                
        speciesYDT[10]=k4*activeSOS*inactiveRas/(k4p+inactiveRas)-k5*Gap*activeRas/(k5p+activeRas);        
                
        speciesYDT[11]=k14*PKAA*inactivePKA/(k14p+inactivePKA)+k15*Cilostamide*inactivePKA/(k15p+inactivePKA)-k16*PKA/(k16p+PKA);
                
        speciesYDT[12]=k8*cRafPP*MEK/(k8p+MEK)-k9*MEKPP/(k9p+MEKPP);
                
        speciesYDT[13]=-k8*cRafPP*MEK/(k8p+MEK)+k9*MEKPP/(k9p+MEKPP);
                
        speciesYDT[14]=k10*MEKPP*ERK/(k10p+ERK)-k17*ERKPP/(k17p+ERKPP);
                
        speciesYDT[15]=-k10*MEKPP*ERK/(k10p+ERK)+k17*ERKPP/(k17p+ERKPP);
                                
        speciesYDT[16]=-k2*EGF*uBEGFR+k2r*bEGFR;  
                                
        speciesYDT[17]=0;
                
        speciesYDT[18]=0;
                        
        speciesYDT[19]=0;                
    }
    if (modelProblem==4)
    {
       double A=sol[0];
       double B=sol[1];
       double C=sol[2];
       double D=sol[3];
       double E=sol[4];
       double F=sol[5];
       double G=sol[6];
       double H=sol[7];
       double I=sol[8];
       double J=sol[9];
       double K=sol[10];
       
        speciesYDT[0]=-k1*A-k3*A-k5*A+k2*C*D;
        speciesYDT[1]=k1*A;
        speciesYDT[2]=-k2*C*D;
        speciesYDT[3]=-k2*C*D;
        speciesYDT[4]=k3*A-k4*E*F;
        speciesYDT[5]=-k4*E*F-k8*F*J;
        speciesYDT[6]=k4*E*F-k10*G;
        speciesYDT[7]=k5*A-k6*H-k7*H;
        speciesYDT[8]=k6*H-k9*I;
        speciesYDT[9]=k7*H-k8*F*J;
        speciesYDT[10]=k8*F*J+k9*I+k10*G;       
       
    }
    if (modelProblem==5)
    {
       double A=sol[0];
       double B=sol[1];
       double C=sol[2];
       double D=sol[3];
       double E=sol[4];
       double F=sol[5];
       double G=sol[6];
       double H=sol[7];
       double I=sol[8];
       double J=sol[9];
       
       speciesYDT[0]=-modelProblem5k1r*A+modelProblem5k1f*B*C-modelProblem5k2*A-modelProblem5k5*A;
       
       speciesYDT[1]=modelProblem5k1r*A-modelProblem5k1f*B*C;
       
       speciesYDT[2]=modelProblem5k1r*A-modelProblem5k1f*B*C;
       
       speciesYDT[3]=modelProblem5k2*A-modelProblem5k3*D*F;
       
       speciesYDT[4]=modelProblem5k5*A-modelProblem5k6*E*G;
       
       speciesYDT[5]=-modelProblem5k3*D*F-modelProblem5k4*F*G;
       
       speciesYDT[6]=-modelProblem5k4*F*G-modelProblem5k6*E*G;
       
       speciesYDT[7]=modelProblem5k6*E*G;
       
       speciesYDT[8]=modelProblem5k4*F*G;
       
       speciesYDT[9]=modelProblem5k3*D*F;       
       
    }
    if (modelProblem==6)
    {
      double A=sol[0];
      double B=sol[1];
      double C=sol[2];
      double D=sol[3];
      double E=sol[4];
      double F=sol[5];
      
      speciesYDT[0]=modelProblem6k1f*B*C-modelProblem6k1r*A-modelProblem6k2*A-modelProblem6k3*A;
      
      speciesYDT[1]=-modelProblem6k1f*B*C+modelProblem6k1r*A;
      
      speciesYDT[2]=-modelProblem6k1f*B*C+modelProblem6k1r*A;
      
      speciesYDT[3]=modelProblem6k2*A-modelProblem6k4*D;
      
      speciesYDT[4]=modelProblem6k3*A-modelProblem6k5*E;
      
      speciesYDT[5]=modelProblem6k4*D+modelProblem6k5*E;        
        
    }
    
    
}




