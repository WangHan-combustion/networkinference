#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>    // defines class IdealGasMix
#include <cantera/equilibrium.h>    // chemical equilibrium
#include <cantera/transport.h>      // transport properties
#include <cantera/numerics.h>
#include <cantera/integrators.h>
#include <cantera/kernel/Array.h>
#include <cantera/kernel/FuncEval.h>
#include "SOFCAnode.h"
#include <time.h>
#include <cantera/Edge.h>

using namespace Cantera;
using namespace std;
using namespace Cantera_CXX;

class SOFCAnodeIntegrator:public Cantera::FuncEval
{
 public:

 SOFCAnodeIntegrator(int,int,int,double*);

 ~SOFCAnodeIntegrator();

void getInitialConditionsIDA(doublereal t0,size_t leny,doublereal* y, doublereal* ydot);

void setID(doublereal* Id);

void printInitialConditions(double* a, double* b);

void reactionRates(double*);

void speciesRates(double*);

int neq();

void eval(doublereal,doublereal*,doublereal*,doublereal*);

void getInitialConditions(doublereal, size_t,doublereal*);

//int nparams();

void initialize(double);

double* advance(int,double*);

 /* Controls. */
 Controls primary;                     //!< Controls.

 int neqvar;
 int m_ntotpar;
 int nReactions;

 // Integrator 
 Integrator* m_integ;
 
 doublereal m_time;
 bool m_init;
 int m_nv; 
 vector_int m_size;
 vector_fp m_atol;
 doublereal m_rtol, m_rtolsens;
 doublereal m_atols, m_atolsens;
 doublereal m_maxstep;
 bool m_verbose;
 vector_int m_nparams;
 vector_int m_connect;
 vector_fp m_ydot;
 double dt;
   
 int modelProblem;
 
 double k1,k2,k2f,k2r,k3,k4,k5,k6,k7,k8,k9,k10;
 
 double k1p, k3p, k4p, k5p, k6p, k7p, k8p, k9p, k10p, k11, k11p, k12, k12p, k13, k13p, k14, k14p, k15, k15p, k16, k16p, k17, k17p, k18, k18p, k19, k19p, k20, k20p, k21, k21p, k22, k22p, k23, k23p, k24, k24p, k25, k25p, k26, k27, k27p, k28, k28p, k29, k29p, k30;
 
 
 double modelProblem5k1f, modelProblem5k1r, modelProblem5k2, modelProblem5k3, modelProblem5k4, modelProblem5k5, modelProblem5k6;
 
 double modelProblem6k1f, modelProblem6k1r, modelProblem6k2, modelProblem6k3, modelProblem6k4, modelProblem6k5;

 double* sol;
 vector_fp reactionYDT;
 vector_fp speciesYDT;

 int delta;
 
 double p1, p2, p2a, p2b, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38;
 double p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56;

};
