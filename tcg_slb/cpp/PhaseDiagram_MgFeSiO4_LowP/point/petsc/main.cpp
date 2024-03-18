#include <math.h>
#include <numeric>
#include <algorithm>
#include "reactions.h"
#include "petsc.h"

// some function headers for PETSc
extern PetscErrorCode FormFunction(SNES snes, Vec x, Vec f,          // petsc snes callback function to form the residual
                                   void* ctx);

extern PetscErrorCode FormJacobian(SNES snes, Vec x, Mat A, Mat B,   // petsc snes callback function to form the jacobian
                                   void* ctx);

extern PetscErrorCode FormIFunction(TS ts, PetscReal tt, Vec u, Vec u_t, 
                                    Vec f,  void* ctx);              // petsc ts callback function to form the residual

extern PetscErrorCode FormRHSFunction(TS ts, PetscReal tt, Vec u,
                                    Vec f,  void* ctx);              // petsc ts callback function to form the residual

extern PetscErrorCode FormIJacobian(TS ts, PetscReal tt, Vec u, Vec u_t, 
                                    PetscReal a, Mat A, Mat B,       // petsc ts callback function to form the jacobian
                                    void* ctx);

extern PetscErrorCode FormRHSJacobian(TS ts, PetscReal tt, Vec u,
                                    Mat A, Mat B,                    // petsc ts callback function to form the jacobian
                                    void* ctx);

extern PetscErrorCode TSCustomPreStep(TS ts);

extern PetscErrorCode TSCustomPreStage(TS ts, PetscReal tt);

extern PetscErrorCode TSCustomPostStep(TS ts);

extern PetscErrorCode SNESCustomMonitor(SNES snes, PetscInt its, PetscReal norm, void* ctx);

extern PetscErrorCode SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu);

//*******************************************************************|************************************************************//
// a class returning the output from a reaction in a useful format
//*******************************************************************|************************************************************//
class PDReactiveODE
{

public:

  MgFeSiO4_all_slb_rx rxn;

  std::size_t I, K, J;
  std::vector<std::size_t> Kis;

  double T, p, Da, eps, rho0;
  std::vector<double> mi, mi0;
  std::vector<std::vector<double> > cik, csik, cik0;

  //*****************************************************************|************************************************************//
  // constructor
  //*****************************************************************|************************************************************//
  PDReactiveODE(double T0, double p0, std::size_t i0, double ci0k0, double eps0, double Da0)
  {
    cik = rxn.zero_C();
    I = cik.size();
    Kis.resize(I,1);
    mi.resize(I,0.);
    for (std::size_t i=0; i < I; ++i)
    {
      Kis[i] = cik[i].size();
    }
    K = std::accumulate(Kis.begin(), Kis.end(), 0.0);

    T = T0;
    p = p0;
    eps = eps0;
    Da = Da0;

    mi0.resize(I,0.);
    mi0[i0] = 1.0;

    cik0 = rxn.zero_C();
    for (std::size_t i=0; i<I; ++i)
    {
      cik0[i][0] = (Kis[i] == 1) ? 1.0 : ci0k0;
      for (std::size_t k=1; k<Kis[i]; ++k)
      {
        cik0[i][k] = (1. - ci0k0)/(Kis[i]-1.);
      }
    }

    mi.assign(mi0.begin(), mi0.end());
    for (std::size_t i=0; i < I; ++i) cik[i].assign(cik0[i].begin(), cik0[i].end());

    csik = rxn.zero_C();
    regularize_cik();
 
    rho0 = 1./v();
  }
  
  //*****************************************************************|************************************************************//
  // update mi
  //*****************************************************************|************************************************************//
  void update_mi(const PetscScalar *uu)
  {
    mi.assign(uu, uu + I);
  }

  //*****************************************************************|************************************************************//
  // update cik
  //*****************************************************************|************************************************************//
  void update_cik(const PetscScalar *uu)
  {
    std::size_t sKi = I;
    for (std::size_t i=0; i < I; ++i)
    {
      for (std::size_t k=0; k < Kis[i]; ++k)
      {
        cik[i][k] = uu[sKi + k];
      }
      sKi += Kis[i];
    }
    regularize_cik();
  }

  //*****************************************************************|************************************************************//
  // regularize cik into csik
  //*****************************************************************|************************************************************//
  void regularize_cik()
  {
    double ctot = 0.0;
    for (std::size_t i=0; i < I; ++i)
    {
      for (std::size_t k=0; k < Kis[i]; ++k)
      {
        csik[i][k] = std::min(std::max(cik[i][k], eps), 1.0-eps);
      }
      ctot = std::accumulate(csik[i].begin(), csik[i].end(), 0.0);
      for (std::size_t k=0; k < Kis[i]; ++k)
      {
        csik[i][k] /= ctot;
      }
    }
  }

  //*****************************************************************|************************************************************//
  // return phase densities
  //*****************************************************************|************************************************************//
  std::vector<double> rhoi()
  {
    std::vector<double> values(I,0.0);

    //Evaluate
    rxn.rho(T,p,csik,values);

    //std::cout << "rhoi: ";
    //for (std::size_t i=0; i < I; ++i) std::cout << values[i] << " ";
    //std::cout << "\n";

    return values;
  }

  //*****************************************************************|************************************************************//
  // return phase density derivatives
  //*****************************************************************|************************************************************//
  std::vector<std::vector<double> > drhoidu()
  {
    std::vector<std::vector<double> > values(I);

    //Evaluate
    std::size_t sKi = 0;
    for (std::size_t i=0; i < I; ++i)
    {
      values[i].resize(I+K, 0.0);
      std::vector<double> drhoidcik = rxn.phases()[i]->drho_dc(T,p,csik[i]);
      for (std::size_t dk=0; dk < Kis[i]; ++dk)
      {
        values[i][I + sKi + dk] = drhoidcik[dk];
      }
      sKi += Kis[i];
    }

    //std::cout << "drhoidu: ";
    //for (std::size_t i=0; i < I; ++i)
    //{
    //  for (std::size_t d=0; d < I + K; ++d) std::cout << values[i][d] << " ";
    //}
    //std::cout << "\n";

    return values;
  }

  //*****************************************************************|************************************************************//
  // return bulk specific volume
  //*****************************************************************|************************************************************//
  double v()
  {
    std::vector<double> work = rhoi();
    std::transform(mi.begin(), mi.end(), work.begin(), work.begin(), std::divides<double>());
    double value = std::accumulate(work.begin(), work.end(), 0.0);
    return value;
  }

  //*****************************************************************|************************************************************//
  // return bulk specific volume derivatives
  //*****************************************************************|************************************************************//
  std::vector<double> dvdu()
  {
    std::vector<double> work(I, 1.0);
    std::vector<double> _rhoi = rhoi();
    std::vector<std::vector<double> > _drhoidu = drhoidu();
    std::vector<double> values(I+K, 0.0);

    std::transform(work.begin(), work.end(), _rhoi.begin(), work.begin(), std::divides<double>());
    std::copy(work.begin(), work.end(), values.begin());

    std::transform(work.begin(), work.end(), work.begin(), work.begin(), std::multiplies<double>());
    std::transform(mi.cbegin(), mi.cend(), work.begin(), work.begin(), std::multiplies<double>());

    std::size_t sKi = I;
    for (std::size_t di=0; di < I; ++di)
    {
      for (std::size_t dk=0; dk < Kis[di]; ++dk)
      {
        for (std::size_t i=0; i < I; ++i)
        {
          values[sKi+dk] -= _drhoidu[i][sKi + dk]*work[i];
        }
      }
      sKi += Kis[di];
    }

    return values;
  }
  
  //*****************************************************************|************************************************************//
  // return phase sources
  //*****************************************************************|************************************************************//
  std::vector<double> Gammai()
  {
    std::vector<double> values(I,0.0);

    //Evaluate
    rxn.Gamma_i(T,p,csik,mi,values);

    //std::cout << "Gammai: ";
    //for (std::size_t i=0; i < I; ++i) std::cout << values[i] << " ";
    //std::cout << "\n";

    return values;
  }

  //*****************************************************************|************************************************************//
  // return phase source derivatives
  //*****************************************************************|************************************************************//
  std::vector<std::vector<double> > dGammaidu()
  {
    std::vector<std::vector<double> > values(I);
    std::vector<std::vector<double> > dGammaidmi(I);
    std::vector<std::vector<std::vector<double> > > dGammaidcik(I);
    
    //Evaluate
    rxn.dGamma_i_dPhi(T,p,csik,mi,dGammaidmi);
    rxn.dGamma_i_dC(T,p,csik,mi,dGammaidcik);
    std::size_t sKi = 0;
    for (std::size_t i=0; i < I; ++i)
    {
      values[i].resize(I+K,0.0);
      std::copy(dGammaidmi[i].begin(), dGammaidmi[i].end(), values[i].begin());
      sKi = I;
      for (int di=0; di < I; di++)
      {
        std::copy(dGammaidcik[i][di].begin(), dGammaidcik[i][di].end(), values[i].begin() + sKi);
        sKi += Kis[di];
      }
    }

    //std::cout << "dGammaidu: ";
    //for (std::size_t i=0; i < I; ++i)
    //{
    //  for (std::size_t d=0; d < I + K; ++d) std::cout << values[i][d] << " ";
    //}
    //std::cout << "\n";

    return values;
  }

  //*****************************************************************|************************************************************//
  // return component sources
  //*****************************************************************|************************************************************//
  std::vector<double> Gammaik()
  {
    std::vector<double> values(K, 0.0);
    std::vector<std::vector<double> > Gammaik(I);

    rxn.Gamma_ik(T,p,csik,mi,Gammaik);
    std::size_t sKi = 0;
    for (std::size_t i=0; i < I; ++i)
    {
      std::copy(Gammaik[i].begin(), Gammaik[i].end(), values.begin() + sKi);
      sKi += Kis[i];
    }

    //std::cout << "Gammaik: ";
    //for (std::size_t i=0; i < K; ++i) std::cout << values[i] << " ";
    //std::cout << "\n";

    return values;
  }

  //*****************************************************************|************************************************************//
  // return component source derivatives
  //*****************************************************************|************************************************************//
  std::vector<std::vector<double> > dGammaikdu()
  {
    std::vector<std::vector<double> > values(K);
    std::vector<std::vector<double> > dGammaikdmii;
    double dGammaikdcikik;

    std::size_t sKi = 0, sdKi = 0;
    for (std::size_t i=0; i<I; ++i)
    {
      rxn.dGamma_ik_dPhi(T,p,csik,mi,i,dGammaikdmii);
      for (std::size_t k=0; k<Kis[i]; ++k)
      {
        values[sKi + k].resize(I+K, 0.0);
        std::copy(dGammaikdmii[k].begin(), dGammaikdmii[k].end(), values[sKi + k].begin());
        sdKi = I;
        for (std::size_t di=0; di<I; ++di)
        {
          for (std::size_t dk=0; dk<Kis[di]; ++dk)
          {
            dGammaikdcikik = rxn.dGamma_ik_dC(T,p,csik,mi,i,k,di,dk);
            values[sKi + k][sdKi + dk] = dGammaikdcikik;
          }
          sdKi += Kis[di];
        }
      }
      sKi += Kis[i];
    }

    //std::cout << "dGammaikdu: ";
    //for (std::size_t i=0; i < K; ++i)
    //{
    //  for (std::size_t d=0; d < I + K; ++d) std::cout << values[i][d] << " ";
    //}
    //std::cout << "\n";

    return values;
  }

};

// a structure used to pass bucket data into SNES and TS callback functions
typedef struct {
  PDReactiveODE* ode;
  std::vector<double> t;
  std::vector<std::vector<double> > y;
  bool print_norms;
} Ctx;

// internal function headers
TS setup_ts(Ctx &ctx, double &idt, double &end);
void run_ts(TS &ts, Ctx &ctx);

//*******************************************************************|************************************************************//
// main
//*******************************************************************|************************************************************//
int main(int argc, char* argv[])
{
  Ctx ctx;

  PetscErrorCode perr;
  perr = PetscInitialize(&argc,&argv,(char*)0,(char*)0);

  PetscBool log = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL,"-l",&log,NULL);

  // log files
  if (log)
  {
    std::ostringstream debug_file;
    debug_file << "ts-log";
    if(std::freopen(debug_file.str().c_str(), "w", stdout) == NULL)
    {
      std::cout << "Failed to redirect stdio for debugging.";
    }
  }

  // do this immediately after the log file has been opened (if requested)
  perr = PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscReal idt = 1.e-6;
  perr = PetscOptionsGetReal(NULL, NULL,"-idt",&idt,NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscReal Da = 1.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-Da",&Da,NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscReal T = 1673.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-T",&T,NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscReal p = 150000.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-p",&p,NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscReal end = 100.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-end",&end,NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscReal eps = 1.e-2;
  perr = PetscOptionsGetReal(NULL, NULL,"-eps",&eps,NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscInt i0 = 0;
  perr = PetscOptionsGetInt(NULL, NULL,"-i0",&i0,NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscReal ci0k0 = 0.9;
  perr = PetscOptionsGetReal(NULL, NULL,"-ci0k0",&ci0k0,NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscBool pnorms = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL, "-print_norms", &pnorms, NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);
  ctx.print_norms = pnorms;

  ctx.ode = new PDReactiveODE(T, p, i0, ci0k0, eps, Da);

  std::size_t& I = (*(ctx.ode)).I;
  std::size_t& K = (*(ctx.ode)).K;
  std::vector<std::size_t>& Kis = (*(ctx.ode)).Kis;

  TS ts = setup_ts(ctx, idt, end);

  PetscInt maxsteps;
  perr = TSGetMaxSteps(ts, &maxsteps); CHKERRABORT(PETSC_COMM_WORLD,perr);

  ctx.t.reserve(maxsteps);
  ctx.y.resize(I+K);
  for (std::size_t i=0; i < I+K; ++i)
  {
    ctx.y[i].reserve(maxsteps);
  }

  run_ts(ts, ctx);

  perr = TSDestroy(&ts); CHKERRABORT(PETSC_COMM_WORLD, perr);

  return 0;
}

//*******************************************************************|************************************************************//
// setup the TS and its child SNES
//*******************************************************************|************************************************************//
TS setup_ts(Ctx &ctx, double& idt, double& end)
{
  PetscErrorCode perr, perr2;

  std::size_t& I = (*(ctx.ode)).I;
  std::size_t& K = (*(ctx.ode)).K;
  std::vector<std::size_t>& Kis = (*(ctx.ode)).Kis;

  TS ts;
  SNES snes;
  KSP ksp;
  PC pc;

  Mat JFmat, JGmat;
  Vec Fres, Gres;

  perr = TSCreate(PETSC_COMM_WORLD, &ts); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetApplicationContext(ts, &ctx); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetType(ts, TSBDF); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetTime(ts, 0.0); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetTimeStep(ts, idt); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetMaxTime(ts, end); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetMaxSteps(ts, 100000000); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetMaxSNESFailures(ts, -1); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRABORT(PETSC_COMM_WORLD,perr);
  
  perr = TSGetSNES(ts, &snes); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = SNESSetTolerances(snes, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 10, PETSC_DEFAULT); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscViewerAndFormat *svf;
  PetscViewer sviewer;
  perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ts),&sviewer); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = PetscViewerAndFormatCreate(sviewer,PETSC_VIEWER_DEFAULT,&svf); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = SNESMonitorSet(snes, (PetscErrorCode (*)(SNES,PetscInt,PetscReal,void*))SNESMonitorDefault, 
                                           svf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = SNESMonitorSet(snes, SNESCustomMonitor, &ctx, PETSC_NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = SNESSetType(snes, SNESVINEWTONRSLS); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = SNESGetKSP(snes, &ksp); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = KSPSetType(ksp, "preonly"); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = KSPGetPC(ksp, &pc); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = PCSetType(pc, "lu"); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = PCFactorSetMatSolverType(pc, "mumps"); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = TSSetFromOptions(ts); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = MatCreate(PETSC_COMM_WORLD, &JFmat); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = MatSetSizes(JFmat, PETSC_DECIDE, PETSC_DECIDE, I+K, I+K); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = MatSetType(JFmat, MATAIJ); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = MatSetUp(JFmat); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = MatCreate(PETSC_COMM_WORLD, &JGmat); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = MatSetSizes(JGmat, PETSC_DECIDE, PETSC_DECIDE, I+K, I+K); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = MatSetType(JGmat, MATAIJ); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = MatSetUp(JGmat); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = MatCreateVecs(JFmat, &Fres, PETSC_NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = MatCreateVecs(JGmat, &Gres, PETSC_NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);
  
  perr = TSSetIFunction(ts, Fres, FormIFunction, (void *) &ctx); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetRHSFunction(ts, Gres, FormRHSFunction, (void *) &ctx); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetIJacobian(ts, JFmat, JFmat, FormIJacobian, (void *) &ctx); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetRHSJacobian(ts, JGmat, JGmat, FormRHSJacobian, (void *) &ctx); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscViewerAndFormat *vf;
  PetscViewer viewer;
  perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ts),&viewer); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = PetscViewerAndFormatCreate(viewer,PETSC_VIEWER_DEFAULT,&vf); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSMonitorSet(ts, (PetscErrorCode (*)(TS,PetscInt,PetscReal,Vec,void*))TSMonitorDefault, 
                                           vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = TSSetPreStep(ts, TSCustomPreStep); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetPreStage(ts, TSCustomPreStage); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = TSSetPostStep(ts, TSCustomPostStep); CHKERRABORT(PETSC_COMM_WORLD,perr);
  //perr = TSSetSolution(ts, ui); CHKERRABORT(PETSC_COMM_WORLD,perr);

  SNESType snestype;
  perr = SNESGetType(snes,&snestype); CHKERRABORT(PETSC_COMM_WORLD,perr);

  if (strcmp(snestype, SNESVINEWTONRSLS) || strcmp(snestype, SNESVINEWTONSSLS))
  {
    std::cout << "applying bounds for snestype: " << snestype << std::endl;
    Vec ub;
    perr = VecDuplicate(Fres, &ub); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = VecSet(ub, 1.0); CHKERRABORT(PETSC_COMM_WORLD,perr);
    Vec lb;
    perr = VecDuplicate(Fres, &lb); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = VecSet(lb, 0.0); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = VecAssemblyBegin(ub); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = VecAssemblyBegin(lb); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = VecAssemblyEnd(ub); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = VecAssemblyEnd(lb); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = SNESVISetVariableBounds(snes, lb, ub); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = SNESVISetComputeVariableBounds(snes, SNESVIDummyComputeVariableBounds); CHKERRABORT(PETSC_COMM_WORLD,perr);
  }

  return ts;
}

//*****************************************************************|************************************************************//
// run the problem using TS
//*****************************************************************|************************************************************//
void run_ts(TS &ts, Ctx &ctx)
{
  PetscErrorCode perr, perr2;
  std::size_t& I = (*(ctx.ode)).I;
  std::size_t& K = (*(ctx.ode)).K;
  std::vector<std::size_t>& Kis = (*(ctx.ode)).Kis;

  Vec u, r;
  perr = TSGetIFunction(ts, &r, PETSC_NULL, PETSC_NULL); CHKERRABORT(PETSC_COMM_WORLD,perr);
  perr = VecDuplicate(r, &u); CHKERRABORT(PETSC_COMM_WORLD,perr);

  PetscScalar *uu;
  perr = VecGetArray(u, &uu); CHKERRABORT(PETSC_COMM_WORLD,perr);
 
  for (std::size_t i=0; i<I; ++i)
  {
    uu[i] = (*(ctx.ode)).mi0[i];
  }
  std::size_t sKi = I;
  for (std::size_t i=0; i<I; ++i)
  {
    for (std::size_t k=0; k<Kis[i]; ++k)
    {
      uu[sKi + k] = (*(ctx.ode)).cik0[i][k];
    }
    sKi += Kis[i];
  }

  perr = VecRestoreArray(u, &uu); CHKERRABORT(PETSC_COMM_WORLD,perr);

  perr = TSSolve(ts, u);
  TSConvergedReason reason;
  perr2 = TSGetConvergedReason(ts, &reason);
  std::cout << "  TSConvergedReason: " << reason << std::endl;
  CHKERRABORT(PETSC_COMM_WORLD, perr);
  CHKERRABORT(PETSC_COMM_WORLD, perr2);
}


//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the residual function
//*******************************************************************|************************************************************//
PetscErrorCode FormIFunction(TS ts, PetscReal tt, Vec u, Vec u_t, Vec f, void* ctx) // petsc ts callback function to form the residual
{
  Ctx *tsctx = (Ctx *)ctx;
  PetscErrorCode perr, perr2;

  std::cout << "  In FormIFunction" << std::endl;

  perr = VecCopy(u_t, f); CHKERRABORT(PETSC_COMM_WORLD, perr);

  if ((*tsctx).print_norms)
  {
    PetscReal norm;

    perr = VecNorm(u,NORM_2,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIFunction: 2-norm u = " << norm << std::endl;

    perr = VecNorm(u,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIFunction: inf-norm u = " << norm << std::endl;

    perr = VecNorm(u_t,NORM_2,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIFunction: 2-norm u_t = " << norm << std::endl;

    perr = VecNorm(u_t,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIFunction: inf-norm u_t = " << norm << std::endl;

    perr = VecNorm(f,NORM_2,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIFunction: 2-norm f = " << norm << std::endl;

    perr = VecNorm(f,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIFunction: inf-norm f = " << norm << std::endl;
  }
  
  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the residual function
//*******************************************************************|************************************************************//
PetscErrorCode FormRHSFunction(TS ts, PetscReal tt, Vec u, Vec f, void* ctx) // petsc ts callback function to form the residual
{
  PetscErrorCode perr, perr2;
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormRHSFunction" << std::endl;

  PetscScalar *ff;
  const PetscScalar *uu;

  perr = VecGetArray(f, &ff); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = VecGetArrayRead(u, &uu); CHKERRABORT(PETSC_COMM_WORLD, perr);

  (*(*tsctx).ode).update_mi(uu);
  (*(*tsctx).ode).update_cik(uu);

  double v = (*(*tsctx).ode).v();
  std::vector<double> Gammai = (*(*tsctx).ode).Gammai();
  std::vector<double> Gammaik = (*(*tsctx).ode).Gammaik();

  std::size_t& I = (*(*tsctx).ode).I;
  std::size_t& K = (*(*tsctx).ode).K;
  std::vector<std::size_t>& Kis = (*(*tsctx).ode).Kis;
  double& Da = (*(*tsctx).ode).Da;
  double& eps = (*(*tsctx).ode).eps;
  double& rho0 = (*(*tsctx).ode).rho0;
  std::vector<double>& mi = (*(*tsctx).ode).mi;
  std::vector<std::vector<double> >& cik = (*(*tsctx).ode).cik;

  double GikcGi;
  std::size_t sKi = 0;
  for (std::size_t i=0; i<I; ++i)
  {
    ff[i] = Da*rho0*Gammai[i]*v;
    for (std::size_t k=0; k<Kis[i]; ++k)
    {
      GikcGi = Gammaik[sKi + k] - cik[i][k]*Gammai[i];
      ff[I + sKi + k] = Da*rho0*GikcGi*v/(mi[i] + eps);
    }
    sKi += Kis[i];
  }

  perr = VecRestoreArray(f, &ff); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = VecRestoreArrayRead(u, &uu); CHKERRABORT(PETSC_COMM_WORLD, perr);

  if ((*tsctx).print_norms)
  {
    PetscReal norm;

    perr = VecNorm(u,NORM_2,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSFunction: 2-norm u = " << norm << std::endl;

    perr = VecNorm(u,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSFunction: inf-norm u = " << norm << std::endl;

    perr = VecNorm(f,NORM_2,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSFunction: 2-norm f = " << norm << std::endl;

    perr = VecNorm(f,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSFunction: inf-norm f = " << norm << std::endl;
  }
  
  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the jacobian function
//*******************************************************************|************************************************************//
PetscErrorCode FormIJacobian(TS ts, PetscReal tt, Vec u, Vec u_t, 
                                    PetscReal a, Mat A, Mat B,       // petsc ts callback function to form the jacobian
                                    void* ctx)
{
  PetscErrorCode perr, perr2;
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormIJacobian" << std::endl;
  std::cout << "    a (shift) = " << a << std::endl;

  Vec diag;
  perr = VecDuplicate(u_t, &diag); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = VecSet(diag, a); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = MatZeroEntries(A); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = MatDiagonalSet(A, diag, INSERT_VALUES); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, perr);

  if ((*tsctx).print_norms)
  {
    PetscReal norm;

    perr = VecNorm(u,NORM_2,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIJacobian: 2-norm u = " << norm << std::endl;

    perr = VecNorm(u,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIJacobian: inf-norm u = " << norm << std::endl;

    perr = VecNorm(u_t,NORM_2,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIJacobian: 2-norm u_t = " << norm << std::endl;

    perr = VecNorm(u_t,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIJacobian: inf-norm u_t = " << norm << std::endl;

    perr = MatNorm(A,NORM_FROBENIUS,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIJacobian: Frobenius norm A = " << norm << std::endl;

    perr = MatNorm(A,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIJacobian: inf-norm A = " << norm << std::endl;

    perr = MatNorm(B,NORM_FROBENIUS,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIJacobian: Frobenius norm B = " << norm << std::endl;

    perr = MatNorm(B,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormIJacobian: inf-norm B = " << norm << std::endl;
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the jacobian function
//*******************************************************************|************************************************************//
PetscErrorCode FormRHSJacobian(TS ts, PetscReal tt, Vec u, Mat A, Mat B,       // petsc ts callback function to form the jacobian
                                    void* ctx)
{
  PetscErrorCode perr, perr2;
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormRHSJacobian" << std::endl;

  const PetscScalar *uu;

  perr = VecGetArrayRead(u, &uu); CHKERRABORT(PETSC_COMM_WORLD, perr);

  (*(*tsctx).ode).update_mi(uu);
  (*(*tsctx).ode).update_cik(uu);

  double v = (*(*tsctx).ode).v();
  std::vector<double> dvdu = (*(*tsctx).ode).dvdu();
  std::vector<double> Gammai = (*(*tsctx).ode).Gammai();
  std::vector<double> Gammaik = (*(*tsctx).ode).Gammaik();
  std::vector<std::vector<double> > dGammaidu = (*(*tsctx).ode).dGammaidu();
  std::vector<std::vector<double> > dGammaikdu = (*(*tsctx).ode).dGammaikdu();

  std::size_t& I = (*(*tsctx).ode).I;
  std::size_t& K = (*(*tsctx).ode).K;
  std::vector<std::size_t>& Kis = (*(*tsctx).ode).Kis;
  double& Da = (*(*tsctx).ode).Da;
  double& eps = (*(*tsctx).ode).eps;
  double& rho0 = (*(*tsctx).ode).rho0;
  std::vector<double>& mi = (*(*tsctx).ode).mi;
  std::vector<std::vector<double> >& cik = (*(*tsctx).ode).cik;

  double GikcGi;
  std::vector<double> dGikcGi(I + K);
  PetscScalar vals[(I+K)*(I+K)];

  std::size_t sKi = 0;
  for (std::size_t i=0; i<I; ++i)
  {
    for (std::size_t d=0; d<I+K; ++d)
    {
      vals[i*(I + K) + d] = Da*rho0*(dGammaidu[i][d]*v + Gammai[i]*dvdu[d]);
    }
    for (std::size_t k=0; k<Kis[i]; ++k)
    {
      GikcGi = Gammaik[sKi + k] - cik[i][k]*Gammai[i];

      //dGammaikdus[sKi + k] - cik[i][k]*dGammaidus[i] - delta_ik*Gammai[i]

      std::transform(dGammaidu[i].begin(), dGammaidu[i].end(), dGikcGi.begin(), 
                     std::bind1st(std::multiplies<double>(), -cik[i][k]));
      std::transform(dGammaikdu[sKi + k].begin(), dGammaikdu[sKi + k].end(), dGikcGi.begin(), dGikcGi.begin(), std::plus<double>());
      dGikcGi[I + sKi + k] -= Gammai[i];
      for (std::size_t d=0; d<I+K; ++d)
      {
        vals[(I + sKi + k)*(I + K) + d] = Da*rho0*(dGikcGi[d]*v/(mi[i] + eps) + GikcGi*dvdu[d]/(mi[i] + eps));
      }
      vals[(I + sKi + k)*(I + K) + i] -= Da*rho0*GikcGi*v/std::pow(mi[i] + eps, 2);
    }
    sKi += Kis[i];
  }

  perr = VecRestoreArrayRead(u, &uu); CHKERRABORT(PETSC_COMM_WORLD, perr);

  PetscInt idx[I+K], jdx[I+K];
  for (std::size_t i=0; i<I+K; ++i) 
  {
    idx[i] = i;
    jdx[i] = i;
  }
  perr = MatSetValues(A, I+K, idx, I+K, jdx, vals, INSERT_VALUES); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, perr);
  perr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, perr);

  if ((*tsctx).print_norms)
  {
    PetscReal norm;

    perr = VecNorm(u,NORM_2,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSJacobian: 2-norm u = " << norm << std::endl;

    perr = VecNorm(u,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSJacobian: inf-norm u = " << norm << std::endl;

    perr = MatNorm(A,NORM_FROBENIUS,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSJacobian: Frobenius norm A = " << norm << std::endl;

    perr = MatNorm(A,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSJacobian: inf-norm A = " << norm << std::endl;

    perr = MatNorm(B,NORM_FROBENIUS,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSJacobian: Frobenius norm B = " << norm << std::endl;

    perr = MatNorm(B,NORM_INFINITY,&norm); CHKERRABORT(PETSC_COMM_WORLD, perr);
    std::cout << "    FormRHSJacobian: inf-norm B = " << norm << std::endl;
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts monitor callback function
//*******************************************************************|************************************************************//
PetscErrorCode TSCustomMonitor(TS ts, PetscInt i, PetscReal tt, Vec u, void* ctx)
{
  PetscErrorCode perr, perr2;
  Ctx *tsctx = (Ctx *)ctx;

  const PetscScalar *uu;

  std::size_t& I = (*(*tsctx).ode).I;
  std::size_t& K = (*(*tsctx).ode).K;
  std::vector<std::size_t>& Kis = (*(*tsctx).ode).Kis;

  perr = VecGetArrayRead(u, &uu); CHKERRABORT(PETSC_COMM_WORLD, perr);

  (*tsctx).t.push_back(tt);
  for (std::size_t i=0; i < I+K; ++i)
  {
    (*tsctx).y[i].push_back(uu[i]);
  }

  perr = VecRestoreArrayRead(u, &uu); CHKERRABORT(PETSC_COMM_WORLD, perr);

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that outputs diagnostics at the end of a timestep
//*******************************************************************|************************************************************//
PetscErrorCode TSCustomPreStep(TS ts)
{
  PetscFunctionReturn(0);
}

//*****************************************************************|************************************************************//
// define the petsc ts callback function that outputs diagnostics at the end of a timestep
//*****************************************************************|************************************************************//
PetscErrorCode TSCustomPreStage(TS ts, PetscReal tt)
{
  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that outputs diagnostics at the end of a timestep
//*******************************************************************|************************************************************//
PetscErrorCode TSCustomPostStep(TS ts)
{
  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc snes monitor callback function
//*******************************************************************|************************************************************//
PetscErrorCode SNESCustomMonitor(SNES snes, PetscInt n, PetscReal norm, void* ctx)
{
  PetscFunctionReturn(0);
}

PetscErrorCode SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu)
{                                                                    
  PetscFunctionReturn(0);
} 


