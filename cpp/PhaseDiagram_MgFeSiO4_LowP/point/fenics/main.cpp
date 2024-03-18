#include <dolfin.h>
#include <math.h>
#include "PhaseDiagramTS.h"
#include "reactions.h"

// Some shared pointer types
typedef std::shared_ptr< dolfin::Form > Form_ptr;
typedef std::shared_ptr< dolfin::Function > Function_ptr;
typedef std::shared_ptr< dolfin::GenericFunction > GenericFunction_ptr;
typedef std::shared_ptr< dolfin::Expression > Expression_ptr;
typedef std::shared_ptr< dolfin::Constant > Constant_ptr;
typedef std::shared_ptr< dolfin::Mesh > Mesh_ptr;
typedef std::shared_ptr< dolfin::FunctionSpace > FunctionSpace_ptr;
typedef std::shared_ptr< dolfin::XDMFFile > XDMFFile_ptr;

class TCGExpression : public dolfin::Expression
{

public:                                                            // available to everyone

  // setup up reactions
  MgFeSiO4_all_slb_rx rxn; 
  
  // number of phases
  int I, K;
  std::vector<int> Kis;
  
  TCGExpression() : initialized_(false)
  {
    init();                                                               // do nothing
    rename("TCGExpression", "TCG");
  }
    
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    std::cerr << "C++ expressions must be called using the eval(values, x, cell) interface.";
    throw std::runtime_error("std::runtime_error thrown");
  }

  void update_T(const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    dolfin::Array<double> T_values(1);
    T_ptr->eval(T_values, x, cell);
    T = T_values[0];
  }
  
  void update_P(const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    dolfin::Array<double> P_values(1);
    P_ptr->eval(P_values, x, cell);
    P = P_values[0];
  }

  void update_mi(const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    dolfin::Array<double> mi_values(1);
    //for (int i=0; i < I-1; i++)
    for (int i=0; i < I; i++)
    {
      mi_ptrs[i]->eval(mi_values, x, cell);
      mi[i] = mi_values[0];
    }
    //mi[I-1] = 1. - std::accumulate(mi.begin(), mi.end()-1, 0.);

    //std::cout << "mi = ";
    //for (std::size_t i = 0; i < I; ++i)
    //{
    //  std::cout << mi[i] << " ";
    //}
    //std::cout << std::endl;
  }
  
  void update_cik(const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    double ctot = 0.0;
    for (int i=0; i < I; i++)
    {
      dolfin::Array<double> cik_values(Kis[i]);
      cik_ptrs[i]->eval(cik_values, x, cell);
      for (int k=0; k < Kis[i]; k++)
      {
        cik[i][k] = std::min(std::max(cik_values[k], eps), 1.0-eps);
      }
      ctot = std::accumulate(cik[i].begin(), cik[i].end(), 0.0);
      for (int k=0; k < Kis[i]; k++)
      {
        cik[i][k] /= ctot;
      }

      //std::cout << "cik[" << i << "] = ";
      //for (std::size_t k = 0; k < Kis[i]; ++k)
      //{
      //  std::cout << cik[i][k] << " ";
      //}
      //std::cout << std::endl;
    }
  }

  void update(const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    update_T(x, cell);
    update_P(x, cell);
    update_cik(x, cell);
  }

  void zero(dolfin::Array<double>& values) const
  {
    for (int i=0; i<values.size(); i++)
    {
      values[i] = 0.0;
    }
  }

  //void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
  //{
  //  std::cerr << "eval method of TCGExpression virtual.";
  //  throw std::runtime_error("std::runtime_error thrown");
  //}
  
  void init()
  {
    if (!initialized_)
    {
      initialized_ = true;
      
      // list the phases and endmembers in these reactions
      //rxn.report();
      //rxn.list_parameters();
      
      // initialize concentration matrices
      cik = rxn.zero_C();
      
      // initialize vector lengths 
      I = cik.size();
      Kis.resize(I,1);
      mi.resize(I,0.);
      for (int i=0; i < I; i++)
      {
        Kis[i] = cik[i].size();
      }
      K = std::accumulate(Kis.begin(), Kis.end(), 0.0);
    }
  }

  void attach(GenericFunction_ptr T_in, GenericFunction_ptr P_in,
              std::vector<Function_ptr> mi_in,
              std::vector<Function_ptr> cik_in,
              double eps_in)
  {
    T_ptr = T_in;
    P_ptr = P_in;
    mi_ptrs = mi_in;
    cik_ptrs = cik_in;
    eps = eps_in;
  }

protected:

  bool initialized_;
  
  GenericFunction_ptr T_ptr, P_ptr;
  std::vector<Function_ptr> mi_ptrs, cik_ptrs;
  
  mutable double T, P;
  mutable std::vector<double> mi;
  mutable std::vector<std::vector<double> > cik;

  double eps;
};

class PhaseDensities : public TCGExpression
{

public:                                                            // available to everyone

  PhaseDensities() : TCGExpression()
  {
    _value_shape.resize(1);
    _value_shape[0] = I;
    _rhoi.resize(I);
    rename("PhaseDensities", "rhoi");
  }
  
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    update(x, cell);
    //update_mi(x, cell);
    zero(values);
    
    //Evaluate
    rxn.rho(T,P,cik,_rhoi);
    for (int i=0; i < I; i++)
    {
      values[i] = _rhoi[i];
    }

    //std::cout << "rhoi: ";
    //for (std::size_t i=0; i < I; ++i) std::cout << values[i] << " ";
    //std::cout << "\n";

    //std::cout << "x: " << std::endl;
    //for (int i=0; i<x.size(); i++)
    //{
    //  std::cout << "  " << x[i];
    //}
    //std::cout << std::endl;
    //std::cout << "rhoi: " << values.str(true) << std::endl;
  }
  
private:                                                           // only available to this class
  
  // vector for values
  mutable std::vector<double> _rhoi;
  
};

class PhaseDensityDerivatives : public TCGExpression
{

public:                                                            // available to everyone

  PhaseDensityDerivatives() : TCGExpression()
  {
    _value_shape.resize(2);
    _value_shape[0] = I;
    _value_shape[1] = I+K;
    rename("PhaseDensityDerivatives", "drhoidus");
  }
  
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    update(x, cell);
    //update_mi(x, cell);
    zero(values);
    
    //Evaluate
    int sKi = 0;
    for (int i=0; i < I; i++)
    {
      std::vector<double> drhoidcik = rxn.phases()[i]->drho_dc(T,P,cik[i]);
      for (int di=0; di < I; di++)
      {
        values[i*(I+K)+di] = 0;
      }
      for (int dk=0; dk < Kis[i]; dk++)
      {
        values[i*(I+K) + I + sKi + dk] = drhoidcik[dk];
      }
      sKi += Kis[i];
    }

    //std::cout << "drhoidu: ";
    //for (std::size_t i=0; i < I; ++i)
    //{
    //  for (std::size_t d=0; d < I + K; ++d) std::cout << values[i*(I+K) + d] << " ";
    //}
    //std::cout << "\n";

    //std::cout << "x: " << std::endl;
    //for (int i=0; i<x.size(); i++)
    //{
    //  std::cout << "  " << x[i];
    //}
    //std::cout << std::endl;
    //std::cout << "drhoidus: " << values.str(true) << std::endl;
  }
  
};

class PhaseSources : public TCGExpression
{

public:                                                            // available to everyone

  PhaseSources() : TCGExpression()
  {
    _value_shape.resize(1);
    _value_shape[0] = I;
    _Gammai.resize(I);
    rename("PhaseSources", "Gammai");
  }
  
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    update(x, cell);
    update_mi(x, cell);
    zero(values);
    
    //Evaluate
    rxn.Gamma_i(T,P,cik,mi,_Gammai);
    //for (int i=0; i < I-1; i++)
    for (int i=0; i < I; i++)
    {
      values[i] = _Gammai[i];
    }

    //std::cout << "Gammai: ";
    //for (std::size_t i=0; i < I; ++i) std::cout << values[i] << " ";
    //std::cout << "\n";

    //std::cout << "x: " << std::endl;
    //for (int i=0; i<x.size(); i++)
    //{
    //  std::cout << "  " << x[i];
    //}
    //std::cout << std::endl;
    //std::cout << "Gammai: " << values.str(true) << std::endl;
  }
  
private:                                                           // only available to this class
  
  // vector for values
  mutable std::vector<double> _Gammai;
  
};

class PhaseSourceDerivatives : public TCGExpression
{

public:                                                            // available to everyone

  PhaseSourceDerivatives() : TCGExpression()
  {
    _value_shape.resize(2);
    _value_shape[0] = I;
    _value_shape[1] = I+K;
    _dGammaidmi.resize(I);
    _dGammaidcik.resize(I);
    rename("PhaseSourceDerivatives", "dGammaidus");
  }
  
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    update(x, cell);
    update_mi(x, cell);
    zero(values);
    
    //Evaluate
    rxn.dGamma_i_dPhi(T,P,cik,mi,_dGammaidmi);
    rxn.dGamma_i_dC(T,P,cik,mi,_dGammaidcik);
    //for (int i=0; i < I-1; i++)
    int sKi = 0;
    for (int i=0; i < I; i++)
    {
      for (int di=0; di < I; di++)
      {
        values[sKi+di] = _dGammaidmi[i][di];
      }
      sKi += I;
      for (int di=0; di < I; di++)
      {
        for (int dk=0; dk < Kis[di]; dk++)
        {
          values[sKi+dk] = _dGammaidcik[i][di][dk];
        }
        sKi += Kis[di];
      }
    }

    //std::cout << "dGammaidu: ";
    //for (std::size_t i=0; i < I; ++i)
    //{
    //  for (std::size_t d=0; d < I + K; ++d) std::cout << values[i*(I+K)+d] << " ";
    //}
    //std::cout << "\n";

    //std::cout << "x: " << std::endl;
    //for (int i=0; i<x.size(); i++)
    //{
    //  std::cout << "  " << x[i];
    //}
    //std::cout << std::endl;
    //std::cout << "dGammaidus: " << values.str(true) << std::endl;
  }
  
private:                                                           // only available to this class
  
  // vector for values
  mutable std::vector<std::vector<double> > _dGammaidmi;
  mutable std::vector<std::vector<std::vector<double> > > _dGammaidcik;
  
};

class ComponentSources : public TCGExpression
{

public:                                                            // available to everyone

  ComponentSources() : TCGExpression()
  {
    _value_shape.resize(1);
    _value_shape[0] = K;
    _Gammaik.resize(I);
    rename("ComponentSources", "Gammaik");
  }
  
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    update(x, cell);
    update_mi(x, cell);
    zero(values);
    
    //Evaluate
    rxn.Gamma_ik(T,P,cik,mi,_Gammaik);
    //for (int i=0; i < I-1; i++)
    int sKi = 0;
    for (int i=0; i < I; i++)
    {
      for (int k=0; k < Kis[i]; k++)
      {
        values[sKi+k] = _Gammaik[i][k];
      }
      sKi += Kis[i];
    }

    //std::cout << "Gammaik: ";
    //for (std::size_t i=0; i < K; ++i) std::cout << values[i] << " ";
    //std::cout << "\n";

    //std::cout << "x: " << std::endl;
    //for (int i=0; i<x.size(); i++)
    //{
    //  std::cout << "  " << x[i];
    //}
    //std::cout << std::endl;
    ////std::cout << "cik: " << std::endl;
    ////for (int i=0; i < I; i++)
    ////{
    ////  for (int k=0; k < Kis[i]; k++)
    ////  {
    ////    std::cout << "  " << cik[i][k];
    ////  }
    ////  std::cout << std::endl;
    ////}
    //std::cout << "Gammaik: " << values.str(true) << std::endl;
  }

private:                                                           // only available to this class
  
  // vector for values
  mutable std::vector<std::vector<double> > _Gammaik;
  
};

class ComponentSourceDerivatives : public TCGExpression
{

public:                                                            // available to everyone

  ComponentSourceDerivatives() : TCGExpression()
  {
    _value_shape.resize(2);
    _value_shape[0] = K;
    _value_shape[1] = I+K;
    rename("ComponentSourceDerivatives", "dGammaikdus");
  }
  
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell &cell) const
  {
    update(x, cell);
    update_mi(x, cell);
    zero(values);
    
    //Evaluate
    //for (int i=0; i < I-1; i++)
    int sKi = 0;
    for (int i=0; i < I; i++)
    {
      rxn.dGamma_ik_dPhi(T,P,cik,mi,i,_dGammaikdmii);
      for (int k=0; k < Kis[i]; k++)
      {
        for (int di=0; di < I; di++)
        {
          values[sKi+di] = _dGammaikdmii[k][di];
        }
        sKi += I;
        for (int di=0; di < I; di++)
        {
          for (int dk=0; dk < Kis[di]; dk++)
          {
            _dGammaikdcikik = rxn.dGamma_ik_dC(T,P,cik,mi,i,k,di,dk);
            values[sKi+dk] = _dGammaikdcikik;
          }
          sKi += Kis[di];
        }
      }
    }

    //std::cout << "dGammaikdu: ";
    //for (std::size_t i=0; i < K; ++i)
    //{
    //  for (std::size_t d=0; d < I + K; ++d) std::cout << values[i*(I+K)+d] << " ";
    //}
    //std::cout << "\n";

    //std::cout << "x: " << std::endl;
    //for (int i=0; i<x.size(); i++)
    //{
    //  std::cout << "  " << x[i];
    //}
    //std::cout << std::endl;
    //std::cout << "dGammaikdus: " << values.str(true) << std::endl;
  }
  
private:                                                           // only available to this class
  
  // vector for values
  mutable std::vector<std::vector<double> > _dGammaikdmii;
  mutable double _dGammaikdcikik;
  
};

class InitialCondition : public TCGExpression
{
public:

  std::size_t i0;
  double ci0k0;

  InitialCondition(std::size_t i00, double ci0k00) : TCGExpression()
  {
    _value_shape.resize(1);
    _value_shape[0] = I+K;
    rename("InitialCondition", "u0");

    i0 = i00;
    ci0k0 = ci0k00;
  }

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    zero(values);
    values[i0] = 1.0;
    std::size_t sKi = I;
    for (std::size_t i=0; i<I; ++i)
    {
      values[sKi] = (Kis[i] == 1) ? 1.0 : ci0k0;
      for (std::size_t k=1; k<Kis[i]; ++k) values[sKi + k] = (1. - ci0k0)/(Kis[i]-1.);
      sKi += Kis[i];
    }
  }
};

// a structure used to pass bucket data into SNES and TS callback functions
typedef struct {
  Form_ptr JF, JG, F, G;
  std::vector< std::shared_ptr<const dolfin::DirichletBC> > bcs;
  Mesh_ptr mesh;
  FunctionSpace_ptr V;
  Function_ptr us_i, us_dot;
  std::vector<Function_ptr> mi_ptrs, cik_ptrs;
  Constant_ptr a;
  double *t;
  XDMFFile_ptr file_u, file_snes;
  double terrl2, terrlinf;
  bool snes_vis_monitor;
  bool print_norms;
  bool ident_zeros;
  std::size_t nx;
  std::map< Function_ptr, Expression_ptr > function_map, dfunction_map;
} Ctx;

// some function headers for PETSc
extern PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec u, Vec u_t, 
                                    Vec f,  void* ctx);              // petsc ts callback function to form the residual

extern PetscErrorCode FormRHSFunction(TS ts, PetscReal t, Vec u,
                                    Vec f,  void* ctx);              // petsc ts callback function to form the residual

extern PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec u, Vec u_t, 
                                    PetscReal a, Mat A, Mat B,       // petsc ts callback function to form the jacobian
                                    void* ctx);

extern PetscErrorCode FormRHSJacobian(TS ts, PetscReal t, Vec u,
                                    Mat A, Mat B,                    // petsc ts callback function to form the jacobian
                                    void* ctx);

extern PetscErrorCode TSCustomPreStep(TS ts);

extern PetscErrorCode TSCustomPreStage(TS ts, PetscReal t);

extern PetscErrorCode TSCustomPostStep(TS ts);

extern PetscErrorCode SNESCustomMonitor(SNES snes, PetscInt its, PetscReal norm, void* ctx);

extern PetscErrorCode SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu);

// internal function headers
TS setup_ts(Ctx &ctx);
void run_ts(TS &ts, Ctx &ctx);
void update_functions(Ctx &ctx);
void update_dfunctions(Ctx &ctx);

//*******************************************************************|************************************************************//
// main
//*******************************************************************|************************************************************//
int main(int argc, char* argv[])
{
  Ctx ctx;
  ctx.ident_zeros = false;

  PetscErrorCode perr, perr2;
  perr = PetscInitialize(&argc,&argv,(char*)0,(char*)0);

  // Create mesh
  PetscInt pnx = 1;
  perr = PetscOptionsGetInt(NULL,NULL,"-nx",&pnx,NULL);
  ctx.nx = pnx;
  ctx.mesh.reset( new dolfin::UnitIntervalMesh(ctx.nx) );

  PetscBool plog = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL,"-l",&plog,NULL);

  // log files
  if (plog)
  {
    std::ostringstream debug_file;
    debug_file << "ts-log." << ctx.nx;
    if(std::freopen(debug_file.str().c_str(), "w", stdout) == NULL)
    {
      std::cout << "Failed to redirect stdio for debugging.";
    }
  }
  std::stringstream buffer;
  buffer.str(""); buffer << "ts-phasediagram_" << ctx.nx << ".xdmf";

  // do this immediately after the log file has been opened (if requested)
  perr = PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  PetscReal pDa = 1.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-Da",&pDa,NULL);
  auto Da = std::make_shared<dolfin::Constant>(pDa);

  PetscReal pT = 1673.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-T",&pT,NULL);
  auto T_ptr = std::make_shared<dolfin::Constant>(pT);

  PetscReal pp = 150000.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-p",&pp,NULL);
  auto P_ptr = std::make_shared<dolfin::Constant>(pp);

  PetscReal peps = 1.e-2;
  perr = PetscOptionsGetReal(NULL, NULL,"-eps",&peps,NULL);
  auto eps = std::make_shared<dolfin::Constant>(peps);

  PetscInt i0 = 0;
  perr = PetscOptionsGetInt(NULL, NULL,"-i0",&i0,NULL);
  PetscReal ci0k0 = 0.9;
  perr = PetscOptionsGetReal(NULL, NULL,"-ci0k0",&ci0k0,NULL);
  auto u0 = std::make_shared<InitialCondition>(i0, ci0k0);

  PetscBool ppnorms = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL, "-print_norms", &ppnorms, NULL);
  ctx.print_norms = ppnorms;

  PetscBool pmonitor = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL, "-snes_vis_monitor", &pmonitor, NULL);
  ctx.snes_vis_monitor = pmonitor;

  PetscBool pfunctions = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL, "-use_functions", &pfunctions, NULL);

  PetscBool pvisualize = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL, "-ts_vis_monitor", &pvisualize, NULL);

  // Create function space
  ctx.V.reset( new PhaseDiagramTS::FunctionSpace(ctx.mesh) );

  ctx.us_i.reset( new dolfin::Function(ctx.V) );
  (*ctx.us_i).rename("Solution", "us_i");
  (*ctx.us_i).interpolate(*u0);

  TCGExpression tcg;

  ctx.mi_ptrs.resize(tcg.I);
  ctx.cik_ptrs.resize(tcg.I);
  for (std::size_t i = 0; i < tcg.I; ++i) 
  {
    ctx.mi_ptrs[i].reset( &(*ctx.us_i)[i], dolfin::NoDeleter() );
    ctx.mi_ptrs[i]->rename(tcg.rxn.phases()[i]->name()+"MassFraction", "m"+tcg.rxn.phases()[i]->abbrev());
    ctx.cik_ptrs[i].reset( &(*ctx.us_i)[tcg.I+i], dolfin::NoDeleter() );
    ctx.cik_ptrs[i]->rename(tcg.rxn.phases()[i]->name()+"ComponentMassFractions", "c"+tcg.rxn.phases()[i]->abbrev());
  }

  auto rhoi_i = std::make_shared<PhaseDensities>();
  rhoi_i->attach(T_ptr, P_ptr, ctx.mi_ptrs, ctx.cik_ptrs, double(*eps));
  auto drhoidus_i = std::make_shared<PhaseDensityDerivatives>();
  drhoidus_i->attach(T_ptr, P_ptr, ctx.mi_ptrs, ctx.cik_ptrs, double(*eps));
  auto Gammai_i = std::make_shared<PhaseSources>();
  Gammai_i->attach(T_ptr, P_ptr, ctx.mi_ptrs, ctx.cik_ptrs, double(*eps));
  auto dGammaidus_i = std::make_shared<PhaseSourceDerivatives>();
  dGammaidus_i->attach(T_ptr, P_ptr, ctx.mi_ptrs, ctx.cik_ptrs, double(*eps));
  auto Gammaik_i = std::make_shared<ComponentSources>();
  Gammaik_i->attach(T_ptr, P_ptr, ctx.mi_ptrs, ctx.cik_ptrs, double(*eps));
  auto dGammaikdus_i = std::make_shared<ComponentSourceDerivatives>();
  dGammaikdus_i->attach(T_ptr, P_ptr, ctx.mi_ptrs, ctx.cik_ptrs, double(*eps));

  std::vector<std::vector<double> > cik0(rhoi_i->I);
  for (std::size_t i=0; i<rhoi_i->I; ++i)
  {
    cik0[i].resize(rhoi_i->Kis[i]);
    cik0[i][0] = (rhoi_i->Kis[i] == 1) ? 1.0 : ci0k0;
    for (std::size_t k=1; k<rhoi_i->Kis[i]; ++k) cik0[i][k] = (1. - ci0k0)/(rhoi_i->Kis[i]-1.);
  }
  std::vector<double> rhoi0(rhoi_i->I, 0.0);
  rhoi_i->rxn.rho(pT,pp,cik0,rhoi0);
  auto rho0 = std::make_shared<dolfin::Constant>(rhoi0[i0]);
  std::cout << "rho0 = " << double(*rho0) << std::endl;
  std::cout << "Da = " << double(*Da) << std::endl;

  ctx.t = new double(0.0);

  // output files
  if (pvisualize)
  {
    ctx.file_u.reset( new dolfin::XDMFFile((*ctx.mesh).mpi_comm(), buffer.str()) );
    (*ctx.file_u).parameters["flush_output"] = true;
    (*ctx.file_u).parameters["functions_share_mesh"] = true;
    (*ctx.file_u).parameters["rewrite_function_mesh"] = false;
  }

  ctx.a.reset( new dolfin::Constant(666.e6) );

  ctx.us_dot.reset( new dolfin::Function(ctx.V) );
  (*ctx.us_dot).rename("TimeDerivative", "us_dot");
  (*(*ctx.us_dot).vector()).zero();

  ctx.JF.reset( new PhaseDiagramTS::Form_JF(ctx.V, ctx.V) );
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JF>(ctx.JF)).a = ctx.a;

  ctx.F.reset( new PhaseDiagramTS::Form_F(ctx.V) );
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_F>(ctx.F)).us_dot = ctx.us_dot;

  ctx.JG.reset( new PhaseDiagramTS::Form_JG(ctx.V, ctx.V) );
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).us_i = ctx.us_i;
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).Da = Da;
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).eps = eps;
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).rho0 = rho0;

  ctx.G.reset( new PhaseDiagramTS::Form_G(ctx.V) );
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).us_i = ctx.us_i;
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).Da = Da;
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).eps = eps;
  (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).rho0 = rho0;

  if (pfunctions)
  {
    auto V_rhoi_i = std::make_shared<PhaseDiagramTS::CoefficientSpace_rhoi_i>(ctx.mesh);
    auto f_rhoi_i = std::make_shared<dolfin::Function>(V_rhoi_i);
    (*f_rhoi_i).rename("PhaseDensitiesFunction", "rhoi_i");
    ctx.function_map[f_rhoi_i] = rhoi_i;
    auto V_drhoidus_i = std::make_shared<PhaseDiagramTS::CoefficientSpace_drhoidus_i>(ctx.mesh);
    auto f_drhoidus_i = std::make_shared<dolfin::Function>(V_drhoidus_i);
    (*f_drhoidus_i).rename("PhaseDensityDerivativesFunction", "drhoidus_i");
    ctx.dfunction_map[f_drhoidus_i] = drhoidus_i;
    auto V_Gammai_i = std::make_shared<PhaseDiagramTS::CoefficientSpace_Gammai_i>(ctx.mesh);
    auto f_Gammai_i = std::make_shared<dolfin::Function>(V_Gammai_i);
    (*f_Gammai_i).rename("PhaseSourcesFunction", "Gammai_i");
    ctx.function_map[f_Gammai_i] = Gammai_i;
    auto V_dGammaidus_i = std::make_shared<PhaseDiagramTS::CoefficientSpace_dGammaidus_i>(ctx.mesh);
    auto f_dGammaidus_i = std::make_shared<dolfin::Function>(V_dGammaidus_i);
    (*f_dGammaidus_i).rename("PhaseSourceDerivativesFunction", "dGammaidus_i");
    ctx.dfunction_map[f_dGammaidus_i] = dGammaidus_i;
    auto V_Gammaik_i = std::make_shared<PhaseDiagramTS::CoefficientSpace_Gammaik_i>(ctx.mesh);
    auto f_Gammaik_i = std::make_shared<dolfin::Function>(V_Gammaik_i);
    (*f_Gammaik_i).rename("ComponentSourcesFunction", "Gammaik_i");
    ctx.function_map[f_Gammaik_i] = Gammaik_i;
    auto V_dGammaikdus_i = std::make_shared<PhaseDiagramTS::CoefficientSpace_dGammaikdus_i>(ctx.mesh);
    auto f_dGammaikdus_i = std::make_shared<dolfin::Function>(V_dGammaikdus_i);
    (*f_dGammaikdus_i).rename("ComponentSourceDerivativesFunction", "dGammaikdus_i");
    ctx.dfunction_map[f_dGammaikdus_i] = dGammaikdus_i;

    update_functions(ctx);
    update_dfunctions(ctx);

    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).rhoi_i = f_rhoi_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).Gammai_i = f_Gammai_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).Gammaik_i = f_Gammaik_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).drhoidus_i = f_drhoidus_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).dGammaidus_i = f_dGammaidus_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).dGammaikdus_i = f_dGammaikdus_i;

    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).rhoi_i = f_rhoi_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).Gammai_i = f_Gammai_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).Gammaik_i = f_Gammaik_i;
  }
  else
  {
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).rhoi_i = rhoi_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).Gammai_i = Gammai_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).Gammaik_i = Gammaik_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).drhoidus_i = drhoidus_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).dGammaidus_i = dGammaidus_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_JG>(ctx.JG)).dGammaikdus_i = dGammaikdus_i;

    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).rhoi_i = rhoi_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).Gammai_i = Gammai_i;
    (*std::dynamic_pointer_cast<PhaseDiagramTS::Form_G>(ctx.G)).Gammaik_i = Gammaik_i;
  }

  TS ts = setup_ts(ctx);

  run_ts(ts, ctx);

  perr = TSDestroy(&ts); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  return 0;
}

//*******************************************************************|************************************************************//
// setup the TS and its child SNES
//*******************************************************************|************************************************************//
TS setup_ts(Ctx &ctx)
{
  TS ts;
  SNES snes;
  KSP ksp;
  PC pc;
  PetscErrorCode perr;

  perr = TSCreate((*ctx.mesh).mpi_comm(), &ts); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetApplicationContext(ts, &ctx); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetType(ts, TSBDF); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetTime(ts, 0.0); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetTimeStep(ts, 1.e-6); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetMaxTime(ts, 1.0); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetMaxSteps(ts, 100000000); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetMaxSNESFailures(ts, -1); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetTolerances(ts, 1.e-9, PETSC_NULL, 1.e-5, PETSC_NULL); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  
  perr = TSGetSNES(ts, &snes); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = SNESSetTolerances(snes, 1.e-6, 1.e-6, PETSC_DEFAULT, 10, PETSC_DEFAULT); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  PetscViewerAndFormat *svf;
  PetscViewer sviewer;
  perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ts),&sviewer); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = PetscViewerAndFormatCreate(sviewer,PETSC_VIEWER_DEFAULT,&svf); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = SNESMonitorSet(snes, (PetscErrorCode (*)(SNES,PetscInt,PetscReal,void*))SNESMonitorDefault, 
                                           svf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  if (ctx.snes_vis_monitor)
  {
    perr = SNESMonitorSet(snes, SNESCustomMonitor,                // set a custom snes monitor
                                          &ctx, PETSC_NULL); 
  }

  perr = SNESSetType(snes, SNESVINEWTONRSLS); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  perr = SNESGetKSP(snes, &ksp); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = KSPSetType(ksp, "preonly"); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  perr = KSPGetPC(ksp, &pc); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = PCSetType(pc, "lu"); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = PCFactorSetMatSolverType(pc, "mumps"); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  perr = TSSetFromOptions(ts); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  dolfin::PETScVector Fres;
  dolfin::Assembler Fassembler;
  Fassembler.assemble(Fres, (*ctx.F));

  perr = TSSetIFunction(ts, Fres.vec(), FormIFunction, (void *) &ctx); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  dolfin::PETScVector Gres;
  dolfin::Assembler Gassembler;
  Gassembler.assemble(Gres, (*ctx.G));

  perr = TSSetRHSFunction(ts, Gres.vec(), FormRHSFunction, (void *) &ctx); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  dolfin::PETScMatrix JFmatrix;
  dolfin::SystemAssembler JFsysassembler(ctx.JF, ctx.F, ctx.bcs);
  JFsysassembler.keep_diagonal = true;
  JFsysassembler.assemble(JFmatrix);
  if (ctx.ident_zeros)
  {
    JFmatrix.ident_zeros();
  }

  perr = TSSetIJacobian(ts, JFmatrix.mat(), JFmatrix.mat(), FormIJacobian, (void *) &ctx); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  dolfin::PETScMatrix JGmatrix;
  dolfin::SystemAssembler JGsysassembler(ctx.JG, ctx.G, ctx.bcs);
  JGsysassembler.keep_diagonal = true;
  JGsysassembler.assemble(JGmatrix);
  if (ctx.ident_zeros)
  {
    JGmatrix.ident_zeros();
  }

  perr = TSSetRHSJacobian(ts, JGmatrix.mat(), JGmatrix.mat(), FormRHSJacobian, (void *) &ctx); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  PetscViewerAndFormat *vf;
  PetscViewer viewer;
  perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ts),&viewer); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = PetscViewerAndFormatCreate(viewer,PETSC_VIEWER_DEFAULT,&vf); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSMonitorSet(ts, (PetscErrorCode (*)(TS,PetscInt,PetscReal,Vec,void*))TSMonitorDefault, 
                                           vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  perr = TSSetPreStep(ts, TSCustomPreStep); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetPreStage(ts, TSCustomPreStage); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  perr = TSSetPostStep(ts, TSCustomPostStep); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  //perr = TSSetSolution(ts, (*std::dynamic_pointer_cast<dolfin::PETScVector>((*ctx.us_i).vector())).vec()); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  SNESType snestype;
  perr = SNESGetType(snes,&snestype); CHKERRABORT((*ctx.mesh).mpi_comm(), perr);

  if (strcmp(snestype, SNESVINEWTONRSLS) || strcmp(snestype, SNESVINEWTONSSLS))
  {
    std::cout << "applying bounds for snestype: " << snestype << std::endl;
    dolfin::PETScVector ub(*std::dynamic_pointer_cast<dolfin::PETScVector>((*ctx.us_i).vector())); 
    std::vector<double> uv(ub.local_size(), 1.0);
    ub.set_local(uv);
    ub.apply("insert");
    dolfin::PETScVector lb(*std::dynamic_pointer_cast<dolfin::PETScVector>((*ctx.us_i).vector())); 
    std::vector<double> lv(lb.local_size(), 0.0);
    lb.set_local(lv);
    lb.apply("insert");

    perr = SNESVISetVariableBounds(snes, lb.vec(), ub.vec());
    CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
    //perr = SNESVISetComputeVariableBounds(snes, SNESVIDummyComputeVariableBounds);
    //CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  }

  return ts;
}

//*******************************************************************|************************************************************//
// run the problem using TS
//*******************************************************************|************************************************************//
void run_ts(TS &ts, Ctx &ctx)
{
  PetscErrorCode perr, perr2;

  if (ctx.file_u)
  {
    for (std::size_t i = 0; i < ctx.mi_ptrs.size(); ++i) 
    {
      (*ctx.file_u).write(*ctx.mi_ptrs[i], *ctx.t);
    }
    for (std::size_t i = 0; i < ctx.cik_ptrs.size(); ++i) 
    {
      (*ctx.file_u).write(*ctx.cik_ptrs[i], *ctx.t);
    }
  }

  dolfin::PETScVector work(*std::dynamic_pointer_cast<dolfin::PETScVector>((*ctx.us_i).vector())); 
  perr = TSSolve(ts, work.vec()); 
  TSConvergedReason reason;
  perr2 = TSGetConvergedReason(ts, &reason); 
  std::cout << "  TSConvergedReason: " << reason << std::endl;
  CHKERRABORT((*ctx.mesh).mpi_comm(), perr);
  CHKERRABORT((*ctx.mesh).mpi_comm(), perr2);
}

//*******************************************************************|************************************************************//
// update any functions we're using
//*******************************************************************|************************************************************//
void update_functions(Ctx &ctx)
{
  for (auto f_it = ctx.function_map.begin(); f_it != ctx.function_map.end(); ++f_it)
  {
    (*f_it->first).interpolate(*f_it->second);
    //std::cout << (*f_it->first).label() << " = ";
    //for (std::size_t i = 0; i < (*f_it->first).vector()->size(); ++i)
    //{
    //  std::cout << (*(*f_it->first).vector())[i] << " ";
    //}
    //std::cout << std::endl;
  }
}

//*******************************************************************|************************************************************//
// update any derivative functions we're using
//*******************************************************************|************************************************************//
void update_dfunctions(Ctx &ctx)
{
  for (auto f_it = ctx.dfunction_map.begin(); f_it != ctx.dfunction_map.end(); ++f_it)
  {
    (*f_it->first).interpolate(*f_it->second);
    //std::cout << (*f_it->first).label() << " = ";
    //for (std::size_t i = 0; i < (*f_it->first).vector()->size(); ++i)
    //{
    //  std::cout << (*(*f_it->first).vector())[i] << " ";
    //}
    //std::cout << std::endl;
  }
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the residual function
//*******************************************************************|************************************************************//
PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec u, Vec u_t, Vec f, void* ctx) // petsc ts callback function to form the residual
{
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormIFunction" << std::endl;

  // update the time so that the bcs are applied at the right time
  // i.e. this function can be called on any time-level
  (*(*tsctx).t) = t;

  dolfin::PETScVector iteratedvec(u), derivativevec(u_t);
  //// can't do this here because the PETSc Vec we get passed isn't
  //// in the right state!
  //for(uint i = 0; i < (*tsctx).bcs.size(); ++i) 
  //{
  //  (*(*tsctx).bcs[i]).apply(iteratedvec);
  //}
  (*(*(*tsctx).us_i).vector()) = iteratedvec;
  (*(*(*tsctx).us_dot).vector()) = derivativevec;
  // need to apply bcs here to ensure the residual gets assembled 
  // with the right bcs
  for(uint i = 0; i < (*tsctx).bcs.size(); ++i) 
  {
    (*(*tsctx).bcs[i]).apply(*(*(*tsctx).us_i).vector());
  }
  //// no point in doing this as we don't know what Vec we've been passed
  //// (local memory etc.?) - though it does mean that we come out
  //// at the end with the wrong solution on the boundaries and need to
  //// correct it poststep
  //iteratedvec = *(*(*tsctx).us_i).vector();

  dolfin::PETScVector rhs(f);
  dolfin::Assembler assembler;
  assembler.assemble(rhs, *(*tsctx).F);
  // apply the bcs non-linearly - this basically sets the change to 0 on the boundaries
  // which means we get the wrong solution at the end but everything is "right" instide the
  // solve (hopefully)
  for(uint i = 0; i < (*tsctx).bcs.size(); ++i)                     // loop over the bcs
  {
    (*(*tsctx).bcs[i]).apply(rhs, (*(*(*tsctx).us_i).vector()));
  }

  if ((*tsctx).print_norms)
  {
    PetscErrorCode perr;
    PetscReal norm;

    perr = VecNorm(u,NORM_2,&norm); CHKERRQ(perr);
    std::cout << "    FormIFunction: 2-norm u = " << norm << std::endl;

    perr = VecNorm(u,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormIFunction: inf-norm u = " << norm << std::endl;

    perr = VecNorm(u_t,NORM_2,&norm); CHKERRQ(perr);
    std::cout << "    FormIFunction: 2-norm u_t = " << norm << std::endl;

    perr = VecNorm(u_t,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormIFunction: inf-norm u_t = " << norm << std::endl;

    perr = VecNorm(f,NORM_2,&norm); CHKERRQ(perr);
    std::cout << "    FormIFunction: 2-norm f = " << norm << std::endl;

    perr = VecNorm(f,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormIFunction: inf-norm f = " << norm << std::endl;
  }
  
  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the residual function
//*******************************************************************|************************************************************//
PetscErrorCode FormRHSFunction(TS ts, PetscReal t, Vec u, Vec f, void* ctx) // petsc ts callback function to form the residual
{
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormRHSFunction" << std::endl;

  // update the time so that the bcs are applied at the right time
  // i.e. this function can be called on any time-level
  (*(*tsctx).t) = t;

  dolfin::PETScVector iteratedvec(u);
  //// can't do this here because the PETSc Vec we get passed isn't
  //// in the right state!
  //for(uint i = 0; i < (*tsctx).bcs.size(); ++i) 
  //{
  //  (*(*tsctx).bcs[i]).apply(iteratedvec);
  //}
  (*(*(*tsctx).us_i).vector()) = iteratedvec;
  // need to apply bcs here to ensure the residual gets assembled 
  // with the right bcs
  for(uint i = 0; i < (*tsctx).bcs.size(); ++i) 
  {
    (*(*tsctx).bcs[i]).apply(*(*(*tsctx).us_i).vector());
  }
  //// no point in doing this as we don't know what Vec we've been passed
  //// (local memory etc.?) - though it does mean that we come out
  //// at the end with the wrong solution on the boundaries and need to
  //// correct it poststep
  //iteratedvec = *(*(*tsctx).us_i).vector();

  update_functions(*tsctx);

  dolfin::PETScVector rhs(f);
  dolfin::Assembler assembler;
  assembler.assemble(rhs, *(*tsctx).G);
  // apply the bcs non-linearly - this basically sets the change to 0 on the boundaries
  // which means we get the wrong solution at the end but everything is "right" instide the
  // solve (hopefully)
  for(uint i = 0; i < (*tsctx).bcs.size(); ++i)                     // loop over the bcs
  {
    (*(*tsctx).bcs[i]).apply(rhs, (*(*(*tsctx).us_i).vector()));
  }

  if ((*tsctx).print_norms)
  {
    PetscErrorCode perr;
    PetscReal norm;

    perr = VecNorm(u,NORM_2,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSFunction: 2-norm u = " << norm << std::endl;

    perr = VecNorm(u,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSFunction: inf-norm u = " << norm << std::endl;

    perr = VecNorm(f,NORM_2,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSFunction: 2-norm f = " << norm << std::endl;

    perr = VecNorm(f,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSFunction: inf-norm f = " << norm << std::endl;
  }
  
  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the jacobian function
//*******************************************************************|************************************************************//
PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec u, Vec u_t, 
                                    PetscReal a, Mat A, Mat B,       // petsc ts callback function to form the jacobian
                                    void* ctx)
{
  Ctx *tsctx = (Ctx *)ctx;

  *(*tsctx).a = a;

  std::cout << "  In FormIJacobian" << std::endl;
  std::cout << "    a (shift) = " << a << std::endl;

  dolfin::PETScVector iteratedvec(u), derivativevec(u_t);
  (*(*(*tsctx).us_i).vector()) = iteratedvec;
  (*(*(*tsctx).us_dot).vector()) = derivativevec;

  dolfin::PETScMatrix matrix(A);
  dolfin::SystemAssembler assembler((*tsctx).JF, (*tsctx).F, (*tsctx).bcs);
  assembler.assemble(matrix);
  if ((*tsctx).ident_zeros)
  {
    matrix.ident_zeros();
  }

  if ((*tsctx).print_norms)
  {
    PetscErrorCode perr;
    PetscReal norm;

    perr = VecNorm(u,NORM_2,&norm); CHKERRQ(perr);
    std::cout << "    FormIJacobian: 2-norm u = " << norm << std::endl;

    perr = VecNorm(u,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormIJacobian: inf-norm u = " << norm << std::endl;

    perr = VecNorm(u_t,NORM_2,&norm); CHKERRQ(perr);
    std::cout << "    FormIJacobian: 2-norm u_t = " << norm << std::endl;

    perr = VecNorm(u_t,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormIJacobian: inf-norm u_t = " << norm << std::endl;

    perr = MatNorm(A,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    std::cout << "    FormIJacobian: Frobenius norm A = " << norm << std::endl;

    perr = MatNorm(A,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormIJacobian: inf-norm A = " << norm << std::endl;

    perr = MatNorm(B,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    std::cout << "    FormIJacobian: Frobenius norm B = " << norm << std::endl;

    perr = MatNorm(B,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormIJacobian: inf-norm B = " << norm << std::endl;
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the jacobian function
//*******************************************************************|************************************************************//
PetscErrorCode FormRHSJacobian(TS ts, PetscReal t, Vec u, Mat A, Mat B,       // petsc ts callback function to form the jacobian
                                    void* ctx)
{
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormRHSJacobian" << std::endl;

  dolfin::PETScVector iteratedvec(u);
  (*(*(*tsctx).us_i).vector()) = iteratedvec;

  update_functions(*tsctx);
  update_dfunctions(*tsctx);

  dolfin::PETScMatrix matrix(A);
  dolfin::SystemAssembler assembler((*tsctx).JG, (*tsctx).G, (*tsctx).bcs);
  assembler.assemble(matrix);
  if ((*tsctx).ident_zeros)
  {
    matrix.ident_zeros();
  }

  if ((*tsctx).print_norms)
  {
    PetscErrorCode perr;
    PetscReal norm;

    perr = VecNorm(u,NORM_2,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSJacobian: 2-norm u = " << norm << std::endl;

    perr = VecNorm(u,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSJacobian: inf-norm u = " << norm << std::endl;

    perr = MatNorm(A,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSJacobian: Frobenius norm A = " << norm << std::endl;

    perr = MatNorm(A,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSJacobian: inf-norm A = " << norm << std::endl;

    perr = MatNorm(B,NORM_FROBENIUS,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSJacobian: Frobenius norm B = " << norm << std::endl;

    perr = MatNorm(B,NORM_INFINITY,&norm); CHKERRQ(perr);
    std::cout << "    FormRHSJacobian: inf-norm B = " << norm << std::endl;
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that outputs diagnostics at the end of a timestep
//*******************************************************************|************************************************************//
PetscErrorCode TSCustomPreStep(TS ts)
{
  PetscErrorCode perr;

  Ctx *tsctx;
  perr = TSGetApplicationContext(ts, &tsctx); CHKERRQ(perr);

  //// can't do this here as the solution vector is not necessarily
  //// the correct bit of memory at this stage!
  //PetscReal pt;
  //perr = TSGetTime(ts, &pt); CHKERRQ(perr);
  //
  //Vec u;
  //perr = TSGetSolution(ts, &u); CHKERRQ(perr);
  //dolfin::PETScVector iteratedvec(u);
  //
  //for(uint i = 0; i < (*tsctx).bcs.size(); ++i) 
  //{
  //  (*(*tsctx).bcs[i]).apply(iteratedvec);
  //}

  if ((*tsctx).snes_vis_monitor)
  {
    // necessary for snes monitor
    PetscInt n;
    perr = TSGetStepNumber(ts, &n); CHKERRQ(perr);
    std::stringstream buffer;
    buffer.str(""); buffer << "ts-phasediagram_" << (*tsctx).nx << "_" << n << "_snes.xdmf";
    (*tsctx).file_snes.reset( new dolfin::XDMFFile((*(*tsctx).mesh).mpi_comm(), buffer.str()) );
    (*(*tsctx).file_snes).parameters["flush_output"] = true;
    (*(*tsctx).file_snes).parameters["functions_share_mesh"] = true;
    (*(*tsctx).file_snes).parameters["rewrite_function_mesh"] = false;
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that outputs diagnostics at the end of a timestep
//*******************************************************************|************************************************************//
PetscErrorCode TSCustomPreStage(TS ts, PetscReal t)
{
  PetscErrorCode perr;

  Ctx *tsctx;
  perr = TSGetApplicationContext(ts, &tsctx); CHKERRQ(perr);

  (*(*tsctx).t) = t;

  //// can't do this here as the solution vector is not necessarily
  //// the correct bit of memory at this stage and GetStages isn't
  //// implemented for all TS methods (it also seems to return the 
  //// wrong/unexpected thing for TSTHETA)
  //Vec *u;
  //PetscInt nstage;
  //perr = TSGetStages(ts, &nstage, &u); CHKERRQ(perr);
  //dolfin::PETScVector iteratedvec(u[nstage-1]);
  //*(*(*tsctx).us_i).vector() = iteratedvec;
  //
  //for(uint i = 0; i < (*tsctx).bcs.size(); ++i) 
  //{
  //  (*(*tsctx).bcs[i]).apply(iteratedvec);
  //  (*(*tsctx).bcs[i]).apply(*(*(*tsctx).us_i).vector());
  //}

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that outputs diagnostics at the end of a timestep
//*******************************************************************|************************************************************//
PetscErrorCode TSCustomPostStep(TS ts)
{
  PetscErrorCode perr, perr2;

  Ctx *tsctx;
  perr = TSGetApplicationContext(ts, &tsctx); CHKERRQ(perr);

  PetscReal pt;
  perr = TSGetTime(ts, &pt); CHKERRQ(perr);
  *(*tsctx).t = pt;

  Vec u;
  perr = TSGetSolution(ts, &u); CHKERRQ(perr);
  dolfin::PETScVector iteratedvec(u);
  // we apply the boundary conditions at the end because this is the 
  // only time we can be certain we're manipulating the n+1 timestep
  // solution
  for(uint i = 0; i < (*tsctx).bcs.size(); ++i) 
  {
    (*(*tsctx).bcs[i]).apply(iteratedvec);
  }
  (*(*(*tsctx).us_i).vector()) = iteratedvec;
  
  if ((*tsctx).file_u)
  {
    for (std::size_t i = 0; i < (*tsctx).mi_ptrs.size(); ++i) 
    {
      (*(*tsctx).file_u).write(*(*tsctx).mi_ptrs[i], *(*tsctx).t);
    }
    for (std::size_t i = 0; i < (*tsctx).cik_ptrs.size(); ++i) 
    {
      (*(*tsctx).file_u).write(*(*tsctx).cik_ptrs[i], *(*tsctx).t);
    }
  }

  PetscFunctionReturn(0);
}

//*******************************************************************|************************************************************//
// define the petsc snes monitor callback function that outputs a visualization file and a convergence file
//*******************************************************************|************************************************************//
PetscErrorCode SNESCustomMonitor(SNES snes, PetscInt n, PetscReal norm, void* ctx)
{
  std::cout << "  In SNESCustomMonitor" << std::endl;

  std::stringstream buffer;                                          // string buffer
  PetscErrorCode perr;                                               // petsc error code

  Ctx *snesctx = (Ctx *)ctx;

  Vec x;
  perr = SNESGetSolution(snes, &x);  CHKERRQ(perr);                  // get the solution vector from snes
  dolfin::PETScVector iteratedvec(x);
  (*(*(*snesctx).us_i).vector()) = iteratedvec;

  Vec rx;
  perr = SNESGetFunction(snes, &rx, PETSC_NULL, PETSC_NULL);         // get the residual function vector from snes
  CHKERRQ(perr);
  dolfin::PETScVector residualvec(rx);
  auto res = std::make_shared< dolfin::Function >((*snesctx).V);
  *(*res).vector() = residualvec;
  (*res).rename("Residual", "res");

  Vec dx;
  perr = SNESGetSolutionUpdate(snes, &dx);  CHKERRQ(perr);           // get the solution update vector from snes
  dolfin::PETScVector updatevec(dx);
  auto du = std::make_shared< dolfin::Function >((*snesctx).V);
  *(*du).vector() = updatevec;
  (*du).rename("Update", "du");

  assert((*snesctx).file_snes);
  for (std::size_t i = 0; i < (*snesctx).mi_ptrs.size(); ++i) 
  {
    (*(*snesctx).file_snes).write(*(*snesctx).mi_ptrs[i], (double) n);
  }
  for (std::size_t i = 0; i < (*snesctx).cik_ptrs.size(); ++i) 
  {
    (*(*snesctx).file_snes).write(*(*snesctx).cik_ptrs[i], (double) n);
  }
  (*(*snesctx).file_snes).write(*res, (double) n);
  (*(*snesctx).file_snes).write(*du, (double) n);

  PetscFunctionReturn(0);
}

PetscErrorCode SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu)
{                                                                    
                                                                     // do nothing
  PetscFunctionReturn(0);
} 


