#include <dolfinx.h>
#include <dolfinx/geometry/utils.h>
#include <dolfinx/fem/petsc.h>
#include <basix/e-lagrange.h>
#include "reactions.h"
#include "petsc.h"
#include "PhaseDiagramTS.h"
#include "TCGCoefficients.h"
#include "TCGPETScTS.h"

//*******************************************************************|************************************************************//
// main
//*******************************************************************|************************************************************//
int main(int argc, char* argv[])
{
  dolfinx::common::subsystem::init_logging(argc, argv);
  dolfinx::common::subsystem::init_petsc(argc, argv);

  Ctx ctx;

  MgFeSiO4_all_slb_rx rxn;
  TCGReaction tcgrxn(&rxn);

  PetscErrorCode perr;
  //perr = PetscInitialize(&argc,&argv,(char*)0,(char*)0);

  // Create mesh
  PetscInt pnx = 1;
  perr = PetscOptionsGetInt(NULL,NULL,"-nx",&pnx,NULL);
  std::size_t nx = pnx;
  auto mesh = std::make_shared<dolfinx::mesh::Mesh>(dolfinx::mesh::create_interval(MPI_COMM_WORLD, 
                                                                                   nx, {0.0,1.0}, 
                                                                                   dolfinx::mesh::GhostMode::none));

  PetscBool plog = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL,"-l",&plog,NULL);

  // log files
  if (plog)
  {
    std::ostringstream debug_file;
    debug_file << "ts-log." << nx;
    if(std::freopen(debug_file.str().c_str(), "w", stdout) == NULL)
    {
      std::cout << "Failed to redirect stdio for debugging.";
    }
  }
  std::stringstream buffer;
  buffer.str(""); buffer << "ts-phasediagram_" << nx << ".xdmf";

  // do this immediately after the log file has been opened (if requested)
  perr = PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(MPI_COMM_WORLD, perr);

  basix::FiniteElement P0_e = basix::element::create_lagrange(dolfinx::mesh::cell_type_to_basix_type(dolfinx::mesh::CellType::interval), 
                                                              0, basix::element::lagrange_variant::equispaced, true);
  auto P0_V = std::make_shared<dolfinx::fem::FunctionSpace>(dolfinx::fem::create_functionspace(mesh, P0_e, 1));

  PetscReal pDa = 1.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-Da",&pDa,NULL);
  auto Da = std::make_shared<dolfinx::fem::Constant<type>>(pDa);

  PetscReal pT = 1673.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-T",&pT,NULL);
  auto T = std::make_shared<dolfinx::fem::Function<type>>(P0_V);
  T->x()->set(pT);

  PetscReal pp = 150000.0;
  perr = PetscOptionsGetReal(NULL, NULL,"-p",&pp,NULL);
  auto p = std::make_shared<dolfinx::fem::Function<type>>(P0_V);
  p->x()->set(pp);

  PetscReal peps = 1.e-2;
  perr = PetscOptionsGetReal(NULL, NULL,"-eps",&peps,NULL);
  auto eps = std::make_shared<dolfinx::fem::Constant<type>>(peps);

  PetscInt pi0 = 0;
  perr = PetscOptionsGetInt(NULL, NULL,"-i0",&pi0,NULL);
  std::size_t i0 = pi0;

  PetscReal pci0k0 = 0.9;
  perr = PetscOptionsGetReal(NULL, NULL,"-ci0k0",&pci0k0,NULL);
  double ci0k0 = pci0k0;

  PetscBool ppnorms = PETSC_FALSE;
  perr = PetscOptionsGetBool(NULL, NULL, "-print_norms", &ppnorms, NULL);
  ctx.print_norms = ppnorms;

  auto V = std::make_shared<dolfinx::fem::FunctionSpace>(
                      dolfinx::fem::create_functionspace(functionspace_form_PhaseDiagramTS_JG,
                                                         "us_i", mesh));
  ctx.us_i.reset(new dolfinx::fem::Function<type>(V));
  ctx.us_i->name = "us_i";

  for (std::size_t i = 0; i < tcgrxn.I; ++i)
  {
    ctx.us_i->sub(i).interpolate(
                             [i,i0](auto&& x)->xt::xtensor<type, 2>
                             { 
                               xt::xtensor<type, 2> ui0({1,x.shape(1)});
                               (i==i0) ? ui0.fill(1.0) : ui0.fill(0.0);
                               return ui0;
                             }
                            );
  }

  for (std::size_t i = 0; i < tcgrxn.I; ++i)
  {
    ctx.us_i->sub(tcgrxn.I+i).interpolate(
                             [tcgrxn,i,ci0k0](auto&& x)->xt::xtensor<type, 2>
                             { 
                               xt::xtensor<type, 2> ui0({tcgrxn.Kis[i],x.shape(1)});
                               (tcgrxn.Kis[i]==1) ? xt::row(ui0,0).fill(1.0) : xt::row(ui0,0).fill(ci0k0);
                               for (std::size_t k=1; k<tcgrxn.Kis[i]; ++k) xt::row(ui0,k).fill((1.-ci0k0)/(tcgrxn.Kis[i]-1.));
                               return ui0;
                             }
                            );
  }

  std::cout << "us_i = ";
  for (std::size_t i = 0; i < ctx.us_i->x()->array().size(); ++i)
  {
    std::cout << ctx.us_i->x()->array()[i] << " ";
  }
  std::cout << std::endl;

  std::vector<std::shared_ptr<dolfinx::fem::Function<type>>> mi_ptrs, cik_ptrs;
  for (std::size_t i = 0; i < tcgrxn.I; ++i) mi_ptrs.push_back(std::make_shared<dolfinx::fem::Function<type>>(ctx.us_i->sub(i)));
  for (std::size_t i = 0; i < tcgrxn.I; ++i) cik_ptrs.push_back(std::make_shared<dolfinx::fem::Function<type>>(ctx.us_i->sub(tcgrxn.I+i)));

  // coefficients
  auto V_rhoi = std::make_shared<dolfinx::fem::FunctionSpace>(
                           dolfinx::fem::create_functionspace(functionspace_form_PhaseDiagramTS_JG,
                                                              "rhoi_i", mesh));
  auto rhoi_i = std::make_shared<dolfinx::fem::Function<type>>(V_rhoi);
  rhoi_i->name = "rhoi";

  auto V_drhoidus = std::make_shared<dolfinx::fem::FunctionSpace>(
                           dolfinx::fem::create_functionspace(functionspace_form_PhaseDiagramTS_JG,
                                                              "drhoidus_i", mesh));
  auto drhoidus_i = std::make_shared<dolfinx::fem::Function<type>>(V_drhoidus);
  drhoidus_i->name = "drhoidus";

  auto V_Gammai = std::make_shared<dolfinx::fem::FunctionSpace>(
                           dolfinx::fem::create_functionspace(functionspace_form_PhaseDiagramTS_JG,
                                                              "Gammai_i", mesh));
  auto Gammai_i = std::make_shared<dolfinx::fem::Function<type>>(V_Gammai);
  Gammai_i->name = "Gammai";

  auto V_dGammaidus = std::make_shared<dolfinx::fem::FunctionSpace>(
                           dolfinx::fem::create_functionspace(functionspace_form_PhaseDiagramTS_JG,
                                                              "dGammaidus_i", mesh));
  auto dGammaidus_i = std::make_shared<dolfinx::fem::Function<type>>(V_dGammaidus);
  dGammaidus_i->name = "dGammaidus";

  auto V_Gammaik = std::make_shared<dolfinx::fem::FunctionSpace>(
                           dolfinx::fem::create_functionspace(functionspace_form_PhaseDiagramTS_JG,
                                                              "Gammaik_i", mesh));
  auto Gammaik_i = std::make_shared<dolfinx::fem::Function<type>>(V_Gammaik);
  Gammaik_i->name = "Gammaik";

  auto V_dGammaikdus = std::make_shared<dolfinx::fem::FunctionSpace>(
                           dolfinx::fem::create_functionspace(functionspace_form_PhaseDiagramTS_JG,
                                                              "dGammaikdus_i", mesh));
  auto dGammaikdus_i = std::make_shared<dolfinx::fem::Function<type>>(V_dGammaikdus);
  dGammaikdus_i->name = "dGammaikdus";

  ctx.tcg = new TCGCoefficients(mi_ptrs, cik_ptrs, 
                                T, p, eps, 
                                rhoi_i, drhoidus_i, 
                                Gammai_i, dGammaidus_i, 
                                Gammaik_i, dGammaikdus_i,
                                tcgrxn);

  ctx.us_dot.reset(new dolfinx::fem::Function<type>(V));
  ctx.us_dot->name = "us_dot";
  ctx.us_dot->x()->set(0.0);

  ctx.a.reset( new dolfinx::fem::Constant<type>(666.e6) );
  auto rho0 = std::make_shared<dolfinx::fem::Constant<type>>(ctx.tcg->rho({{0.5},{0.5},{0.5}})(0));
  std::cout << "rho0 = " << rho0->value[0] << std::endl;
  std::cout << "Da = " << Da->value[0] << std::endl;

  ctx.tcg->interpolate_all_coeffs();

  std::map<std::string, std::shared_ptr<const dolfinx::fem::Function<type>>> 
     coefficients = {{"us_i", ctx.us_i},
                     {"us_dot", ctx.us_dot},
                     {"rhoi_i", rhoi_i},
                     {"drhoidus_i", drhoidus_i},
                     {"Gammai_i", Gammai_i},
                     {"dGammaidus_i", dGammaidus_i},
                     {"Gammaik_i", Gammaik_i},
                     {"dGammaikdus_i", dGammaikdus_i}};

  std::map<std::string, std::shared_ptr<const dolfinx::fem::Constant<type>>> 
     constants = {{"rho0", rho0},
                  {"Da", Da},
                  {"eps", eps},
                  {"a", ctx.a}};

  ctx.JF.reset( new dolfinx::fem::Form<type>(dolfinx::fem::create_form<type>(
                         *form_PhaseDiagramTS_JF, {V,V}, coefficients, constants, {})) );

  ctx.F.reset( new dolfinx::fem::Form<type>(dolfinx::fem::create_form<type>(
                        *form_PhaseDiagramTS_F, {V}, coefficients, constants, {})) );

  ctx.JG.reset( new dolfinx::fem::Form<type>(dolfinx::fem::create_form<type>(
                         *form_PhaseDiagramTS_JG, {V,V}, coefficients, constants, {})) );

  ctx.G.reset( new dolfinx::fem::Form<type>(dolfinx::fem::create_form<type>(
                        *form_PhaseDiagramTS_G, {V}, coefficients, constants, {})) );


  TS ts = setup_ts(ctx);

  run_ts(ts, ctx);

  dolfinx::common::subsystem::finalize_petsc();
  return 0;
}


