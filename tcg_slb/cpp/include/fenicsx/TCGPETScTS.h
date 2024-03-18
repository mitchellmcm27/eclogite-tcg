#include "petsc.h"

// a structure used to pass bucket data into SNES and TS callback functions
typedef struct {
  TCGCoefficients* tcg;
  std::shared_ptr<dolfinx::fem::Form<type>> JF, JG, F, G;
  std::shared_ptr<dolfinx::fem::Function<type>> us_i, us_dot;
  std::shared_ptr<dolfinx::fem::Constant<type>> a;
  std::vector<double> t;
  std::vector<std::vector<double> > y;
  bool print_norms;
} Ctx;

// some function headers for PETSc
PetscErrorCode FormIFunction(TS ts, PetscReal tt, Vec u, Vec u_t, 
                                    Vec f,  void* ctx);              // petsc ts callback function to form the residual

PetscErrorCode FormRHSFunction(TS ts, PetscReal tt, Vec u,
                                    Vec f,  void* ctx);              // petsc ts callback function to form the residual

PetscErrorCode FormIJacobian(TS ts, PetscReal tt, Vec u, Vec u_t, 
                                    PetscReal a, Mat A, Mat B,       // petsc ts callback function to form the jacobian
                                    void* ctx);

PetscErrorCode FormRHSJacobian(TS ts, PetscReal tt, Vec u,
                                    Mat A, Mat B,                    // petsc ts callback function to form the jacobian
                                    void* ctx);

//extern PetscErrorCode TSCustomPreStep(TS ts);
//
//extern PetscErrorCode TSCustomPreStage(TS ts, PetscReal tt);
//
//extern PetscErrorCode TSCustomPostStep(TS ts);
//
//extern PetscErrorCode SNESCustomMonitor(SNES snes, PetscInt its, PetscReal norm, void* ctx);

PetscErrorCode SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu);

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

  perr = TSCreate(MPI_COMM_WORLD, &ts); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetApplicationContext(ts, &ctx); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetType(ts, TSBDF); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetTime(ts, 0.0); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetTimeStep(ts, 1.e-6); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetMaxTime(ts, 1.0); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetMaxSteps(ts, 100000000); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetMaxSNESFailures(ts, -1); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSSetTolerances(ts, 1.e-9, PETSC_NULL, 1.e-5, PETSC_NULL); CHKERRABORT(MPI_COMM_WORLD, perr);
  
  perr = TSGetSNES(ts, &snes); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = SNESSetTolerances(snes, 1.e-6, 1.e-6, PETSC_DEFAULT, 10, PETSC_DEFAULT); CHKERRABORT(MPI_COMM_WORLD, perr);

  PetscViewerAndFormat *svf;
  PetscViewer sviewer;
  perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ts),&sviewer); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = PetscViewerAndFormatCreate(sviewer,PETSC_VIEWER_DEFAULT,&svf); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = SNESMonitorSet(snes, (PetscErrorCode (*)(SNES,PetscInt,PetscReal,void*))SNESMonitorDefault, 
                                           svf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); CHKERRABORT(MPI_COMM_WORLD, perr);

  //if (ctx.snes_vis_monitor)
  //{
  //  perr = SNESMonitorSet(snes, SNESCustomMonitor,                // set a custom snes monitor
  //                                        &ctx, PETSC_NULL); 
  //}

  perr = SNESSetType(snes, SNESVINEWTONRSLS); CHKERRABORT(MPI_COMM_WORLD, perr);

  perr = SNESGetKSP(snes, &ksp); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = KSPSetType(ksp, "preonly"); CHKERRABORT(MPI_COMM_WORLD, perr);

  perr = KSPGetPC(ksp, &pc); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = PCSetType(pc, "lu"); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = PCFactorSetMatSolverType(pc, "umfpack"); CHKERRABORT(MPI_COMM_WORLD, perr);

  perr = TSSetFromOptions(ts); CHKERRABORT(MPI_COMM_WORLD, perr);

  Mat JFmat = dolfinx::fem::petsc::create_matrix(*ctx.JF);
  perr = TSSetIJacobian(ts, JFmat, JFmat, FormIJacobian, (void *) &ctx); CHKERRABORT(MPI_COMM_WORLD, perr);
  
  Mat JGmat = dolfinx::fem::petsc::create_matrix(*ctx.JG);
  perr = TSSetRHSJacobian(ts, JGmat, JGmat, FormRHSJacobian, (void *) &ctx); CHKERRABORT(MPI_COMM_WORLD, perr);

  Vec Fvec = dolfinx::la::petsc::create_vector(*ctx.F->function_spaces()[0]->dofmap()->index_map,
                                               ctx.F->function_spaces()[0]->dofmap()->index_map_bs());
  perr = TSSetIFunction(ts, Fvec, FormIFunction, (void *) &ctx); CHKERRABORT(MPI_COMM_WORLD, perr);

  Vec Gvec = dolfinx::la::petsc::create_vector(*ctx.G->function_spaces()[0]->dofmap()->index_map,
                                               ctx.G->function_spaces()[0]->dofmap()->index_map_bs());
  perr = TSSetRHSFunction(ts, Gvec, FormRHSFunction, (void *) &ctx); CHKERRABORT(MPI_COMM_WORLD, perr);

  PetscViewerAndFormat *vf;
  PetscViewer viewer;
  perr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ts),&viewer); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = PetscViewerAndFormatCreate(viewer,PETSC_VIEWER_DEFAULT,&vf); CHKERRABORT(MPI_COMM_WORLD, perr);
  perr = TSMonitorSet(ts, (PetscErrorCode (*)(TS,PetscInt,PetscReal,Vec,void*))TSMonitorDefault, 
                                           vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy); CHKERRABORT(MPI_COMM_WORLD, perr);

  //perr = TSSetPreStep(ts, TSCustomPreStep); CHKERRABORT(MPI_COMM_WORLD, perr);
  //perr = TSSetPreStage(ts, TSCustomPreStage); CHKERRABORT(MPI_COMM_WORLD, perr);
  //perr = TSSetPostStep(ts, TSCustomPostStep); CHKERRABORT(MPI_COMM_WORLD, perr);
  ////perr = TSSetSolution(ts, (*std::dynamic_pointer_cast<dolfin::PETScVector>((*ctx.us_i).vector())).vec()); CHKERRABORT(MPI_COMM_WORLD, perr);

  SNESType snestype;
  perr = SNESGetType(snes,&snestype); CHKERRABORT(MPI_COMM_WORLD, perr);

  if (strcmp(snestype, SNESVINEWTONRSLS) || strcmp(snestype, SNESVINEWTONSSLS))
  {
    std::cout << "applying bounds for snestype: " << snestype << std::endl;
    Vec ub;
    perr = VecDuplicate(Fvec, &ub); CHKERRABORT(PETSC_COMM_WORLD,perr);
    perr = VecSet(ub, 1.0); CHKERRABORT(PETSC_COMM_WORLD,perr);
    Vec lb;
    perr = VecDuplicate(Fvec, &lb); CHKERRABORT(PETSC_COMM_WORLD,perr);
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

//*******************************************************************|************************************************************//
// run the problem using TS
//*******************************************************************|************************************************************//
void run_ts(TS &ts, Ctx &ctx)
{
  PetscErrorCode perr, perr2;

  Vec work = dolfinx::la::petsc::create_vector_wrap(*(ctx.us_i->x()));
  perr = TSSolve(ts, work); 
  TSConvergedReason reason;
  perr2 = TSGetConvergedReason(ts, &reason); 
  std::cout << "  TSConvergedReason: " << reason << std::endl;
  CHKERRABORT(MPI_COMM_WORLD, perr);
  CHKERRABORT(MPI_COMM_WORLD, perr2);
}

//*******************************************************************|************************************************************//
// define the petsc ts callback function that assembles the residual function
//*******************************************************************|************************************************************//
PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec u, Vec u_t, Vec f, void* ctx) // petsc ts callback function to form the residual
{
  PetscErrorCode perr;
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormIFunction" << std::endl;

  PetscInt n;
  perr = VecGetLocalSize(u, &n);

  const type *uu, *uu_t;
  perr = VecGetArrayRead(u, &uu);
  perr = VecGetArrayRead(u_t, &uu_t);
  xtl::span<const type> xu(uu, n);
  xtl::span<const type> xu_t(uu_t, n);
  std::copy(xu.begin(), xu.end(), (*tsctx).us_i->x()->mutable_array().begin());
  std::copy(xu_t.begin(), xu_t.end(), (*tsctx).us_dot->x()->mutable_array().begin());
  perr = VecRestoreArrayRead(u, &uu);
  perr = VecRestoreArrayRead(u_t, &uu_t);
  
  perr = VecZeroEntries(f);
  type *ff;
  perr = VecGetArray(f, &ff);
  xtl::span<type> xf(ff, n);
  dolfinx::fem::assemble_vector<type>(xf, *(*tsctx).F);
  perr = VecRestoreArray(f, &ff);

  if ((*tsctx).print_norms)
  {
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
  PetscErrorCode perr;
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormRHSFunction" << std::endl;

  tsctx->tcg->interpolate_coeffs();

  PetscInt n;
  perr = VecGetLocalSize(u, &n);

  const type *uu;
  perr = VecGetArrayRead(u, &uu);
  xtl::span<const type> xu(uu, n);
  std::copy(xu.begin(), xu.end(), (*tsctx).us_i->x()->mutable_array().begin());
  perr = VecRestoreArrayRead(u, &uu);
  
  perr = VecZeroEntries(f);
  type *ff;
  perr = VecGetArray(f, &ff);
  xtl::span<type> xf(ff, n);
  dolfinx::fem::assemble_vector<type>(xf, *(*tsctx).G);
  perr = VecRestoreArray(f, &ff);

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
  PetscErrorCode perr;
  Ctx *tsctx = (Ctx *)ctx;

  (*tsctx).a->value[0] = a;

  std::cout << "  In FormIJacobian" << std::endl;
  std::cout << "    a (shift) = " << a << std::endl;

  PetscInt n;
  perr = VecGetLocalSize(u, &n);

  const type *uu, *uu_t;
  perr = VecGetArrayRead(u, &uu);
  perr = VecGetArrayRead(u_t, &uu_t);
  xtl::span<const type> xu(uu, n);
  xtl::span<const type> xu_t(uu_t, n);
  std::copy(xu.begin(), xu.end(), (*tsctx).us_i->x()->mutable_array().begin());
  std::copy(xu_t.begin(), xu_t.end(), (*tsctx).us_dot->x()->mutable_array().begin());
  perr = VecRestoreArrayRead(u, &uu);
  perr = VecRestoreArrayRead(u_t, &uu_t);
  
  MatZeroEntries(A);
  dolfinx::fem::assemble_matrix(dolfinx::la::petsc::Matrix::set_block_fn(A, ADD_VALUES), *(*tsctx).JF, {});
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  if ((*tsctx).print_norms)
  {
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
  PetscErrorCode perr;
  Ctx *tsctx = (Ctx *)ctx;

  std::cout << "  In FormRHSJacobian" << std::endl;

  tsctx->tcg->interpolate_all_coeffs();

  PetscInt n;
  perr = VecGetLocalSize(u, &n);

  const type *uu, *uu_t;
  perr = VecGetArrayRead(u, &uu);
  xtl::span<const type> xu(uu, n);
  std::copy(xu.begin(), xu.end(), (*tsctx).us_i->x()->mutable_array().begin());
  perr = VecRestoreArrayRead(u, &uu);
  
  MatZeroEntries(A);
  dolfinx::fem::assemble_matrix(dolfinx::la::petsc::Matrix::set_block_fn(A, ADD_VALUES), *(*tsctx).JG, {});
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  if ((*tsctx).print_norms)
  {
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

////*******************************************************************|************************************************************//
//// define the petsc ts monitor callback function
////*******************************************************************|************************************************************//
//PetscErrorCode TSCustomMonitor(TS ts, PetscInt i, PetscReal tt, Vec u, void* ctx)
//{
//  PetscErrorCode perr, perr2;
//  Ctx *tsctx = (Ctx *)ctx;
//
//  const PetscScalar *uu;
//
//  std::size_t& I = (*(*tsctx).tcg).I;
//  std::size_t& K = (*(*tsctx).tcg).K;
//  std::vector<std::size_t>& Kis = (*(*tsctx).tcg).Kis;
//
//  perr = VecGetArrayRead(u, &uu); CHKERRABORT(PETSC_COMM_WORLD, perr);
//
//  (*tsctx).t.push_back(tt);
//  for (std::size_t i=0; i < I+K; ++i)
//  {
//    (*tsctx).y[i].push_back(uu[i]);
//  }
//
//  perr = VecRestoreArrayRead(u, &uu); CHKERRABORT(PETSC_COMM_WORLD, perr);
//
//  PetscFunctionReturn(0);
//}
//
PetscErrorCode SNESVIDummyComputeVariableBounds(SNES snes, Vec xl, Vec xu)
{                                                                    
  PetscFunctionReturn(0);
} 

