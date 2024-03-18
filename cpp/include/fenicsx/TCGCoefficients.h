#include <dolfinx.h>
#include "Reaction.h"

using type = PetscScalar;

class TCGReaction
{

public:
  Reaction* rxn;

  std::vector<double> mi;
  std::vector<std::vector<double>> cik;

  std::size_t I, K, J;
  std::vector<std::size_t> Kis;

  TCGReaction(Reaction *rxn_in)
  {
    rxn = rxn_in;

    cik = rxn->zero_C();
    I = cik.size();
    Kis.resize(I,1);
    mi.resize(I,0.);
    for (std::size_t i=0; i < I; ++i)
    {
      Kis[i] = cik[i].size();
    }
    K = std::accumulate(Kis.begin(), Kis.end(), 0.0);
  }

};


//*******************************************************************|************************************************************//
// a class returning the output from a reaction in a useful format
//*******************************************************************|************************************************************//
class TCGCoefficients : public TCGReaction
{

public:
  std::vector<std::shared_ptr<dolfinx::fem::Function<type>>> mi_ptrs, cik_ptrs;

  std::shared_ptr<dolfinx::fem::Function<type>> T_ptr, p_ptr;
  std::shared_ptr<dolfinx::fem::Constant<type>> eps_ptr;

  std::shared_ptr<dolfinx::fem::Function<type>> rhoi_ptr, drhoidus_ptr, Gammai_ptr, dGammaidus_ptr, Gammaik_ptr, dGammaikdus_ptr;

  double eps, T, p;

  xt::xtensor<type, 2> T_values, p_values;
  std::vector<xt::xtensor<type, 2>> mi_values, cik_values;

  std::shared_ptr<const dolfinx::mesh::Mesh> mesh;
  std::shared_ptr<const dolfinx::fem::FiniteElement> element;

  xt::xtensor<double, 2> default_x;
  std::vector<std::int32_t> default_cells;

  bool params_evald;

  //*****************************************************************|************************************************************//
  // constructor
  //*****************************************************************|************************************************************//
  TCGCoefficients(std::vector<std::shared_ptr<dolfinx::fem::Function<type>>> mi_ptrs_in, 
                  std::vector<std::shared_ptr<dolfinx::fem::Function<type>>> cik_ptrs_in,
                  std::shared_ptr<dolfinx::fem::Function<type>> T_ptr_in, 
                  std::shared_ptr<dolfinx::fem::Function<type>> p_ptr_in, 
                  std::shared_ptr<dolfinx::fem::Constant<type>> eps_ptr_in,
                  std::shared_ptr<dolfinx::fem::Function<type>> rhoi_ptr_in,
                  std::shared_ptr<dolfinx::fem::Function<type>> drhoidus_ptr_in,
                  std::shared_ptr<dolfinx::fem::Function<type>> Gammai_ptr_in,
                  std::shared_ptr<dolfinx::fem::Function<type>> dGammaidus_ptr_in,
                  std::shared_ptr<dolfinx::fem::Function<type>> Gammaik_ptr_in,
                  std::shared_ptr<dolfinx::fem::Function<type>> dGammaikdus_ptr_in,
                  TCGReaction &tcgrxn) : TCGReaction(tcgrxn),
                  mi_ptrs(mi_ptrs_in), cik_ptrs(cik_ptrs_in), 
                  T_ptr(T_ptr_in), p_ptr(p_ptr_in), eps_ptr(eps_ptr_in),
                  rhoi_ptr(rhoi_ptr_in), drhoidus_ptr(drhoidus_ptr_in), 
                  Gammai_ptr(Gammai_ptr_in), dGammaidus_ptr(dGammaidus_ptr_in), 
                  Gammaik_ptr(Gammaik_ptr_in), dGammaikdus_ptr(dGammaikdus_ptr_in),
                  params_evald(false)
  {
    mi_values.resize(I);
    cik_values.resize(I);

    mesh = rhoi_ptr->function_space()->mesh();
    element = rhoi_ptr->function_space()->element();

    const int tdim = mesh->topology().dim();
    std::int32_t num_cells = mesh->topology().index_map(tdim)->size_local() + mesh->topology().index_map(tdim)->num_ghosts();

    // Get mesh geometry data and the element coordinate map
    const std::size_t gdim = mesh->geometry().dim();
    const dolfinx::graph::AdjacencyList<std::int32_t>& x_dofmap = mesh->geometry().dofmap();
    const std::size_t num_dofs_g = x_dofmap.num_links(0);
    xtl::span<const double> x_g = mesh->geometry().x();
    
    const dolfinx::fem::CoordinateElement& cmap = mesh->geometry().cmap();
    
    // Get the interpolation points on the reference cells
    const xt::xtensor<double, 2>& X = element->interpolation_points();
    const xt::xtensor<double, 2> phi = xt::view(cmap.tabulate(0, X), 0, xt::all(), xt::all(), 0);
      
    // Push reference coordinates (X) forward to the physical coordinates (x) for each cell
    xt::xtensor<double, 2> coordinate_dofs = xt::zeros<double>({num_dofs_g, gdim});
    std::array<std::size_t, 2> shape = {3, num_cells * X.shape(0)};

    default_x.resize(shape);
    default_x.fill(0.0);
    default_cells.resize(num_cells * X.shape(0));
    for (std::int32_t c = 0; c < num_cells; ++c)
    { 
      // Get geometry data for current cell
      auto x_dofs = x_dofmap.links(c);
      for (std::size_t i = 0; i < x_dofs.size(); ++i)
      {
        std::copy_n(std::next(x_g.begin(), 3 * x_dofs[i]), gdim,
                    std::next(coordinate_dofs.begin(), i * gdim));
      }

      // Push forward coordinates (X -> x)
      for (std::size_t p = 0; p < X.shape(0); ++p)
      {
        for (std::size_t j = 0; j < gdim; ++j)
        {
          double acc = 0;
          for (std::size_t k = 0; k < num_dofs_g; ++k)
            acc += phi(p, k) * coordinate_dofs(k, j);
          default_x(j, c * X.shape(0) + p) = acc;
        }
        default_cells[c * X.shape(0) + p] = c;
      }
    }

  }

  //*
  void interpolate_all_coeffs()
  {
    interpolate_coeffs();
    params_evald = true;
    interpolate_derivative_coeffs();
  }

  //*
  void interpolate_coeffs()
  {
    eval_params(default_x);
    params_evald = true;
    rhoi_ptr->interpolate(rhoi());
    Gammai_ptr->interpolate(Gammai());
    Gammaik_ptr->interpolate(Gammaik());

    //std::cout << "rhoi_ptr = ";
    //for (std::size_t i = 0; i < rhoi_ptr->x()->array().size(); ++i)
    //{
    //  std::cout << rhoi_ptr->x()->array()[i] << " ";
    //}
    //std::cout << std::endl;

    //std::cout << "Gammai_ptr = ";
    //for (std::size_t i = 0; i < Gammai_ptr->x()->array().size(); ++i)
    //{
    //  std::cout << Gammai_ptr->x()->array()[i] << " ";
    //}
    //std::cout << std::endl;

    //std::cout << "Gammaik_ptr = ";
    //for (std::size_t i = 0; i < Gammaik_ptr->x()->array().size(); ++i)
    //{
    //  std::cout << Gammaik_ptr->x()->array()[i] << " ";
    //}
    //std::cout << std::endl;

    params_evald = false;
  }

  void interpolate_derivative_coeffs()
  {
    eval_params(default_x);
    params_evald = true;
    drhoidus_ptr->interpolate(drhoidus());
    dGammaidus_ptr->interpolate(dGammaidus());
    dGammaikdus_ptr->interpolate(dGammaikdus());

    //std::cout << "drhoidus_ptr = ";
    //for (std::size_t i = 0; i < drhoidus_ptr->x()->array().size(); ++i)
    //{
    //  std::cout << drhoidus_ptr->x()->array()[i] << " ";
    //}
    //std::cout << std::endl;

    //std::cout << "dGammaidus_ptr = ";
    //for (std::size_t i = 0; i < dGammaidus_ptr->x()->array().size(); ++i)
    //{
    //  std::cout << dGammaidus_ptr->x()->array()[i] << " ";
    //}
    //std::cout << std::endl;

    //std::cout << "dGammaikdus_ptr = ";
    //for (std::size_t i = 0; i < dGammaikdus_ptr->x()->array().size(); ++i)
    //{
    //  std::cout << dGammaikdus_ptr->x()->array()[i] << " ";
    //}
    //std::cout << std::endl;

    params_evald = false;
  }

  //*
  std::vector<std::int32_t> get_interpolation_cells(const xt::xtensor<double, 2>& x)
  {
    // assume that x coordinates will match if shape matches!
    if ((x.shape(0) == default_x.shape(0)) and (x.shape(1) == default_x.shape(1)))
    {
      return default_cells;
    }
    else
    {
      const int tdim = mesh->topology().dim();
      std::int32_t num_cells = mesh->topology().index_map(tdim)->size_local() + mesh->topology().index_map(tdim)->num_ghosts();
      std::vector<std::int32_t> cell_indices(num_cells);
      std::iota(cell_indices.begin(), cell_indices.end(), 0);
      dolfinx::geometry::BoundingBoxTree midpoint_tree = 
              dolfinx::geometry::create_midpoint_tree(*mesh, tdim, 
                                                      xtl::span<const std::int32_t>(cell_indices.data(), cell_indices.size()));
      dolfinx::geometry::BoundingBoxTree tree(*mesh, tdim);
      return dolfinx::geometry::compute_closest_entity(tree, midpoint_tree, *mesh, xt::transpose(x));
    }
  }

  void eval_params(const xt::xtensor<double, 2>& x)
  {
    if (!params_evald)
    {
      std::vector<std::int32_t> cells = get_interpolation_cells(x);
      xtl::span<const std::int32_t> xcells(cells.data(), cells.size());
      xt::xtensor<double, 2> xT = xt::transpose(x);
      
      eval_eps();
      eval_T(xT, xcells);
      eval_p(xT, xcells);
      eval_mi(xT, xcells);
      eval_cik(xT, xcells);
    }
  }

  void eval_eps()
  {
    eps = eps_ptr->value[0];
  }

  void eval_T(const xt::xtensor<double, 2>& xT, const xtl::span<const std::int32_t>& cells)
  {
    T_values.resize({xT.shape(0), 1});
    T_ptr->eval(xT, cells, T_values);
  }

  void eval_p(const xt::xtensor<double, 2>& xT, const xtl::span<const std::int32_t>& cells)
  {
    p_values.resize({xT.shape(0), 1});
    p_ptr->eval(xT, cells, p_values);
  }

  void eval_mi(const xt::xtensor<double, 2>& xT, const xtl::span<const std::int32_t>& cells)
  {
    for (std::size_t i=0; i < I; ++i)
    {
      mi_values[i].resize({xT.shape(0), 1});
      mi_ptrs[i]->eval(xT, cells, mi_values[i]);
    }
  }

  void eval_cik(const xt::xtensor<double, 2>& xT, const xtl::span<const std::int32_t>& cells)
  {
    for (std::size_t i=0; i < I; ++i)
    {
      cik_values[i].resize({xT.shape(0), Kis[i]});
      cik_ptrs[i]->eval(xT, cells, cik_values[i]);
    }
  }

  void update_params_x(std::size_t& xi)
  {
    update_T_x(xi);
    update_p_x(xi);
    update_mi_x(xi);
    update_cik_x(xi);
  }

  void update_T_x(std::size_t& xi)
  {
    T = T_values.at(xi, 0);
  }
  
  void update_p_x(std::size_t& xi)
  {
    p = p_values.at(xi, 0);
  }

  void update_mi_x(std::size_t& xi)
  {
    for (std::size_t i=0; i < I; i++)
    {
      mi[i] = mi_values[i].at(xi, 0);
    }

    //std::cout << "mi = ";
    //for (std::size_t i = 0; i < I; ++i)
    //{
    //  std::cout << mi[i] << " ";
    //}
    //std::cout << std::endl;
  }
  
  void update_cik_x(std::size_t& xi)
  {
    double ctot = 0.0;
    for (std::size_t i=0; i < I; i++)
    {
      for (std::size_t k=0; k < Kis[i]; k++)
      {
        cik[i][k] = std::min(std::max(cik_values[i].at(xi,k), eps), 1.0-eps);
      }
      ctot = std::accumulate(cik[i].begin(), cik[i].end(), 0.0);
      for (std::size_t k=0; k < Kis[i]; k++)
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

  //*****************************************************************|************************************************************//
  // return bulk density
  //*****************************************************************|************************************************************//
  xt::xarray<type> rho(const xt::xtensor<double, 2> &x)
  {
    xt::xarray<type> xvals(x.shape(1));
    xt::xarray<type> work = rhoi()(x);
    for (std::size_t xi=0; xi < x.shape(1); ++xi)
    {
      update_params_x(xi);
      auto workx = xt::col(work,xi);
      std::transform(mi.begin(), mi.end(), workx.begin(), workx.begin(), std::divides<double>());
      xvals(xi) = 1./std::accumulate(workx.begin(), workx.end(), 0.0);
    }
    return xvals;
  }

  //*****************************************************************|************************************************************//
  // return phase densities
  //*****************************************************************|************************************************************//
  const std::function<xt::xarray<type>(const xt::xtensor<double, 2>&)> rhoi()
  {
    return [&](const xt::xtensor<double, 2> &x)->xt::xarray<type>
    {
      eval_params(x);

      std::vector<double> values(I,0.0);
      xt::xarray<type>::shape_type shape = {I, x.shape(1)};
      xt::xarray<type> xvals(shape);

      for (std::size_t xi=0; xi < x.shape(1); ++xi)
      {
        update_params_x(xi);

        //Evaluate
        rxn->rho(T,p,cik,values);

        std::copy(values.begin(), values.end(), xt::col(xvals,xi).begin());

        //std::cout << "rhoi: ";
        //for (std::size_t i=0; i < I; ++i) std::cout << values[i] << " ";
        //std::cout << "\n";
      }

      return xvals;
    };
  }

  //*****************************************************************|************************************************************//
  // return phase density derivatives
  //*****************************************************************|************************************************************//
  const std::function<xt::xarray<type>(const xt::xtensor<double, 2>&)> drhoidus()
  {
    return [&](const xt::xtensor<double, 2> &x)->xt::xarray<type>
    {
      eval_params(x);

      xt::xarray<type>::shape_type shape = {I*(I+K), x.shape(1)};
      xt::xarray<type> xvals(shape);
      xvals.fill(0.);

      std::size_t sKi;
      for (std::size_t xi=0; xi < x.shape(1); ++xi)
      {
        update_params_x(xi);

        //Evaluate
        sKi = 0;
        for (std::size_t i=0; i < I; ++i)
        {
          std::vector<double> drhoidcik = rxn->phases()[i]->drho_dc(T,p,cik[i]);
          std::copy(drhoidcik.begin(), drhoidcik.end(), xt::col(xvals,xi).begin() + i*(I+K) + I + sKi);
          sKi += Kis[i];
        }

        //std::cout << "drhoidus: ";
        //for (std::size_t i=0; i < I; ++i)
        //{
        //  for (std::size_t d=0; d < I + K; ++d) std::cout << values[i][d] << " ";
        //}
        //std::cout << "\n";
      }

      return xvals;
    };
  }

  //*****************************************************************|************************************************************//
  // return phase sources
  //*****************************************************************|************************************************************//
  const std::function<xt::xarray<type>(const xt::xtensor<double, 2>&)> Gammai()
  {
    return [&](const xt::xtensor<double, 2> &x)->xt::xarray<type>
    {
      eval_params(x);

      std::vector<double> values(I,0.0);
      xt::xarray<type>::shape_type shape = {I, x.shape(1)};
      xt::xarray<type> xvals(shape);

      for (std::size_t xi=0; xi < x.shape(1); ++xi)
      {
        update_params_x(xi);

        //Evaluate
        rxn->Gamma_i(T,p,cik,mi,values);

        std::copy(values.begin(), values.end(), xt::col(xvals,xi).begin());

        //std::cout << "Gammai: ";
        //for (std::size_t i=0; i < I; ++i) std::cout << values[i] << " ";
        //std::cout << "\n";
      }

      return xvals;
    };
  }

  //*****************************************************************|************************************************************//
  // return phase source derivatives
  //*****************************************************************|************************************************************//
  const std::function<xt::xarray<type>(const xt::xtensor<double, 2>&)> dGammaidus()
  {
    return [&](const xt::xtensor<double, 2> &x)->xt::xarray<type>
    {
      eval_params(x);

      std::vector<std::vector<double> > dGammaidmi(I);
      std::vector<std::vector<std::vector<double> > > dGammaidcik(I);

      xt::xarray<type>::shape_type shape = {I*(I+K), x.shape(1)};
      xt::xarray<type> xvals(shape);
      xvals.fill(0.);
    
      std::size_t sKi;
      for (std::size_t xi=0; xi < x.shape(1); ++xi)
      {
        update_params_x(xi);

        //Evaluate
        rxn->dGamma_i_dPhi(T,p,cik,mi,dGammaidmi);
        rxn->dGamma_i_dC(T,p,cik,mi,dGammaidcik);
        for (std::size_t i=0; i < I; ++i)
        {
          std::copy(dGammaidmi[i].begin(), dGammaidmi[i].end(), xt::col(xvals,xi).begin() + i*(I+K));
          sKi = I;
          for (std::size_t di=0; di < I; di++)
          {
            std::copy(dGammaidcik[i][di].begin(), dGammaidcik[i][di].end(), xt::col(xvals,xi).begin() + i*(I+K) + sKi);
            sKi += Kis[di];
          }
        }

        //std::cout << "dGammaidus: ";
        //for (std::size_t i=0; i < I; ++i)
        //{
        //  for (std::size_t d=0; d < I + K; ++d) std::cout << values[i][d] << " ";
        //}
        //std::cout << "\n";

      }

      return xvals;
    };
  }

  //*****************************************************************|************************************************************//
  // return component sources
  //*****************************************************************|************************************************************//
  const std::function<xt::xarray<type>(const xt::xtensor<double, 2>&)> Gammaik()
  {
    return [&](const xt::xtensor<double, 2> &x)->xt::xarray<type>
    {
      eval_params(x);

      std::vector<std::vector<double> > Gammaik(I);

      xt::xarray<type>::shape_type shape = {K, x.shape(1)};
      xt::xarray<type> xvals(shape);

      std::size_t sdKi;
      for (std::size_t xi=0; xi < x.shape(1); ++xi)
      {
        update_params_x(xi);

        rxn->Gamma_ik(T,p,cik,mi,Gammaik);
        sdKi = 0;
        for (std::size_t i=0; i < I; ++i)
        {
          std::copy(Gammaik[i].begin(), Gammaik[i].end(), xt::col(xvals,xi).begin() + sdKi);
          sdKi += Kis[i];
        }

        //std::cout << "Gammaik: ";
        //for (std::size_t i=0; i < K; ++i) std::cout << values[i] << " ";
        //std::cout << "\n";

      }

      return xvals;
    };
  }

  //*****************************************************************|************************************************************//
  // return component source derivatives
  //*****************************************************************|************************************************************//
  const std::function<xt::xarray<type>(const xt::xtensor<double, 2>&)> dGammaikdus()
  {
    return [&](const xt::xtensor<double, 2> &x)->xt::xarray<type>
    {
      eval_params(x);

      xt::xarray<type>::shape_type shape = {K*(I+K), x.shape(1)};
      xt::xarray<type> xvals(shape);
      xvals.fill(0.);
    
      std::vector<std::vector<double> > dGammaikdmii;
      double dGammaikdcikik;

      std::size_t sKi, sdKi;
      for (std::size_t xi=0; xi < x.shape(1); ++xi)
      {
        update_params_x(xi);

        sKi = 0;
        sdKi = 0;
        for (std::size_t i=0; i<I; ++i)
        {
          rxn->dGamma_ik_dPhi(T,p,cik,mi,i,dGammaikdmii);
          for (std::size_t k=0; k<Kis[i]; ++k)
          {
            std::copy(dGammaikdmii[k].begin(), dGammaikdmii[k].end(), xt::col(xvals,xi).begin() + (sKi + k)*(I + K));
            sdKi = I;
            for (std::size_t di=0; di<I; ++di)
            {
              for (std::size_t dk=0; dk<Kis[di]; ++dk)
              {
                dGammaikdcikik = rxn->dGamma_ik_dC(T,p,cik,mi,i,k,di,dk);
                xt::col(xvals,xi).at((sKi + k)*(I + K) + sdKi + dk) = dGammaikdcikik;
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
      }

      return xvals;
    };
  }

};


