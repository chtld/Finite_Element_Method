/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2006 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Ivan Christov, Wolfgang Bangerth, Texas A&M University, 2006
 */
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
namespace Step25
{
  using namespace dealii;
  template <int dim>
  class SineGordonProblem
  {
  public:
    SineGordonProblem();
    void run();
  private:
    void         make_grid_and_dofs();
    void         assemble_system();
    void         compute_nl_term(const Vector<double> &old_data,
                                 const Vector<double> &new_data,
                                 Vector<double> &      nl_term) const;
    void         compute_nl_matrix(const Vector<double> &old_data,
                                   const Vector<double> &new_data,
                                   SparseMatrix<double> &nl_matrix) const;
    unsigned int solve();
    void         output_results(const unsigned int timestep_number) const;
    Triangulation<dim> triangulation;
    FE_Q<dim>          fe;
    DoFHandler<dim>    dof_handler;
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    const unsigned int n_global_refinements;
    double       time;
    const double final_time, time_step;
    const double theta;
    Vector<double> solution, solution_update, old_solution;
    Vector<double> M_x_velocity;
    Vector<double> system_rhs;
    const unsigned int output_timestep_skip;
  };
  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution(const unsigned int n_components = 1, const double time = 0.)
      : Function<dim>(n_components, time)
    {}
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };
  template <int dim>
  double ExactSolution<dim>::value(const Point<dim> &p,
                                   const unsigned int /*component*/) const
  {
    double t = this->get_time();
    switch (dim)
      {
        case 1:
          {
            const double m  = 0.5;
            const double c1 = 0.;
            const double c2 = 0.;
            return -4. * std::atan(m / std::sqrt(1. - m * m) *
                                   std::sin(std::sqrt(1. - m * m) * t + c2) /
                                   std::cosh(m * p[0] + c1));
          }
        case 2:
          {
            const double theta  = numbers::PI / 4.;
            const double lambda = 1.;
            const double a0     = 1.;
            const double s      = 1.;
            const double arg    = p[0] * std::cos(theta) +
                               std::sin(theta) * (p[1] * std::cosh(lambda) +
                                                  t * std::sinh(lambda));
            return 4. * std::atan(a0 * std::exp(s * arg));
          }
        case 3:
          {
            const double theta = numbers::PI / 4;
            const double phi   = numbers::PI / 4;
            const double tau   = 1.;
            const double c0    = 1.;
            const double s     = 1.;
            const double arg   = p[0] * std::cos(theta) +
                               p[1] * std::sin(theta) * std::cos(phi) +
                               std::sin(theta) * std::sin(phi) *
                                 (p[2] * std::cosh(tau) + t * std::sinh(tau));
            return 4. * std::atan(c0 * std::exp(s * arg));
          }
        default:
          Assert(false, ExcNotImplemented());
          return -1e8;
      }
  }
  template <int dim>
  class InitialValues : public Function<dim>
  {
  public:
    InitialValues(const unsigned int n_components = 1, const double time = 0.)
      : Function<dim>(n_components, time)
    {}
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };
  template <int dim>
  double InitialValues<dim>::value(const Point<dim> & p,
                                   const unsigned int component) const
  {
    return ExactSolution<dim>(1, this->get_time()).value(p, component);
  }
  template <int dim>
  SineGordonProblem<dim>::SineGordonProblem()
    : fe(1)
    , dof_handler(triangulation)
    , n_global_refinements(6)
    , time(-5.4414)
    , final_time(2.7207)
    , time_step(10 * 1. / std::pow(2., 1. * n_global_refinements))
    , theta(0.5)
    , output_timestep_skip(1)
  {}
  template <int dim>
  void SineGordonProblem<dim>::make_grid_and_dofs()
  {
    GridGenerator::hyper_cube(triangulation, -10, 10);
    triangulation.refine_global(n_global_refinements);
    std::cout << "   Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "   Total number of cells: " << triangulation.n_cells()
              << std::endl;
    dof_handler.distribute_dofs(fe);
    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);
    system_matrix.reinit(sparsity_pattern);
    mass_matrix.reinit(sparsity_pattern);
    laplace_matrix.reinit(sparsity_pattern);
    MatrixCreator::create_mass_matrix(dof_handler,
                                      QGauss<dim>(fe.degree + 1),
                                      mass_matrix);
    MatrixCreator::create_laplace_matrix(dof_handler,
                                         QGauss<dim>(fe.degree + 1),
                                         laplace_matrix);
    solution.reinit(dof_handler.n_dofs());
    solution_update.reinit(dof_handler.n_dofs());
    old_solution.reinit(dof_handler.n_dofs());
    M_x_velocity.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }
  template <int dim>
  void SineGordonProblem<dim>::assemble_system()
  {
    system_matrix.copy_from(mass_matrix);
    system_matrix.add(std::pow(time_step * theta, 2), laplace_matrix);
    SparseMatrix<double> tmp_matrix(sparsity_pattern);
    compute_nl_matrix(old_solution, solution, tmp_matrix);
    system_matrix.add(std::pow(time_step * theta, 2), tmp_matrix);
    system_rhs = 0.;
    Vector<double> tmp_vector(solution.size());
    mass_matrix.vmult(system_rhs, solution);
    laplace_matrix.vmult(tmp_vector, solution);
    system_rhs.add(std::pow(time_step * theta, 2), tmp_vector);
    mass_matrix.vmult(tmp_vector, old_solution);
    system_rhs.add(-1.0, tmp_vector);
    laplace_matrix.vmult(tmp_vector, old_solution);
    system_rhs.add(std::pow(time_step, 2) * theta * (1 - theta), tmp_vector);
    system_rhs.add(-time_step, M_x_velocity);
    compute_nl_term(old_solution, solution, tmp_vector);
    system_rhs.add(std::pow(time_step, 2) * theta, tmp_vector);
    system_rhs *= -1.;
  }
  template <int dim>
  void SineGordonProblem<dim>::compute_nl_term(const Vector<double> &old_data,
                                               const Vector<double> &new_data,
                                               Vector<double> &nl_term) const
  {
    nl_term = 0;
    const QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_values | update_JxW_values |
                              update_quadrature_points);
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    Vector<double>                       local_nl_term(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<double>                  old_data_values(n_q_points);
    std::vector<double>                  new_data_values(n_q_points);
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        local_nl_term = 0;
        fe_values.reinit(cell);
        fe_values.get_function_values(old_data, old_data_values);
        fe_values.get_function_values(new_data, new_data_values);
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            local_nl_term(i) +=
              (std::sin(theta * new_data_values[q_point] +
                        (1 - theta) * old_data_values[q_point]) *
               fe_values.shape_value(i, q_point) * fe_values.JxW(q_point));
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          nl_term(local_dof_indices[i]) += local_nl_term(i);
      }
  }
  template <int dim>
  void SineGordonProblem<dim>::compute_nl_matrix(
    const Vector<double> &old_data,
    const Vector<double> &new_data,
    SparseMatrix<double> &nl_matrix) const
  {
    QGauss<dim>   quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_JxW_values |
                              update_quadrature_points);
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();
    FullMatrix<double> local_nl_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    std::vector<double>                  old_data_values(n_q_points);
    std::vector<double>                  new_data_values(n_q_points);
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        local_nl_matrix = 0;
        fe_values.reinit(cell);
        fe_values.get_function_values(old_data, old_data_values);
        fe_values.get_function_values(new_data, new_data_values);
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              local_nl_matrix(i, j) +=
                (std::cos(theta * new_data_values[q_point] +
                          (1 - theta) * old_data_values[q_point]) *
                 fe_values.shape_value(i, q_point) *
                 fe_values.shape_value(j, q_point) * fe_values.JxW(q_point));
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            nl_matrix.add(local_dof_indices[i],
                          local_dof_indices[j],
                          local_nl_matrix(i, j));
      }
  }
  template <int dim>
  unsigned int SineGordonProblem<dim>::solve()
  {
    SolverControl solver_control(1000, 1e-12 * system_rhs.l2_norm());
    SolverCG<>    cg(solver_control);
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
    cg.solve(system_matrix, solution_update, system_rhs, preconditioner);
    return solver_control.last_step();
  }
  template <int dim>
  void SineGordonProblem<dim>::output_results(
    const unsigned int timestep_number) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u");
    data_out.build_patches();
    const std::string filename =
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtk";
    std::ofstream output(filename);
    data_out.write_vtk(output);
  }
  template <int dim>
  void SineGordonProblem<dim>::run()
  {
    make_grid_and_dofs();
    {
      AffineConstraints<double> constraints;
      constraints.close();
      VectorTools::project(dof_handler,
                           constraints,
                           QGauss<dim>(fe.degree + 1),
                           InitialValues<dim>(1, time),
                           solution);
    }
    output_results(0);
    unsigned int timestep_number = 1;
    for (time += time_step; time <= final_time;
         time += time_step, ++timestep_number)
      {
        old_solution = solution;
        std::cout << std::endl
                  << "Time step #" << timestep_number << "; "
                  << "advancing to t = " << time << "." << std::endl;
        double initial_rhs_norm = 0.;
        bool   first_iteration  = true;
        do
          {
            assemble_system();
            if (first_iteration == true)
              initial_rhs_norm = system_rhs.l2_norm();
            const unsigned int n_iterations = solve();
            solution += solution_update;
            if (first_iteration == true)
              std::cout << "    " << n_iterations;
            else
              std::cout << '+' << n_iterations;
            first_iteration = false;
          }
        while (system_rhs.l2_norm() > 1e-6 * initial_rhs_norm);
        std::cout << " CG iterations per nonlinear step." << std::endl;
        Vector<double> tmp_vector(solution.size());
        laplace_matrix.vmult(tmp_vector, solution);
        M_x_velocity.add(-time_step * theta, tmp_vector);
        laplace_matrix.vmult(tmp_vector, old_solution);
        M_x_velocity.add(-time_step * (1 - theta), tmp_vector);
        compute_nl_term(old_solution, solution, tmp_vector);
        M_x_velocity.add(-time_step, tmp_vector);
        if (timestep_number % output_timestep_skip == 0)
          output_results(timestep_number);
      }
  }
} // namespace Step25
int main()
{
  try
    {
      using namespace dealii;
      using namespace Step25;
      SineGordonProblem<1> sg_problem;
      sg_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
