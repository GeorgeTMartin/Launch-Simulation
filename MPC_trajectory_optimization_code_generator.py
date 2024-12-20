from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, builders
from MPC_rocket_model import rocket_model
import scipy.linalg
import numpy as np

def acados_generator(rebuild_flag):
    Tf = 600  # prediction horizon
    N = 600  # number of discretization steps

    # create render arguments
    ocp = AcadosOcp()

    # export model
    model, constraint = rocket_model()

    # define acados ODE
    model_ac = AcadosModel()
    model_ac.f_impl_expr = model.f_impl_expr
    model_ac.f_expl_expr = model.f_expl_expr
    model_ac.x = model.x
    model_ac.xdot = model.xdot
    model_ac.u = model.u
    model_ac.z = model.z
    model_ac.p = model.p
    model_ac.name = model.name
    ocp.model = model_ac

    # define constraint
    model_ac.con_h_expr = constraint.expr

    # set dimensions
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx

    ocp.solver_options.N_horizon = N
    ns = 1
    nsh = 1

    # set costs
    Q = np.diag([0, 0, 0, 0, 0])      # Should be of Size = nx
    R = np.eye(nu)                    # Should be of Size = nu 
    R[0, 0] = 0                       # Should be of Size = nu 
    Qe = np.diag([0, 1e1, 0, 1e1, 0]) # Should be of Size = nx

    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"
    unscale = N / Tf

    ocp.cost.W = unscale * scipy.linalg.block_diag(Q, R)
    ocp.cost.W_e = Qe / unscale

    Vx = np.zeros((ny, nx))
    Vx[:nx, :nx] = np.eye(nx)
    ocp.cost.Vx = Vx

    Vu = np.zeros((ny, nu))
    Vu[5, 0] = 0.1
    ocp.cost.Vu = Vu

    Vx_e = np.zeros((ny_e, nx))
    Vx_e[:nx, :nx] = np.eye(nx)
    ocp.cost.Vx_e = Vx_e

    ocp.cost.zl = 100 * np.ones((ns,))
    ocp.cost.Zl = 0.1 * np.ones((ns,))
    ocp.cost.zu = 100 * np.ones((ns,))
    ocp.cost.Zu = 0.1 * np.ones((ns,))

    # set initial references
    ocp.cost.yref = np.array([0, 0, 0, 0, 0, 0])       # Size: nx + nu
    ocp.cost.yref_e = np.array([0, 1e1, 0, 1e1, 0])        # Size: nx 

    # setting constraints
    ocp.constraints.lbu = np.array([model.theta_min])
    ocp.constraints.ubu = np.array([model.theta_max])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.lh = np.array([model.theta_min])
    ocp.constraints.uh = np.array([model.theta_max])
    ocp.constraints.lsh = np.zeros(nsh)
    ocp.constraints.ush = np.ones(nsh)*0.01
    ocp.constraints.idxsh = np.array([0])

    # set intial condition
    ocp.constraints.x0 = model.x0

    # set QP solver and integration
    ocp.solver_options.tf = Tf
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    # ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    # ocp.solver_options.nlp_solver_type = "SQP_RTI"
    ocp.solver_options.nlp_solver_max_iter = 5
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "ERK"
    ocp.solver_options.sim_method_num_stages = 4
    ocp.solver_options.sim_method_num_steps = 1


    ocp.solver_options.qp_solver_tol_stat = 1e-2
    ocp.solver_options.qp_solver_tol_eq = 1e-2
    ocp.solver_options.qp_solver_tol_ineq = 1e-2
    ocp.solver_options.qp_solver_tol_comp = 1e-2

    cmake_template = builders.CMakeBuilder()
    # cmake_template.build_targets = "GENERIC"
    cmake_template.options_on = ["BUILD_ACADOS_OCP_SOLVER_LIB"]

    # create solver
    if rebuild_flag == False:
        acados_solver = AcadosOcpSolver(ocp, build=False, generate=False, json_file="acados_ocp.json", cmake_builder=cmake_template)
    else:
        acados_solver = AcadosOcpSolver(ocp, json_file="acados_ocp.json", cmake_builder=cmake_template)
    return constraint, model, acados_solver

if __name__ == "__main__":
    acados_generator(rebuild_flag=True)