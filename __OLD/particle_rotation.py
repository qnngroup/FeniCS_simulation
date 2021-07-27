from leopart import (
    particles,
    advect_rk3,
    PDEStaticCondensation,
    FormsPDEMap,
    RandomCircle,
    SlottedDisk,
    assign_particle_values,
    l2projection,
)

from dolfin import (
Expression,
Point,
VectorFunctionSpace,
Mesh,
Constant,
FunctionSpace,
assemble,
dx,
refine,
Function,
assign,
DirichletBC,
)

from mpi4py import MPI as pyMPI
import numpy as np

comm = pyMPI.COMM_WORLD

(x0, y0) = (0.0, 0.0)
(xc, yc) = (-0.15, 0.0)
(r, rdisk) = (0.5, 0.2)
rwidth = 0.05
(lb, ub) = (0.0, 1.0)

# Mesh
mesh = Mesh("circle_0.xml")
mesh = refine(refine(refine(mesh)))
store_step = 5

psi0_expr = SlottedDisk(radius=r, center=[xc, yc], width=rwidth, depth=0.0, degree=3, lb=lb, ub=ub)

Tend = 2.0
dt = Constant(0.02)
num_steps = np.rint(Tend / float(dt))

W = FunctionSpace(mesh, "DG", 1)
T = FunctionSpace(mesh, "DG", 0)
Wbar = FunctionSpace(mesh, "DGT", 1)

bc = DirichletBC(Wbar, Constant(0.0), "on_boundary")

(psi_h, psi_h0, psi_h00) = (Function(W), Function(W), Function(W))
psibar_h = Function(Wbar)

V = VectorFunctionSpace(mesh, "DG", 3)
uh = Function(V)
uh.assign(Expression(("-Uh*x[1]", "Uh*x[0]"), Uh=np.pi, degree=3))

x = RandomCircle(Point(x0, y0), r).generate([750, 750])
s = assign_particle_values(x, psi0_expr)
p = particles(x, [s], mesh)

ap = advect_rk3(p, V, uh, "closed")
lstsq_psi = l2projection(p, W, 1)

lstsq_psi.project(psi_h0, lb, ub)
assign(psi_h00, psi_h0)

step = 0
t = 0.0
area_0 = assemble(psi_h0 * dx)
(psi_h_min, psi_h_max) = (0.0, 0.0)

while step < num_steps:
    step += 1
    t += float(dt)

    if comm.Get_rank() == 0:
        print(f"Step {step}")

    ap.do_step(float(dt))

    lstsq_psi.project(psi_h, lb, ub)
    psi_h_min = min(psi_h_min, psi_h.vector().min())
    psi_h_max = max(psi_h_max, psi_h.vector().max())
    if comm.rank == 0:
        print("Min max phi {} {}".format(psi_h_min, psi_h_max))


    assign(psi_h0, psi_h)

area_end = assemble(psi_h * dx)
num_part = p.number_of_particles()
l2_error = np.sqrt(abs(assemble((psi_h00 - psi_h) * (psi_h00 - psi_h) * dx)))
if comm.Get_rank() == 0:
    print(f"Num particles {num_part}")
    print(f"Area error {abs(area_end - area_0)}")
    print(f"L2-Error {l2_error}")