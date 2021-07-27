
from __future__ import print_function
from fenics import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np


x_gap = 0.5
y_gap = 1
h = 1
w = 2
l =2

tol = 1e-8
INITIAL_GUESS_ZERO = True

emitter_coords = [Point(-x_gap/2, 0), Point(-x_gap/2-h, w/2), Point(-x_gap/2-h-l, w/2), Point(-x_gap/2-h-l, -w/2), Point(-x_gap/2-h, -w/2)]
emitter = Polygon(emitter_coords)
collector_coords = [Point(x_gap/2, 0), Point(x_gap/2+h, -w/2), Point(x_gap/2+h+l, -w/2), Point(x_gap/2+h+l, w/2), Point(x_gap/2+h, w/2)]
collector = Polygon(collector_coords)
domain = Rectangle(Point(-x_gap/2-h-l-1,-w/2-3), Point(x_gap/2+h+l+1,w/2+3)) - emitter - collector

mesh = generate_mesh(domain, 20)

V = FunctionSpace(mesh, 'CG', 1)

# Define boundary condition
def out_boundary(x, on_boundary):
    return on_boundary and (x[0] <= -x_gap/2-h-l-1+tol or x[0] >= x_gap/2+h+l+1-tol or x[1]<= -w/2-3 +tol or x[1]>= w/2+3-tol)

def bemitter(x, on_boundary):
    return on_boundary and x[0] <=0 and not out_boundary(x, on_boundary)

def bcollector(x, on_boundary):
    return on_boundary and x[0] >=0 and not out_boundary(x, on_boundary)

meshfile = File('poisson/mesh.pvd')
meshfile << mesh

u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = Constant(0.0)*v*dx

collector_voltages = np.linspace(0, 10, 2)
for i, collector_voltage in enumerate(collector_voltages):

  out = DirichletBC(V, Constant(0), out_boundary)
  em = DirichletBC(V, Constant(0), bemitter)
  col = DirichletBC(V, Constant(collector_voltage), bcollector)

  bcs = [out, em, col]
 
  A, b = assemble_system(a, L, bcs)
  solver = KrylovSolver('cg', 'ilu')
  u = Function(V)
  U = u.vector()
  solver.solve(A, U, b)

  # # Compute solution
  # if INITIAL_GUESS_ZERO:
  #   u = Function(V)
  #   U = u.vector()
  #   solver.solve(A, U, b)
  #   INITIAL_GUESS_ZERO = False
  # else:
  #   solver.solve(A, U, b, solver_parameters={'nonzero_initial_guess': True})



  # Save solution in VTK format
  f = File("voltage_sweep/no{}_Vcoll{}.pvd".format(i, collector_voltage))
  f << u


  plt.figure()
  plot(u)
  plt.savefig("voltage_sweep/no{}_Vcoll{}.png".format(i, collector_voltage))
  plt.close()