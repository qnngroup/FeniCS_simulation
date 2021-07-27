
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

emitter_coords = [Point(-x_gap/2, 0), Point(-x_gap/2-h, w/2), Point(-x_gap/2-h-l, w/2), Point(-x_gap/2-h-l, -w/2), Point(-x_gap/2-h, -w/2)]
emitter = Polygon(emitter_coords)
collector_coords = [Point(x_gap/2, 0), Point(x_gap/2+h, -w/2), Point(x_gap/2+h+l, -w/2), Point(x_gap/2+h+l, w/2), Point(x_gap/2+h, w/2)]
collector = Polygon(collector_coords)
domain = Rectangle(Point(-x_gap/2-h-l-1,-w/2-3), Point(x_gap/2+h+l+1,w/2+3)) - emitter - collector
#domain = Rectangle(Point(-x_gap/2-h-l+1,-w/2-3), Point(x_gap/2+h+l+1,w/2+3))

#domain.set_subdomain(1, emitter)
#domain.set_subdomain(2, collector)
mesh = generate_mesh(domain, 20)

V = FunctionSpace(mesh, 'CG', 1)

# Define boundary condition
def out_boundary(x, on_boundary):
    return on_boundary and (x[0] <= -x_gap/2-h-l-1+tol or x[0] >= x_gap/2+h+l+1-tol or x[1]<= -w/2-3 +tol or x[1]>= w/2+3-tol)

def bemitter(x, on_boundary):
    return on_boundary and x[0] <=0 and not out_boundary(x, on_boundary)

def bcollector(x, on_boundary):
    return on_boundary and x[0] >=0 and not out_boundary(x, on_boundary)

out = DirichletBC(V, Constant(0), out_boundary)
em = DirichletBC(V, Constant(0), bemitter)
col = DirichletBC(V, Constant(30), bcollector)

bcs = [out, em, col]
# Define subdomain markers and integration measure
# markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1, mesh.domains())
# dx = Measure('dx', domain=mesh, subdomain_data=markers)

materials = MeshFunction('size_t', mesh, mesh.topology().dim()-1, mesh.domains())

u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = Constant(0.0)*v*dx
A, b = assemble_system(a, L, bcs)
u = Function(V)
U = u.vector()

# Compute solution
solve(A, U, b)

  # Save solution in VTK format
f = File("problem.pvd")
f << u
u_array = np.array(U)
print(u_array.shape, u_array)
  # Plot solution

meshfile = File('poisson/mesh.pvd')
meshfile << mesh


# f= Constant(0.0)
# a = dot(grad(u), grad(v))*dx
# L = f*v*dx

# u = Function(V)
# solve(a == L, u, bcs)
# outfile = File('poisson/e-potential.pvd')
# outfile << u

plt.figure()
plot(u)
plt.savefig("potential_sol_pol.png")