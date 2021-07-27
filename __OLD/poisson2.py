
from __future__ import print_function
from fenics import *
from mshr import *
import matplotlib.pyplot as plt


x_gap = 0.5
y_gap = 1
h = 1
w = 2
l =2

emitter_coords = [Point(-x_gap/2, 0), Point(-x_gap/2-h, w/2), Point(-x_gap/2-h-l, w/2), Point(-x_gap/2-h-l, -w/2), Point(-x_gap/2-h, -w/2)]
emitter = Polygon(emitter_coords)
collector_coords = [Point(x_gap/2, 0), Point(x_gap/2+h, -w/2), Point(x_gap/2+h+l, -w/2), Point(x_gap/2+h+l, w/2), Point(x_gap/2+h, w/2)]
collector = Polygon(collector_coords)
domain = Rectangle(Point(-x_gap/2-h-l,-w/2-3), Point(x_gap/2+h+l,w/2+3))

domain.set_subdomain(1, emitter)
domain.set_subdomain(2, collector)
mesh = generate_mesh(domain, 32)

V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
bc1 = DirichletBC(V, Constant(0), 'on_boundary')
class Emitter(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class collector(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

     
# Define subdomain markers and integration measure
markers = MeshFunction('size_t', mesh, 2, mesh.domains())
dx = Measure('dx', domain=mesh, subdomain_data=markers)
plt.figure()
plot(markers)
plt.savefig("electrodes.png")  

plt.close()
plt.figure()
plot(mesh)
plt.savefig("oldmesh.png")

print(markers)
meshfile = File('poisson/mesh.pvd')
meshfile << mesh

plot(mesh)
plt.show()
