from dolfin import *
import mshr
import matplotlib.pyplot as plt

# Domain
domain = mshr.Rectangle(Point(0,0), Point(4,4))\
          - mshr.Rectangle(Point(1,1), Point(3,3))

# Mesh
mesh = mshr.generate_mesh(domain, 20)

class lower_edge(SubDomain):
  def inside(self, x, on_boundary):
    return True if x[0]<= 3.0 and x[0]>= 1.0 and x[1] == 1.0 else False

class top_edge(SubDomain):
  def inside(self, x, on_boundary):
    return True if x[0]<= 3.0 and x[0]>= 1.0 and x[1] == 3.0 else False

class left_edge(SubDomain):
  def inside(self, x, on_boundary):
    return True if x[0]==1.0 and x[1]<= 3.0 and x[1]>= 1.0 else False

class right_edge(SubDomain):
  def inside(self, x, on_boundary):
    return True if x[0]==3.0 and x[1]<= 3.0 and x[1]>= 1.0 else False

lower_edge = lower_edge()
top_edge = top_edge()
left_edge = left_edge()
right_edge = right_edge()

# Subdomain marker
mf = MeshFunction("size_t", mesh, 2)
mf.set_all(0)
lower_edge.mark(mf, 1)
top_edge.mark(mf, 2)
left_edge.mark(mf, 3)
right_edge.mark(mf, 4)
plot(mesh)
plt.show()
# Define facet function over the new mesh


# Extract boundary mesh
bmesh = BoundaryMesh(mesh, "exterior", True)
#plot(bmesh, interactive=True)

t = 0.0
for i in range(10):
  t += 0.1
  for x in bmesh.coordinates():

    if lower_edge.inside(x, True):
      x[0] += 0.0
      x[1] += 0.05

    if top_edge.inside(x, True):
      x[0] += 0.0
      x[1] += 0.05

    if left_edge.inside(x, True):
      x[0] += 0.0
      x[1] += 0.05

    if right_edge.inside(x, True):
      x[0] += 0.0
      x[1] += 0.05

  ALE.move(mesh, bmesh)
plot(mesh)
plt.show()