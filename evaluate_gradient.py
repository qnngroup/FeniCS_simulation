import scipy.interpolate as spint
import numpy as np
import matplotlib.pyplot as plt
from Geometry_emitter import Geometry_emitter
from fenics import *
from mshr import *

def intersect_surface(outline_points, point_star, tol):
    x_star, y_star = point_star
    is_on = False
    for i in range(len(outline_points[:, 0]) -1):
        if abs(outline_points[i, 0] - outline_points[i+1, 0]) < tol:
            if abs(x_star - outline_points[i, 0]) < tol and ((y_star <= outline_points[i,1] and y_star >= outline_points[i+1,1]) or (y_star >= outline_points[i,1] and y_star <= outline_points[i+1,1])):
                is_on = True
        elif abs(outline_points[i, 1] - outline_points[i+1, 1]) < tol:
            if abs(y_star - outline_points[i, 1]) < tol and ((x_star <= outline_points[i,0] and x_star >= outline_points[i+1,0]) or (x_star >= outline_points[i,0] and x_star <= outline_points[i+1,0])):
                is_on = True
        elif (x_star <= outline_points[i,0] and x_star >= outline_points[i+1,0]) or (x_star >= outline_points[i,0] and x_star <= outline_points[i+1,0]):
            m = (outline_points[i,1] - outline_points[i+1,1])/(outline_points[i,0] - outline_points[i+1,0])
            b = outline_points[i,1] - m*outline_points[i,0]
            if abs(y_star - m*x_star - b) < tol:
                is_on = True
    return is_on

x_gap = 20                 # gap between emitter and collector in the x direction 
y_gap = 50                # gap between 2 gates in the y direction
w = 100                    # width of emitter and collector
wg = 80                    # width of gates
l = 100                     # length of emitter
lg = 100                   # length of gate
h = 50                     # height of triangle
hg = 50                     # height of gate triangle 

emitter = np.array([[-x_gap/2, 0], [-x_gap/2-h, w/2], [-x_gap/2-h-l, w/2], [-x_gap/2-h-l, -w/2], [-x_gap/2-h, -w/2], [-x_gap/2, 0]])
x_points = [-x_gap/2, -x_gap/2-h, -x_gap/2-h-l, -x_gap/2-h-l, -x_gap/2-h, -x_gap/2]
y_points = [ 0, w/2, w/2, -w/2, -w/2, 0]

geom_object = Geometry_emitter(emitter, round_corner=True)
x_coords = geom_object.emitter_points[:, 0]
y_coords = geom_object.emitter_points[:, 1]
tck2, u2 = spint.splprep([x_coords, y_coords], s=0)
unew = np.arange(0, 1.005, 0.005)
out2 = spint.splev(unew, tck2)
#plt.plot(out2[0], out2[1], '.-', out2[0][0], out2[1][0], '*')
out2[0][-1] = out2[0][0]
out2[1][-1] = out2[1][0]
emitter_coords = [Point(out2[0][i], out2[1][i]) for i in range(len(out2[0]))]
minx = -x_gap/2-h-l-10
maxx = x_gap 
miny = -w/2-10
maxy = w/2+10
emitter = Polygon(emitter_coords)
domain = Rectangle(Point(minx, miny), Point(maxx, maxy)) - emitter
mesh = generate_mesh(domain, 20)
print("Number of mesh vertices: ", mesh.num_vertices())
tree = mesh.bounding_box_tree()
function_space =  FunctionSpace(mesh, 'CG', 1)
print("function space dim:", function_space.dim())
trial_fxn = TrialFunction(function_space )
test_fxn = TestFunction(function_space )
tol = 1e-2

def out_boundary(x, on_boundary):
    #return on_boundary and (abs(x[0] - minx) <= tol or abs(x[0] - maxx) <= tol or abs(x[1] - miny) <= tol or abs(x[1] - maxy) <= tol) 
    return abs(x[0] - maxx) <= tol 

def on_emitter(x, on_boundary):
    return on_boundary and intersect_surface(np.transpose(out2), x, tol)

meshfile = File('poisson/mesh.pvd')
meshfile << mesh

a = inner(grad(trial_fxn), grad(test_fxn))*dx
L = Constant(0.0)*test_fxn*dx

out = DirichletBC(function_space, Constant(0), out_boundary)
em = DirichletBC(function_space, Constant(15), on_emitter)

bcs = [out, em]
parameters['reorder_dofs_serial'] = False
A, b = assemble_system(a, L, bcs)
solver = KrylovSolver('cg', 'ilu')
u = Function(function_space)
U = u.vector()
solver.solve(A, U, b)


bmesh = BoundaryMesh(mesh, "exterior", True)

plt.plot(bmesh.coordinates()[:,0], bmesh.coordinates()[:,1])
plt.scatter(bmesh.coordinates()[:,0], bmesh.coordinates()[:,1], c=range(len(bmesh.coordinates()[:,1])), cmap='gray')
plot(mesh)
#plt.plot(bmesh.coordinates()[96:,0],bmesh.coordinates()[96:,1], '-.', bmesh_pert.coordinates()[96:,0],bmesh_pert.coordinates()[96:,1], '.r')
plt.show()


for x in bmesh.coordinates()[12:12+len(out2[0]), :]:
    cell_index = tree.compute_first_entity_collision(Point(*x))
    cell_global_dofs = function_space.dofmap().cell_dofs(cell_index)
    cell = Cell(mesh, cell_index)
    values = np.array([0,])
    print(function_space.element().evaluate_basis_derivatives_all(1,x,cell.get_vertex_coordinates(),cell.orientation()))
    print(function_space.element().evaluate_basis_derivatives_all(1,x,cell.get_vertex_coordinates(),cell.orientation()).shape)

print(vertex_to_dof_map(function_space), len(vertex_to_dof_map(function_space)))
print(dof_to_vertex_map(function_space), len(dof_to_vertex_map(function_space)))





#print("Figure out what U is", abs(A.array()*U - b.get_local()))
print("shape of u:", U.get_local().shape)
print("shape of b:", b.get_local().shape)
print("shape of A:", A.array().shape)
print("number of mesh nodes:", mesh.num_vertices())
print("number of mesh cells:", mesh.num_cells())
#print("Content partialg_p:", partialg_p)

