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
    
def get_per_points(x, y, epsilon):
    """ finds better points to calculate fields at """
    x_new = []
    y_new = []
    for i in range(len(x)-1):
        x_vec = x[i+1] - x[i]
        y_vec = y[i+1] - y[i]
        norm_factor = epsilon/np.sqrt(x_vec**2 + y_vec**2) 
        x_new.append( x[i] + norm_factor*(y[i+1] - y[i]))
        y_new.append(y[i] - norm_factor*(x[i+1] + x[i]))
    return x_new, y_new

def fowler_nordheim_emission(E_strength, area, phi):
    afn = 1.5414e-6 # eV v^-2
    bfn = 6.8309 # eV^(-3/2)V/nm
    nu = 1
    J = afn*np.power(E_strength, 2)/phi*np.exp(-nu*bfn*np.power(phi, 3/2)/E_strength)
    return J*area

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

emitter_coords = [Point(out2[0][i], out2[1][i]) for i in range(len(out2[0]))]
minx = -x_gap/2-h-l-10
maxx = x_gap 
miny = -w/2-10
maxy = w/2+10
emitter = Polygon(emitter_coords)
domain = Rectangle(Point(minx, miny), Point(maxx, maxy)) - emitter
mesh = generate_mesh(domain, 20)
print("Number of mesh vertices: ", mesh.num_vertices())

function_space =  FunctionSpace(mesh, 'CG', 1)
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

A, b = assemble_system(a, L, bcs)
solver = KrylovSolver('cg', 'ilu')
u = Function(function_space)
U = u.vector()
solver.solve(A, U, b)

# get del\delp 
epsilon = 0.001
x_pert, y_pert = get_per_points(x_coords, y_coords, epsilon)
partialg_p = np.zeros((len(b), len(x_pert)))
for j in range(len(x_pert)):
    x = x_coords
    x[j] = x_pert[j]
    y = y_coords
    y[j] = y_pert[j]
    tck, u = spint.splprep([x, y], s=0)
    unew = np.arange(0, 1.005, 0.005)
    out_pert = spint.splev(unew, tck)
    emitter_coords_pert = [Point(out_pert[0][i], out_pert[1][i]) for i in range(len(out_pert[0]))]
    emitter_pert = Polygon(emitter_coords_pert)
    domain_pert = Rectangle(Point(minx, miny), Point(maxx, maxy)) - emitter_pert
    mesh_pert = generate_mesh(domain_pert, 20)
    bmesh = BoundaryMesh(mesh, "exterior", True)
    bmesh_pert = BoundaryMesh(mesh_pert, "exterior", True)
    print("Length comparison. \nOld boundary: {}\nNew bounary: {}".format(len(bmesh), len(bmesh_pert)))

    # plt.figure()
    # plt.plot(out2[0], out2[1], out_pert[0], out_pert[1]) #x_points, y_points, 'b'
    # plt.show()
    # print("Number of pert mesh vertices: ", mesh_pert.num_vertices())
    # function_space_pert =  FunctionSpace(mesh_pert, 'CG', 1)
    # trial_fxn_pert = TrialFunction(function_space_pert)
    # test_fxn_pert = TestFunction(function_space_pert)

    # def on_emitter_pert(x, on_boundary):
    #     return on_boundary and intersect_surface(np.transpose(out_pert), x, tol)

    # a_pert = inner(grad(trial_fxn_pert), grad(test_fxn_pert))*dx
    # L_pert = Constant(0.0)*test_fxn_pert*dx

    # out = DirichletBC(function_space_pert, Constant(0), out_boundary)
    # em_pert = DirichletBC(function_space_pert, Constant(15), on_emitter_pert)

    # bcs_pert = [out, em_pert]
    # A_pert, b_pert = assemble_system(a_pert, L_pert, bcs_pert)
    # partialg_p[:,j] = np.matmul((A_pert.array() - A.array())/epsilon, U.get_local()) - (b_pert.get_local() - b.get_local())/epsilon


#print("Figure out what U is", abs(A.array()*U - b.get_local()))
print("shape of u:", U.get_local().shape)
print("shape of b:", b.get_local().shape)
print("shape of A:", A.array().shape)
print("shape of Au:", np.matmul(A.array(), U.get_local()).shape)
print("shape of partialg_p:", partialg_p.shape)

