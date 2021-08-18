from scipy import interpolate
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
x = geom_object.emitter_points[:, 0]
y = geom_object.emitter_points[:, 1]
x[25]+=1
y[25]-=5
tck, u = interpolate.splprep([x_points, y_points], s=0)
tck2, u2 = interpolate.splprep([x, y], s=0)
unew = np.arange(0, 1.005, 0.005)
out = interpolate.splev(unew, tck)
out2 = interpolate.splev(unew, tck2)
plt.figure()
plt.plot(x_points, y_points, 'x', out[0], out[1], out2[0], out2[1]) #x_points, y_points, 'b'
plt.legend(['Points', 'Cubic Spline1', 'Cubic spline2', 'Linear'])
plt.title('Spline of parametrically-defined curve')
plt.show()

mitter_coords = [Point(out2[0][i], out2[1][i]) for i in range(len(out2[0]))]
minx = -x_gap/2-h-l-10
maxx = x_gap 
miny = -w/2-10
maxy = w/2+10
emitter = Polygon(emitter_coords)
domain = Rectangle(Point(minx, miny), Point(maxx, maxy)) - emitter
mesh = generate_mesh(domain, 20)

function_space =  FunctionSpace(mesh, 'CG', 1)
trial_fxn = TrialFunction(function_space )
test_fxn = TestFunction(function_space )
tol = 1e-2

def out_boundary(x, on_boundary):
    return on_boundary and (abs(x[0] - minx) <= tol or abs(x[0] - maxx) <= tol or abs(x[1] - miny) <= tol or abs(x[1] - maxy) <= tol) 

def on_emitter(x, on_boundary):
    return on_boundary and intersect_surface(np.transpose(out2), x, tol)

meshfile = File('poisson/mesh.pvd')
meshfile << mesh

u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = Constant(0.0)*v*dx

out = DirichletBC(V, Constant(0), out_boundary)
em = DirichletBC(V, Constant(10), bemitter)

bcs = [out, em]

A, b = assemble_system(a, L, bcs)
solver = KrylovSolver('cg', 'ilu')
u = Function(V)
U = u.vector()
solver.solve(A, U, b)
f = File("poisson/solution.pvd")
f << u
plt.figure()
plot(u)
plt.savefig("poisson/solution.png")
plt.close()
