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
x = geom_object.emitter_points[:, 0]
y = geom_object.emitter_points[:, 1]
x[3]+=2
y[3]+=2
x[4]+=7
y[4]+=5
x[5]+=2
y[5]+=2
tck, u = spint.splprep([x_points, y_points], s=0)
tck2, u2 = spint.splprep([x, y], s=0)
unew = np.arange(0, 1.005, 0.005)
out = spint.splev(unew, tck)
out2 = spint.splev(unew, tck2)
plt.figure()
plt.plot(x_points, y_points, 'x', out[0], out[1], out2[0], out2[1]) #x_points, y_points, 'b'
plt.legend(['Points', 'Cubic Spline1', 'Cubic spline2', 'Linear'])
plt.title('Spline of parametrically-defined curve')
plt.show()

emitter_coords = [Point(out2[0][i], out2[1][i]) for i in range(len(out2[0]))]
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
    #return on_boundary and (abs(x[0] - minx) <= tol or abs(x[0] - maxx) <= tol or abs(x[1] - miny) <= tol or abs(x[1] - maxy) <= tol) 
    return abs(x[0] - maxx) <= tol 

def on_emitter(x, on_boundary):
    return on_boundary and intersect_surface(np.transpose(out2), x, tol)

meshfile = File('poisson/mesh.pvd')
meshfile << mesh

a = inner(grad(trial_fxn), grad(test_fxn))*dx
L = Constant(0.0)*test_fxn*dx

out = DirichletBC(function_space, Constant(0), out_boundary)
em = DirichletBC(function_space, Constant(10), on_emitter)

bcs = [out, em]

A, b = assemble_system(a, L, bcs)
solver = KrylovSolver('cg', 'ilu')
u = Function(function_space)
U = u.vector()
solver.solve(A, U, b)
f = File("poisson/solution.pvd")
f << u

V_g = VectorFunctionSpace(mesh, 'CG', 1)
v = TestFunction(V_g)
w = TrialFunction(V_g)
a = inner(w, v)*dx
L = inner(grad(u), v)*dx
grad_u = Function(V_g)
solve(a == L, grad_u)
plt.figure()
plot(u)
#plot(grad(u))
x_coord = [-104, -45, -32, -22, -18, -11, -22, -35, -48, 0, 10, -25, 3, -2]
y_coord = [54, 43, 30, 28, 10, 0, -15, -27, -40, 0, 0, 22, 38, -44]
U1 = []
V1 = []
U2 = []
V2 = []
for i in range(len(x_coord)):
    grad_at_x = grad_u(Point(x_coord[i], y_coord[i]))
    U1.append(grad_at_x[0])
    V1.append(grad_at_x[1])
    U2.append(x_coord[i]+ grad_at_x[0])
    V2.append(y_coord[i] + grad_at_x[1])
plt.quiver(x_coord, y_coord, U1, V1)
plt.show()
plt.savefig("poisson/potential_and_field1.png")
plt.close()

plt.figure()
plot(u)
plt.quiver(x_coord, y_coord, U2, V2)
plt.show()
plt.savefig("poisson/potential_and_field2.png")
plt.close()
# plot(grad_u)
# plt.show()
# plt.savefig("poisson/grad2.png")
print("potential at 0,0:", u(Point(0,0)))
print("E-field at 0,0:", grad_u(Point(0,0)))
print("E-field at -3,0:", grad_u(Point(-3,0)))
print("E-field at -10,0:", grad_u(Point(-10,0)))
print("E-field at -34,-27:", grad_u(Point(-34,-27)))

f2 = File("poisson/grad.pvd")
f2 << grad(u)
