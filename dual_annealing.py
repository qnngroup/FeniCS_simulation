# dual annealing global optimization for the ackley multimodal objective function
from scipy.optimize import dual_annealing
from numpy.random import rand
import scipy.interpolate as spint
import numpy as np
import matplotlib.pyplot as plt
from Geometry_emitter import Geometry_emitter
from fenics import *
from mshr import *


def intersect_surface(outline_points, point_star, tol):
    """ used for determaining points on the inner boundary for boundary conditions. 
    Takes in a test point and an array of points which define a surface. 
    The code linearly iterates through the outline array and checks whether the test point is on the line connecting point i with point i+1. 
    Returns a boolean value """
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
    

def find_ppoints(x, y):
    """ The code takes in an array of x values and an array of y values. Together they describe an outline of a shape. 
    The code uses the cross product to generate the outward pointing normal for every point on the outline. 
    The outward normal is used to find a new point that is in the simulation domian, just outside the shape. 
    This is used to find a point just outside the emitter to test the e field in because the gradient isn't defined on the boundary (apperantly)   """
    x_new = []
    y_new = []
    for i in range(len(x)-1):
        x_vec = x[i+1] - x[i]
        y_vec = y[i+1] - y[i]
        x_new.append( x[i] + y[i+1] - y[i])
        y_new.append(y[i] -x[i+1] + x[i])
    return x_new, y_new

def fowler_nordheim_emission(E_strength, area, phi):
    """ Takes in an electric fielsd magnitude, area, and work function. Determains the current due to fowler nordheim emission in that area
    """
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
x_coords[3]+=2
y_coords[3]+=2
x_coords[4]+=7
y_coords[4]+=5
x_coords[5]+=2
y_coords[5]+=2

minx = -x_gap/2-h-l-10
maxx = x_gap 
miny = -w/2-10
maxy = w/2+10

outline_points = np.concatenate((x_coords, y_coords))    
def objective(outline_points):
    """ Objective function solves the poisson equation and calculates the gradient to find the e field strength"""
    x_coords = outline_points[0:int(len(outline_points)/2)]
    y_coords = outline_points[int(len(outline_points)/2):len(outline_points)]
    tck2, u2 = spint.splprep([x_coords, y_coords], s=0)
    unew = np.arange(0, 1.005, 0.005)
    out2 = spint.splev(unew, tck2)
    plt.figure()
    plt.plot(x_coords, y_coords, 'x', out2[0], out2[1]) #x_points, y_points, 'b'
    plt.legend(['Points', 'Cubic Spline1', 'Cubic spline2', 'Linear'])
    plt.title('Spline of parametrically-defined curve')
    plt.show()

    emitter_coords = [Point(out2[0][i], out2[1][i]) for i in range(len(out2[0]))]
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
    em = DirichletBC(function_space, Constant(15), on_emitter)

    bcs = [out, em]

    A, b = assemble_system(a, L, bcs)
    solver = KrylovSolver('cg', 'ilu')
    u = Function(function_space)
    U = u.vector()
    solver.solve(A, U, b)
    f = File("poisson/solution.pvd")
    f << u

    #print("Figure out what U is", abs(A.array()*U - b.get_local()))
    # print("shape of u:", U.get_local().shape)
    # print("shape of b:", b.get_local().shape)
    # print("shape of A:", A.array().shape)
    # print("shape of Au:", np.matmul(A.array(), U.get_local()).shape)
    # plt.figure()
    # plt.hist(np.matmul(A.array(), U.get_local()) - b.get_local(), 20)
    # plt.show()
    #print("Au -b :", np.matmul(A.array(), U.get_local()) - b.get_local())
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
    U1 = []
    V1 = []
    magnitude = []
    area = []
    x_coords = out2[0]
    y_coords = out2[1]
    x,y = find_ppoints(x_coords, y_coords)
    for i in range(len(x)):
        grad_at_x = grad_u(Point(x[i], y[i]))
        U1.append(-grad_at_x[0])
        V1.append(-grad_at_x[1])
        magnitude.append(np.sqrt((grad_at_x[0])**2 + (grad_at_x[1])**2)*1e9)
        if i<len(x) - 1:
            seg1 = np.sqrt((x[i] - x[i-1])**2 + (y[i] -y[i-1])**2)
            seg2 =  np.sqrt((x[i+1] - x[i])**2 + (y[i+1] -y[i])**2)
            area.append(1e-9*(seg1+seg2)/2)
        else:
            seg1 = np.sqrt((x[i] - x[i-1])**2 + (y[i] -y[i-1])**2)
            seg2 =  np.sqrt((x[0] - x[i])**2 + (y[0] -y[i])**2)
            area.append(1e-9*(seg1+seg2)/2)

    plt.quiver(x, y, U1, V1)
    #plt.scatter(x,y)
    total_current = np.sum(fowler_nordheim_emission(magnitude, area, 5.3))*40e-9 * 1e6 # get current as micro amps
    print("Total current:", total_current)
    plt.text(-120, 5, 'Total current:\n {} \mu A'.format(total_current), fontsize=12)
    plt.show()
    plt.savefig("poisson/potential_and_field.png")
    plt.close()
    return 1/total_current




# define the bounds on the search
bounds = []
for x in outline_points[0:int(len(outline_points)/2)]:
    bounds.append([minx+5, maxx-5])
for y in outline_points[int(len(outline_points)/2):len(outline_points)]:
    if y>10:
        bounds.append([5, maxy])
    elif y<-10:
        bounds.append([miny, -5])
    else:
        bounds.append([-10, 10])

# perform the dual annealing search
result = dual_annealing(objective, bounds, visit=1.3, initial_temp=5000, x0=outline_points)
# summarize the result
print('Status : %s' % result['message'])
print('Total Evaluations: %d' % result['nfev'])
# evaluate solution
solution = result['x']
evaluation = objective(solution)
print('Solution: f(%s) = %.5f' % (solution, evaluation))