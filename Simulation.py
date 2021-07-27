from __future__ import print_function
from Geometry import Geometry
from fenics import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np

class Simulation(object):

    def __init__(self, emitter, collector, gate, Ve, Vc, Vg, work_fxn):
        self.geometry = Geometry(emitter, collector, gate, round_corner=True)
        self.Ve = Ve
        self.Vc = Vc
        self.Vg = Vg
        self.work_fxn = work_fxn
        self.minx = min([np.min(self.geometry.emitter_points[:, 0]), np.min(self.geometry.collector_points[:, 0]), np.min(self.geometry.gate_points[:, 0])]) - 10
        self.miny = min([np.min(self.geometry.emitter_points[:, 1]), np.min(self.geometry.collector_points[:, 1]), np.min(self.geometry.gate_points[:, 1])]) - 10
        self.maxx = max([np.max(self.geometry.emitter_points[:, 0]), np.max(self.geometry.collector_points[:, 0]), np.max(self.geometry.gate_points[:, 0])]) + 10
        self.maxy = max([np.max(self.geometry.emitter_points[:, 1]), np.max(self.geometry.collector_points[:, 1]), np.max(self.geometry.gate_points[:, 1])]) +10
        emitter_coords = [Point(self.geometry.emitter_points[i, 0], self.geometry.emitter_points[i, 1]) for i in range(len(self.geometry.emitter_points))]
        collector_coords = [Point(self.geometry.collector_points[i, 0], self.geometry.collector_points[i, 1]) for i in range(len(self.geometry.collector_points))]
        gate_coords = [Point(self.geometry.gate_points[i, 0], self.geometry.gate_points[i, 1]) for i in range(len(self.geometry.gate_points))]
        emitter = Polygon(emitter_coords)
        collector = Polygon(collector_coords)
        gate = Polygon(gate_coords)
        self.domain = Rectangle(Point(self.minx, self.miny), Point(self.maxx, self.maxy)) - emitter - collector - gate
        self.mesh = generate_mesh(self.domain, 20)
        self.function_space =  FunctionSpace(self.mesh, 'CG', 1)
        self.trial_fxn = TrialFunction(self.function_space )
        self.test_fxn = TestFunction(self.function_space )
        self.electric_solution = None
        self.collector_current = np.zeros((len(Vc), len(Vg)))
        self.gate_current = np.zeros((len(Vc), len(Vg)))
        self.num_iter = 0

    def run_electric_sim(self, emitter_voltage, collector_voltage, gate_voltage):
        a = inner(grad(self.trial_fxn), grad(self.test_fxn))*dx
        L = Constant(0.0)*self.test_fxn*dx
        out = DirichletBC(self.function_space, Constant(0), out_boundary)
        emitter = DirichletBC(self.function_space, Constant(emitter_voltage), on_emitter)
        collector = DirichletBC(self.function_space, Constant(collector_voltage), on_collector)
        gate = DirichletBC(self.function_space, Constant(gate_voltage), self.on_gate)
        bcs = [out, emitter, collector, gate]
        A, b = assemble_system(a, L, bcs)
        solver = KrylovSolver('cg', 'ilu')
        if self.electric_solution is None:
            INITIAL_GUESS_ZERO = True
            u = Function(self.function_space)
            U = u.vector()
        else:
            INITIAL_GUESS_ZERO = False
            U = self.electric_solution
        solver.solve(A, U, b)
        self.electric_solution = U



    def run_particle_trajectory(self):
        pass

    def save_geometry(self, fname):
        meshfile = File('{}.pvd'.format(fname))
        meshfile << mesh

    def save_solutions(self, fname):
        plt.figure()
        plot(self.electric_solution)
        plt.savefig("{}.png".format(fname))
        plt.close()

    def out_boundary(self, x, on_boundary):
        return on_boundary and (abs(x[0] - self.minx) <= self.tol or abs(x[0] - self.maxx) <= self.tol or abs(x[1] - self.miny) <= self.tol or abs(x[1] - self.maxy) <= self.tol) 

    def on_emitter(self, x, on_boundary):
        return on_boundary and self.intersect_surface(self.geometry.emitter_points, x)

    def on_collector(self, x, on_boundary):
        return on_boundary and self.intersect_surface(self.geometry.collector_points, x)

    def on_gate(self, x, on_boundary):
        return on_boundary and self.intersect_surface(self.geometry.gate_points, x)

    def intersect_surface(self, outline_points, point_star):
        x_star, y_star = point_star
        is_on = False
        for i in range(len(outline_points[:, 0]) -1):
            if abs(outline_points[i, 0] - outline_points[i+1, 0]) < self.tol:
                if abs(x_star - outline_points[i, 0]) < self.tol and ((y_star <= outline_points[i,1] and y_star >= outline_points[i+1,1]) or (y_star >= outline_points[i,1] and y_star <= outline_points[i+1,1])):
                    is_on = True
            elif abs(outline_points[i, 1] - outline_points[i+1, 1]) < self.tol:
                if abs(y_star - outline_points[i, 1]) < self.tol and ((x_star <= outline_points[i,0] and x_star >= outline_points[i+1,0]) or (x_star >= outline_points[i,0] and x_star <= outline_points[i+1,0])):
                    is_on = True
            elif (x_star <= outline_points[i,0] and x_star >= outline_points[i+1,0]) or (x_star >= outline_points[i,0] and x_star <= outline_points[i+1,0]):
                m = (outline_points[i,1] - outline_points[i+1,1])/(outline_points[i,0] - outline_points[i+1,0])
                b = outline_points[i,1] - m*outline_points[i,0]
                if abs(y_star - m*x_star - b) < self.tol:
                    is_on = True
        return is_on
    
    def cost_function(self, alpha, beta):
        Vce = (self.Vc - self.Ve).T
        # Collector_current is size (len(Vc), len(Vg))
        current_slope = self.collector_current[1:len(collector_current), :] - collector_current[0:len(collector_current)-1, :]
        for i in range(len(collector_current[0, :])):
            current_slope[:, i] = np.divide(np.multiplyy(Vce[0:len(Vce)-1], current_slope[:, i]), (Vce[1:len(Vce)] - Vce[0:len(Vce-1)]))
        
        for j in range(len(current_slope[:, 0])):
            current_slope[j, :] = np.multiply(current_slope[j, :], Vg)

        cost = np.sum(current_slope) + alpha*np.sum(gate_current) + beta*1/np.sum(collector_current)
        return cost

minx = -150
miny = -150
maxx = 150
maxy = 150
tol = 1e-2
def out_boundary(x, on_boundary):
    return on_boundary and (abs(x[0] - minx) <= tol or abs(x[0] - maxx) <= tol or abs(x[1] - miny) <= tol or abs(x[1] - maxy) <= tol) 

def on_emitter(x, on_boundary):
    geometry = Geometry(emitter, collector, gate, round_corner=True)
    return on_boundary and intersect_surface(geometry.emitter_points, x)

def on_collector(self, x, on_boundary):
    geometry = Geometry(emitter, collector, gate, round_corner=True)
    return on_boundary and self.intersect_surface(geometry.collector_points, x)