
import numpy as np 
import matplotlib.pyplot as plt

class Geometry(object):
    segment_length = 5
    tol = 0.01

    def __init__(self, emitter_vertices, collector_vertices, gate_vertices, round_corner=True):
        # input geometry
        self.emitter_vertices = emitter_vertices
        self.collector_vertices = collector_vertices
        self.gate_vertices = gate_vertices
        
        self.emitter_points = self.create_points(emitter_vertices, round_corner)
        self.collector_points = self.create_points(collector_vertices, round_corner)
        self.gate_points = self.create_points(gate_vertices, round_corner)
    
    def create_points(self, vertices, round_corners):
        points = []
        indecies_of_corners = []
        corner_ind = 0
        for i in range(len(vertices[:, 0])-1):
            seg_points = self.get_segments_between_points(vertices[i, :], vertices[i+1, :], self.segment_length)
            points.append(seg_points)
            indecies_of_corners.append(corner_ind)
            corner_ind += len(seg_points)
        points = np.concatenate(points)
        points = np.append(points, [vertices[-1, :]], axis=0)
        indecies_of_corners.append(corner_ind)
        if not round_corners:
            return points
        else:
            #round corners
            for j in range(1, len(indecies_of_corners)-1):
                ind1 = indecies_of_corners[j] -3
                ind2 = indecies_of_corners[j] +3
                pts, xs, ys = self.round_corners(points, ind1, ind2, self.segment_length)
                points[ind1+1:ind2, :] = pts
            pts, xs, ys = self.round_corners(points, len(points)-3, 2, 3)
            points_pr = points[2:len(points)-2,:]
            points = np.concatenate((points_pr, pts, [points[2, :]]))
            return points

                                                                                
    def get_segments_between_points(self, p1, p2, seglen):
        x1, y1 = p1
        x2, y2 = p2
        Nsegs = round(np.sqrt((x1-x2)**2 +(y1-y2)**2)/seglen).astype(int)
        if Nsegs < 2:
            points = np.array([[x1, y1], [x2, y2]])
        else:
            points = np.zeros((Nsegs, 2))
            points[0, :] = [x1, y1]
            for i in range(1, Nsegs):
                x_t = (1- i/Nsegs)*x1 + i/Nsegs*x2
                y_t = (1- i/Nsegs) *y1 + i/Nsegs*y2
                points[i, :] = [x_t, y_t]
        return points

    def round_corners(self, coords, ind1, ind2, N):
        x1, y1 = coords[ind1, :]
        x2, y2 = coords[ind2, :]
        if abs(y1 - coords[ind1+1, 1]) > self.tol and abs(x1 - coords[ind1+1, 0]) > self.tol:
            m1 = (y1 - coords[ind1+1, 1])/(x1 - coords[ind1+1, 0])
            m1_perp = -1/m1
            b1 = y1 - m1_perp*x1

        if abs(y2 - coords[ind2-1, 1]) > self.tol and abs(x2 - coords[ind2-1, 0]) > self.tol:
            m2 = (y2 - coords[ind2-1, 1])/(x2 - coords[ind2-1, 0])
            m2_perp = -1/m2
            b2 = y2 - m2_perp*x2

        if abs(y1 - coords[ind1+1, 1]) < self.tol and abs(x2 - coords[ind2-1, 0]) < self.tol:
            xc = x1
            yc = y2
            if y2 < y1: 
                phi_initial = np.pi/2
                phi_final = np.pi
            else:
                phi_initial = -np.pi/2
                phi_final = 0

        elif abs(x1 - coords[ind1+1, 0]) < self.tol and abs(y2 - coords[ind2-1, 1]) < self.tol:
            xc = x2
            yc = y1
            if y1 < y2:
                phi_initial = 0
                phi_final = np.pi/2
            else: 
                phi_initial = np.pi
                phi_final = 3*np.pi/2

        elif abs(y1 - coords[ind1+1, 1]) < self.tol:
            xc = x1
            yc = m2_perp*xc +b2
            if y1 > yc:
                phi_initial = np.pi/2
                phi_final = np.arctan((y2 - yc)/(x2 - xc)) +np.pi
            else:
                phi_initial = -np.pi/2
                phi_final = np.arctan((y2 - yc)/(x2 - xc))

        elif abs(y2 -coords[ind2-1, 1]) < self.tol:
            xc = x2
            yc = m1_perp*xc +b1
            if y2 > yc:
                phi_initial = np.arctan((y1 - yc)/(x1 - xc))
                phi_final = np.pi/2
            else:
                phi_initial = np.arctan((y1 - yc)/(x1 - xc)) +np.pi
                phi_final = 3*np.pi/2

        elif abs(x1 - coords[ind1+1, 0]) < self.tol:
            yc = y1
            xc = (yc - b2)/m2_perp
            if x1 > xc:
                phi_initial = 0
            else:
                phi_initial = np.pi
            if x2 > xc:
                phi_final = np.arctan((y2 - yc)/(x2 - xc))
            else:
                phi_final = np.arctan((y2 - yc)/(x2 - xc)) + np.pi
    
        elif abs(x2 - coords[ind2+1, 0]) < self.tol:
            yc = y2
            xc = (yc - b1)/m1_perp
            if x1 > xc: 
                phi_initial = np.arctan((y1 - yc)/(x1 - xc))
            else:
                phi_initial = np.arctan((y1 - yc)/(x1 - xc))
            if x2 > xc:
                phi_final = 0
            else:
                phi_final = np.pi
        
        else:
            xc = (b1 - b2)/(m2_perp - m1_perp)
            yc = m1_perp*xc +b1
            if (x1 - xc) < 0:
                phi_initial = np.arctan((x1 - xc)/(y1 - yc)) + np.pi
            else:
                phi_initial = np.arctan((y1 - yc)/(x1 - xc))
            if (x2 - xc) < 0:
                phi_final = np.arctan((y2 - yc)/(x2 - xc)) + np.pi
            else:
                phi_final = np.arctan((y2 - yc)/(x2 - xc))

        R = np.sqrt((xc-x1)**2 + (yc-y1)**2)
        new_coords = np.zeros((N, 2))

        if phi_final < phi_initial:
            phi_final = phi_final +2*np.pi

        for i in range(1, N+1):
            phi_i = phi_initial + i*(phi_final - phi_initial)/(N+1)
            xi = xc + R*np.cos(phi_i)
            yi = yc + R*np.sin(phi_i)
            new_coords[i-1, :] = [xi, yi]

        phi = np.linspace(0, 2*np.pi)
        xs = xc + R*np.cos(phi)
        ys = yc + R*np.sin(phi)
        return new_coords, xs, ys

    def draw_geom(self):
        plt.plot(self.emitter_points[:, 0], self.emitter_points[:, 1], '.-r')
        plt.plot(self.collector_points[:, 0], self.collector_points[:, 1], '.-b')
        plt.plot(self.gate_points[:, 0], self.gate_points[:, 1], '.-g')
        plt.show()
    



