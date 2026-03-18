
import ipctk
import meshio
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
import networkx as nx
import scipy.spatial.distance as dst
import matplotlib.pyplot as plt
import math
import random
import sys

pio.renderers.default = "browser"

def plot_wireframe(points, faces):
    points = np.asarray(points)
    faces = np.asarray(faces)
    if points.shape[1] == 2:
        points = np.column_stack([points, np.zeros(len(points))])
    edges = set()
    for tri in faces:
        pairs = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]
        for a, b in pairs:
            if a > b:
                a, b = b, a
            edges.add((a, b))
    edges = np.array(list(edges))
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    xs, ys, zs = [], [], []
    for a, b in edges:
        xs += [x[a], x[b], None]
        ys += [y[a], y[b], None]
        zs += [z[a], z[b], None]
    fig = go.Figure(
        go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="lines",
        )
    )
    fig.update_layout(scene=dict(aspectmode="data"))
    fig.show()

def plot_quad_wireframe(points, edges_path, faces, thickness=8):
    edges = []
    with open(edges_path) as file:
        for line in file.readlines():
            e = list(map(int, line.split(" ")))
            edges.append(e)
    xs = []
    ys = []
    zs = []
    figs = []
    for e in edges:
        if len(e) == 2:
            a, b = e
            xs += [points[a-1][0], points[b-1][0], None]
            ys += [points[a-1][1], points[b-1][1], None]
            zs += [points[a-1][2], points[b-1][2], None]
        elif len(e) == 4:
            _xs = []
            _ys = []
            _zs = []
            a, b, c, d = e
            _xs = [points[a-1][0], points[b-1][0], None,points[b-1][0], points[c-1][0], None,points[c-1][0], points[d-1][0], None]
            _ys = [points[a-1][1], points[b-1][1], None,points[b-1][1], points[c-1][1], None,points[c-1][1], points[d-1][1], None]
            _zs = [points[a-1][2], points[b-1][2], None,points[b-1][2], points[c-1][2], None,points[c-1][2], points[d-1][2], None]
            figs.append(go.Scatter3d(x=_xs,y=_ys,z=_zs, mode="lines", line=dict(width=thickness,color=f'rgb({random.randint(0,255)}, {random.randint(0,255)}, {random.randint(0,255)})')))
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
    i = faces[:,0]
    j = faces[:,1]
    k = faces[:,2]
    fig = go.Figure(
        data = [go.Mesh3d(x=x,y=y,z=z,i=i,j=j,k=k,color="gray",opacity=0.2), go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="lines",
            line=dict(
                color="red",
                width=thickness
            )
        )] + figs
    )
    fig.update_layout(scene=dict(aspectmode="data"))
    fig.show()

def triangulated_meshio_to_graph(mesh, triangle_key="triangle", weight_attr="weight"):
    pts = np.asarray(mesh.points, dtype=float)
    tris = mesh.cells_dict[triangle_key]
    tris = np.asarray(tris, dtype=np.int64)
    G = nx.Graph()
    G.add_nodes_from(range(len(pts)))
    edges = np.vstack([
        tris[:, [0, 1]],
        tris[:, [1, 2]],
        tris[:, [2, 0]],
    ])
    edges.sort(axis=1)
    edges = np.unique(edges, axis=0)
    u = edges[:, 0]
    v = edges[:, 1]
    w = np.linalg.norm(pts[u] - pts[v], axis=1)
    for a, b, d in zip(u, v, w):
        G.add_edge(int(a), int(b))
    return G

def all_pairs_hop_matrix(G, dtype=np.int32):
    n = G.number_of_nodes()
    D = np.full((n, n), -1, dtype=dtype)   # -1 = unreachable
    np.fill_diagonal(D, 0)
    for s, dist in nx.all_pairs_shortest_path_length(G):
        for t, d in dist.items():
            D[s, t] = d
    return D

def spring_energy_gradient(V, M, true_dist, eps=1e-12):
    V = np.asarray(V, dtype=float)
    M = np.asarray(M, dtype=float)
    iu, ju = np.triu_indices(M.shape[0], k=1)
    L0 = M[iu, ju]
    mask = L0 > 0
    iu = iu[mask]
    ju = ju[mask]
    L0 = L0[mask]
    L0_truedist = true_dist[iu,ju]
    if L0.size == 0:
        return np.zeros_like(V)
    d = V[iu] - V[ju]
    dist = np.linalg.norm(d, axis=1)
    inv_dist = np.where(dist > eps, 1.0 / dist, 0.0)
    scale = (dist - L0) * inv_dist * (1 / (2 ** (L0_truedist / 2.0)))
    g = d * scale[:, None]
    grad = np.zeros_like(V)
    np.add.at(grad, iu, g)
    np.add.at(grad, ju, -g)
    return grad

def spring_energy(V, M,true_dist):
    dists = dst.squareform(dst.pdist(V))
    y = 1/(2 ** (true_dist / 2.0)) * ((dists - M) ** 2)
    return np.mean(y)

def average_dist(V):
    dists = dst.squareform(dst.pdist(V))
    return np.mean(dists)

def average_graph_dist(M):
    return np.mean(M)

def desired_dist(V, M):
    return average_dist(V)/average_graph_dist(M)

def read_distance_matrix(path):
    distances = []
    with open(path) as file:
        for line in file:
            distances.append(int(line))
        n = len(distances)
        v = int(math.sqrt(n))
        distances = np.array(distances)
        distances = np.reshape(distances, (v, v))
        distances = np.asarray(distances, dtype=np.float32)
        return distances

distances = read_distance_matrix(sys.argv[2])
true_dist = distances.copy()
mesh = meshio.read(sys.argv[1])
rest_positions = np.asarray(mesh.points, dtype=np.float64)
faces = mesh.cells_dict["triangle"]
edges = ipctk.edges(faces)
collision_mesh = ipctk.CollisionMesh(rest_positions, edges, faces,)
distances = desired_dist(rest_positions, distances) * distances
vertices = collision_mesh.rest_positions.copy()

def iteration(vertices, collision_mesh, dhat, stiffness, initial_step_size, distances, spring_energies, spring_weight=1):
    collisions = ipctk.NormalCollisions()
    collisions.build(collision_mesh, vertices, dhat)
    B = ipctk.BarrierPotential(dhat, stiffness)
    barrier_potential = B(collisions, collision_mesh, vertices)
    barrier_grad = B.gradient(collisions, collision_mesh, vertices)
    barrier_grad = np.reshape(barrier_grad, (int(barrier_grad.shape[0] / 3), 3))
    spring_grad = spring_energy_gradient(vertices, distances, true_dist)
    #electric_grad = coulomb_repulsion_grad(vertices)
    grad = barrier_grad + spring_grad #+ electric_grad
    new_vertices = vertices.copy() - initial_step_size * grad
    max_step_size = ipctk.compute_collision_free_stepsize(collision_mesh, vertices, new_vertices)*0.8
    collision_free_vertices = (new_vertices - vertices) * max_step_size + vertices
    assert(ipctk.is_step_collision_free(collision_mesh, vertices, collision_free_vertices))
    return collision_free_vertices
    #return new_vertices

v = vertices
spring_energies = []
for i in range(int(sys.argv[3])):
    if i % 50 == 0:
        print("=========== " + str(i) + " ===========")
    v = iteration(v, collision_mesh, 1e-3, 1, 1e-3, distances, spring_energies, spring_weight=2)
mesh.points = v
meshio.write(sys.argv[4], mesh)
