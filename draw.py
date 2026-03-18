import meshio
import plotly.graph_objects as go
import plotly.io as pio
import random
import sys
import numpy as np
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

mesh = meshio.read(sys.argv[1])
rest_positions = np.asarray(mesh.points, dtype=np.float64)
faces = mesh.cells_dict["triangle"]
plot_quad_wireframe(rest_positions, sys.argv[1] + ".edges", faces)
