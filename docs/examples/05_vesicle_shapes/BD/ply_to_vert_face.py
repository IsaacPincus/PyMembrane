import open3d as o3d
import numpy as np

mesh = o3d.io.read_triangle_mesh('sphere_R1_l3.ply')

vertices = np.asarray(mesh.vertices)
triangles = np.asarray(mesh.triangles)

with open('vertices_sphere_N1024.inp', 'w') as file:
    for index, vertex in enumerate(vertices):
        file.write(str(index) + ' ' + ' '.join(map(str, vertex)) + ' 1 \n')


with open('faces_sphere_N1024.inp', 'w') as file:
    for index, face in enumerate(triangles):
        file.write(str(index) + ' ' + ' '.join(map(str, face)) + ' 1 1 \n')


mesh = o3d.io.read_triangle_mesh('ellipse.ply')

vertices = np.asarray(mesh.vertices)
triangles = np.asarray(mesh.triangles)

with open('vertices_ellipse.inp', 'w') as file:
    for index, vertex in enumerate(vertices):
        file.write(str(index) + ' ' + ' '.join(map(str, vertex)) + ' 1 \n')


with open('faces_ellipse.inp', 'w') as file:
    for index, face in enumerate(triangles):
        file.write(str(index) + ' ' + ' '.join(map(str, face)) + ' 1 1 \n')