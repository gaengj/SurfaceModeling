#!usr/bin/env python

"""
Developability Algorithm
"""

from math import *
import pymesh
import numpy as np

def g(mesh, idx_vertice, border_vertices_list):
    """
    Compute the vertex developability detect function g(q) at each vertex q on the given mesh patches
    idx_vertice is the index of the vertex q
    """
    mesh.enable_connectivity()
    idx_vertice = int(idx_vertice)
    if idx_vertice not in border_vertices_list:
        list_idx_adjacent_faces = mesh.get_vertex_adjacent_faces(idx_vertice)
        somme_angles = 0
        for idx_face in list_idx_adjacent_faces:
            face = mesh.faces[idx_face]
            if face[0] == idx_vertice:
                edge1 = mesh.vertices[face[1]] - mesh.vertices[face[0]]
                edge2 = mesh.vertices[face[2]] - mesh.vertices[face[0]]
            elif face[1] == idx_vertice:
                edge1 = mesh.vertices[face[0]] - mesh.vertices[face[1]]
                edge2 = mesh.vertices[face[2]] - mesh.vertices[face[1]]
            else:
                edge1 = mesh.vertices[face[0]] - mesh.vertices[face[2]]
                edge2 = mesh.vertices[face[1]] - mesh.vertices[face[2]]
            dot_product = np.dot(edge1, edge2)/(np.linalg.norm(edge1)*np.linalg.norm(edge2))
            if (dot_product > 1):
                dot = 1.0
            elif (dot_product < -1):
                dot = -1.0
            else:
                dot = dot_product
            angle = acos(dot)
            somme_angles += angle
        return(2.0*pi - somme_angles)
    else:
        return(0)


def T(mesh, mesh_normal, border_vertices_list, delta, idx_vertice):
    """
    Compute the function T(delta)
    """
    mesh.enable_connectivity()
    somme = 0
    idx_vertice = int(idx_vertice)
    list_idx_adjacent_vertices = mesh.get_vertex_adjacent_vertices(idx_vertice)
    for qi in list_idx_adjacent_vertices:
        somme += g(mesh, qi, border_vertices_list) ** 2

    #We have to copy the vertices to modify it
    vertices_copy = mesh.vertices.copy()
    delta_normal = [delta * mesh_normal[3*idx_vertice], delta * mesh_normal[3*idx_vertice + 1], delta * mesh_normal[3*idx_vertice+2]]
    vertices_copy[idx_vertice] = mesh.vertices[idx_vertice] + delta_normal
    new_mesh = pymesh.form_mesh(vertices_copy, mesh.faces)
    somme += g(new_mesh, idx_vertice, border_vertices_list) **2
    return(somme)

def T_point(mesh, mesh_normal, border_vertices_list, delta, idx_vertice, h):
    """
    Compute the derivative of T with the formula : (T(delta+h) - T(delta))/h
    """
    return (T(mesh, mesh_normal, border_vertices_list, delta + h, idx_vertice)-T(mesh, mesh_normal, border_vertices_list, delta, idx_vertice))/h




def main():
    mesh = pymesh.load_mesh("/models/models/front_dress_anim/meshes/mesh_00040.off")
    # mesh = pymesh.load_mesh("/models/models/test_plane.obj")
    mesh.enable_connectivity()

    N_max = 1000
    epsilon = 0.001

    # 0 Find the edge border of the mesh

    edge_list = []
    border_vertices_list=[]

    for [x,y,z] in mesh.faces:
        edge_list.append((x,y))
        edge_list.append((y,x))
        edge_list.append((y,z))
        edge_list.append((z,y))
        edge_list.append((z,x))
        edge_list.append((x,z))

    compte = {}.fromkeys(set(edge_list),0)
    for value in edge_list:
        compte[value] += 1

    for key in compte.keys():
        if compte.get(key) == 1:
            border_vertices_list.append(key[0])
            border_vertices_list.append(key[1])

    border_vertices_list = set(border_vertices_list)


    # 1 Compute the vertex developability detect function g(q)
    # at each vertex q on the given mesh patches

    list_g = np.array(len(mesh.vertices)*[0.0])

    for qi in range (len(mesh.vertices)):
        list_g[qi] = g(mesh, qi, border_vertices_list)

    # 2 Compute the unit normal n of each vertex q on O
    mesh.add_attribute("vertex_normal")
    mesh_normal = mesh.get_attribute("vertex_normal")

    # 3 Place all vertices in a maximum heap H keyed on the [g(...)^2]
    # measure - the vertex with the maximum [g(...)^2] is placed at the top
    # of H
    g_2 = []
    for i in range(len(list_g)):
        g_2.append(list_g[i]*list_g[i])

    #4 -> #11
    j = 0
    delta_zero = 0.000001
    delta = delta_zero
    h = 0.01

    g_2_vertex_max = max(g_2)
    idx_vertex_max = g_2.index(g_2_vertex_max)

    while(j< N_max and g_2_vertex_max > epsilon):
        print("j = ", j)

        #Updating the normal vectors
        mesh.add_attribute("vertex_normal")
        mesh_normal = mesh.get_attribute("vertex_normal")

        mesh.add_attribute("vertex_gaussian_curvature")
        mesh_gauss = mesh.get_attribute("vertex_gaussian_curvature")
        print("idx_vertex_max before =", idx_vertex_max, ", g_2_vertex_max before=", g_2_vertex_max)
        print("Gaussian curvature before =", mesh_gauss[idx_vertex_max])

        # Moving vertex_max by delta*n(vertex_max) according to equation 17
        T_delta = T(mesh, mesh_normal, border_vertices_list, delta, idx_vertex_max)
        T_point_delta = T_point(mesh, mesh_normal, border_vertices_list, delta, idx_vertex_max, h)
        print("T(delta) =", T_delta, ", T_point(delta) =", T_point_delta)
        if ((abs(T_point_delta) > 0.1) and (abs(T_delta/T_point_delta) < 0.1)):
            delta = delta - T_delta/T_point_delta
            print("delta =", delta)

            #We have to copy the vertices to modify it
            vertices_copy = mesh.vertices.copy()
            idx_vertex_max = int(idx_vertex_max)
            delta_normal = [delta * mesh_normal[3*idx_vertex_max], delta * mesh_normal[3*idx_vertex_max + 1], delta * mesh_normal[3*idx_vertex_max+2]]
            vertices_copy[idx_vertex_max] = mesh.vertices[idx_vertex_max] + delta_normal
            print("vertex_max position =", mesh.vertices[idx_vertex_max])
            print("delta * vertex_max_normal =", delta_normal)
            print("new position =", vertices_copy[idx_vertex_max])

            #Updating the mesh
            mesh = pymesh.form_mesh(vertices_copy, mesh.faces)
            mesh.enable_connectivity()

            #Compute the new g_2 for the vertex and its adjacent vertices
            g_2[idx_vertex_max] = g(mesh, idx_vertex_max, border_vertices_list) **2
            list_idx_adjacent_vertices = mesh.get_vertex_adjacent_vertices(idx_vertex_max)
            for qi in list_idx_adjacent_vertices:
                g_2[qi] = g(mesh, qi, border_vertices_list) ** 2

        else:
            #If T_point is too close to 0, we switch vertex
            g_2[idx_vertex_max] = 0.0

        mesh.add_attribute("vertex_gaussian_curvature")
        mesh_gauss = mesh.get_attribute("vertex_gaussian_curvature")
        print("Gaussian curvature after =", mesh_gauss[idx_vertex_max])

        #Updating the g_2 max and its vertex
        g_2_vertex_max = max(g_2)
        idx_vertex_max = g_2.index(g_2_vertex_max)
        print("idx_vertex_max after =", idx_vertex_max, ", g_2_vertex_max after =", g_2_vertex_max)
        print("\n")

        name = "test_mesh_" + str(j) + ".obj"
        pymesh.save_mesh("./test/" + name, mesh)
        j+=1

if __name__ == "__main__":
    main()
