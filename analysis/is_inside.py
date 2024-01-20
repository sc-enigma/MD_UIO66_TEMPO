import numpy as np

def is_inside(verts, pt):
    center =  np.zeros(3)
    for vert in verts:
        center += vert / len(verts)
    res = True
    for comb in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]:
        norm = np.cross(verts[comb[1]] - verts[comb[0]], verts[comb[2]] - verts[comb[0]])
        res &= np.dot(norm, pt - verts[comb[0]]) * np.dot(norm, center - verts[comb[0]]) > 0.0
    return res
        

verts = []
verts.append(np.array([0.0, 0.0, 0.0]))
verts.append(np.array([1.0, 0.0, 0.0]))
verts.append(np.array([0.0, 1.0, 0.0]))
verts.append(np.array([0.0, 0.0, 1.0]))

pt = np.array([0.25, 0.25, 0.25])

print(is_inside(verts, pt))