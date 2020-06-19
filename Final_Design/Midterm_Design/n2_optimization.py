import numpy as np
from itertools import permutations
from tqdm import tqdm
import scipy.sparse as sparse
from matplotlib import pyplot as plt


def create_matrix(keys, data):
    mat = np.zeros((len(keys), len(keys)), dtype=int)
    for row, subsystem in enumerate(keys):
        for col, output_sys in enumerate(keys):
            if output_sys in data[subsystem]:
                mat[row, col] = 1
    return mat


def score_matrix(weights, matrix):
    # print(weights)
    FS = []
    BS = []
    BR = []
    FR = []
    for i, row in enumerate(matrix):
        # print("row: ", row)
        interfaces = np.sum(row[i+1:])
        # print(interfaces)
        cells = len(row)-i-1
        if interfaces == 0 or cells == 0:
            FS.append(0)
        else:
            # print("Here! (", cells-1, ",", interfaces-1, ") ", weights[cells-1, interfaces-1])
            FS.append(weights[cells-1, interfaces-1])
        interfaces = np.sum(row[:i])
        cells = i
        if interfaces == 0 or cells == 0:
            BS.append(0)
        else:
            BS.append(weights[cells-1, interfaces-1])
    for i, col in enumerate(matrix.T):
        interfaces = np.sum(col[i+1:])
        cells = len(row)-i-1
        if interfaces == 0 or cells == 0:
            BR.append(0)
        else:
            BR.append(weights[cells-1, interfaces-1])
        interfaces = np.sum(col[:i])
        cells = i
        if interfaces == 0 or cells == 0:
            FR.append(0)
        else:
            FR.append(weights[cells-1, interfaces-1])
    # print(FS)
    # print(FR)
    # print(BS)
    # print(BR)
    rowscore = np.asarray(FS) + np.asarray(BS) + np.asarray(BR) + np.asarray(FR)
    return np.sum(rowscore)


outputs = {1: [2,3,4,5,7,10,11], 2: [1,3,4,5,7,10,11], 3: [1,2,4], 4: [6,7,11],
        5: [1,2,4,6,7], 6: [9], 7: [1,2,8], 8: [9], 9: [3,4,5,10],
        10: [1,2,3,4,5], 11: [3,4,5,6,7,9,10]}

weights = np.array([[1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ],
                    [2  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ],
                    [3  ,3  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ],
                    [4  ,6  ,4  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ],
                    [5  ,10 ,10 ,5  ,1  ,0  ,0  ,0  ,0  ,0  ],
                    [6  ,15 ,20, 15 ,6  ,1  ,0  ,0  ,0  ,0  ],
                    [7  ,21 ,35 ,35 ,21 ,7  ,1  ,0  ,0  ,0  ],
                    [8  ,28 ,56 ,70 ,56 ,28 ,8  ,1  ,0  ,0  ],
                    [9  ,36 ,84 ,126,126,84 ,36 ,9  ,1  ,0  ],
                    [10 ,45 ,120,210,252,210,120,45 ,10 ,1  ]])

# bestkeys = list(outputs.keys())
# matrix = create_matrix(bestkeys, outputs)
# bestscore = score_matrix(weights, matrix)
# for keys in tqdm(permutations(list(outputs.keys()))):
#     matrix = create_matrix(keys, outputs)
#     score = score_matrix(weights, matrix)
#     if score > bestscore:
#         bestkeys = keys
#         bestscore = score
# print(bestscore)
# print(bestkeys)

matrix1 = create_matrix(list(outputs.keys()), outputs)
keys = [3, 1, 2, 5, 10, 11, 7, 4, 9, 8, 6]
matrix2 = create_matrix(keys, outputs)
plt.spy(matrix1, marker="o", markersize=10, color="b")
plt.spy(np.diag(np.ones(11)), marker="s", markersize=15, color="black")
ax = plt.gca()
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
plt.show()
plt.spy(matrix2, marker="o", markersize=10, color="b")
plt.spy(np.diag(np.ones(11)), marker="s", markersize=15, color="black")
ax = plt.gca()
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
plt.show()

