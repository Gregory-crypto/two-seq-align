
import numpy as np
import pandas as pd

def align(seqi, seqii):

    seqi = np.array(list(seqi))
    seqii = np.array(list(seqii))


    # parameters

    match = 2
    mismatch = -1
    open_gap = -2
    gap = -1
    ni = len(seqi)
    nii = len(seqii)

    matrix = [[2, -6, -6, -6], [-6, 2, -6, -6], 
              [-6, -6, 2, -6], [-6, -6, -6, 2]]


    M = np.zeros((ni + 1, nii + 1))
    M[0, 1] += open_gap
    M[1, 0] += open_gap


    M[1:, 0] = np.linspace(-2, -2 + (ni - 1) * gap, ni)
    M[0, 1:] = np.linspace(-2, -2 + (nii - 1) * gap, nii)

    # computing NW matrix

    for i in range(1, ni + 1):
        for j in range(1, nii + 1):
            diag = M[i - 1, j - 1]
            ver = M[i - 1, j]
            hor = M[i, j - 1]
            diag += match if (seqi[i - 1] == seqii[j - 1]) else mismatch
            ver += open_gap if (i == 1) else gap
            hor += open_gap if (j == 1) else gap
            M[i, j] = max([diag, ver, hor])


    # finding the optimal way

    al_seqi = []
    al_seqii = []
    al_seqi.append(seqi[-1])
    al_seqii.append(seqii[-1])


    i = ni - 1
    j = nii - 1


    while i > 0 and j > 0:
        diag = M[i - 1, j - 1]
        ver = M[i - 1, j]
        hor = M[i, j - 1]
        


        if diag >= ver and diag >= hor:
            i -= 1
            j -= 1
            al_seqi.append(seqi[i])
            al_seqii.append(seqii[j])

        elif hor > diag and hor > ver:
            j -= 1
            al_seqii.append(seqii[j])
            al_seqi.append('-')


        elif ver > diag and ver > hor:
            i -= 1
            al_seqi.append(seqi[i])
            al_seqii.append('-')
            
            
        
    al_seqi = ''.join(al_seqi)[::-1]
    al_seqii = ''.join(al_seqii)[::-1]
            

    # Columns & rows

    return [al_seqi, al_seqii]
