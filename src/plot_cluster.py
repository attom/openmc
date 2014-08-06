#! /usr/bin/env python

'''
Andrew Till
Summer 2014

Plot clustering from OpenMC
'''

#stdlib
import os
#tpl
import numpy as np
from matplotlib import pyplot as plt

def do_all():
    dirr = '.'
    files = os.listdir(dirr)
    nuclideNames = [name[4:13] for name in files if name.endswith('.txt') and name.startswith('obs')]
    for nuclideName in nuclideNames:
        filename = 'obs_{0}.txt'.format(nuclideName)
        filePath = os.path.join(dirr, filename)
        observations = np.loadtxt(filePath)
        numXS = observations.shape[1]
        #
        filename = 'cen_{0}.txt'.format(nuclideName)
        filePath = os.path.join(dirr, filename)
        centroids = np.loadtxt(filePath)
        if centroids.shape[0] == 10:
            continue
        #
        plt.figure(1)
        names = ['total', 'elastic scattering', 'fission']
        for dim1 in range(0, numXS):
            for dim2 in range(dim1+1, numXS):
                if dim1 == 1 and dim2 == 2:
                    dim1 = 2
                    dim2 = 1
                plt.clf()
                plt.loglog(observations[:,dim1], observations[:,dim2], marker='o', linestyle='', markersize=10, color='b', rasterized=True)
                plt.loglog(centroids[:,dim1], centroids[:, dim2], marker='s', linestyle='', markersize=10, color='r')
                plt.xlabel('{0} (b)'.format(names[dim1]))
                plt.ylabel('{0} (b)'.format(names[dim2]))
                plt.title('Clustering for {0}'.format(nuclideName))
                plt.axes().set_aspect('equal')
                filename = 'p_{0}_{1}_{2}.pdf'.format(nuclideName, dim1, dim2)
                filePath = os.path.join(dirr, filename)
                plt.savefig(filePath)

if __name__ == '__main__':
    do_all()
