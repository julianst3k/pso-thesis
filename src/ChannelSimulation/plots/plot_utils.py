from numpy import genfromtxt
import numpy as np

def open_file(path, shape):
    arr = genfromtxt(path, delimiter=",")
    max_size = np.prod(shape[:-1])
    arr_placeholder = np.zeros(shape)
    cnt = len(shape)-1
    for i in range(shape[0]):
        if len(shape)==1:
            return arr
        else:
            arr_placeholder[i, :] = recursive_opening(shape[1:],arr,i)

    return arr_placeholder
def recursive_opening(shape, arr, j):
    placeholder_arr = np.zeros(shape)
    if len(shape)==1:
        return arr[j*shape[0]:(j+1)*shape[0]]
    else:
        for i in range(shape[0]):
            placeholder_arr[i,:] = recursive_opening(shape[1:], arr, i*j+i)
    return placeholder_arr
def do_exponential_form(arr, intensity, time):
    return np.exp(-intensity*arr*time)

def pablo_H(x, y, xt, yt, xr, yr, zj):
    y_ij = yt-yr
    x_ij = xt-xr
    dij = 2*np.sqrt(y_ij**2+x_ij**2)
    arr = (y_ij**2+x_ij**2+(x-xr)**2+(y-yr)**2-((x-xt)**2+(y-yt)**2))/dij+zj
    print(arr)
    return arr

def paper_H(x, y, xt, yt, xr, yr, zj, zi):
    y_ij = yt-yr
    x_ij = xt-xr
    dij = 2*np.sqrt(y_ij**2+x_ij**2)
    cot = 2*(zi-zj)/dij
    arr = (y_ij**2+x_ij**2+(x-xr)**2+(y-yr)**2-((x-xt)**2+(y-yt)**2))/dij*cot+zj
    print(arr)
    return arr

def my_H(x,y,xt,yt,xr,yr,zj,zi):
    y_ij = yt - yr
    x_ij = xt - xr
    dij = 2 * np.sqrt(y_ij ** 2 + x_ij ** 2)
    cot = 2 * (zi - zj) / dij
    try:
        alpha = y_ij/x_ij
    except:
        alpha = y_ij/(x_ij+0.001)
    arr = (y_ij ** 2 + x_ij ** 2 + ((y - yr)/alpha) ** 2 + (y - yr) ** 2 - (((y - yt)/alpha) ** 2 + (y - yt) ** 2)) / dij * cot + zj
    return arr
def H_wrapper(x, y, xt, yt, xr, yr, zj, zi):
    return pablo_H(x, y, xt, yt, xr, yr, zj), paper_H(x, y, xt, yt, xr, yr, zj, zi), my_H(x, y, xt, yt, xr, yr, zj, zi)