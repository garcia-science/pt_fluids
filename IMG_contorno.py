import cv2
import os
import numpy as np
from matplotlib import pyplot as plt
from PIL import Image

def auto_canny(image, sigma):
	v = np.median(image)
	lower = int(max(0, (1.0 - sigma) * v))
	upper = int(min(255, (1.0 + sigma) * v))
	edged = cv2.Canny(image, lower, upper)
	return edged
IMGs = os.listdir('C:\IMG') #lista de nombres de archivos en la carpeta indicada
sigma = 6
if __name__ == '__main__':
    im = cv2.imread("C:\IMG\cam0001.jpg")
    fromCenter = False
    r = cv2.selectROI(im)
    r = list(r)
    imCrop = im[r[1]:(r[1]+r[3]), r[0]:(r[0]+r[2])]
    imBlur = cv2.GaussianBlur(imCrop, (9,9), 6, 6, cv2.BORDER_DEFAULT)
    #edges = cv2.Canny(imBlur,10,200)
    edges = auto_canny(imBlur,sigma)
    pixels = np.asarray(edges)
    cv2.imwrite(os.path.join("C:\IIMG", IMGs[0]), edges)
    for i in range(1,len(IMGs)):
        im = cv2.imread('C:\IMG'+'\\'+IMGs[i])
        imCrop = im[r[1]:(r[1]+r[3]), r[0]:(r[0]+r[2])]
        imBlur = cv2.GaussianBlur(imCrop, (5,5), 6, 6, cv2.BORDER_DEFAULT)
        #edges = cv2.Canny(imBlur,10,200)
        edges = auto_canny(imBlur,sigma)
        pix = np.asarray(edges)
        pixels = np.concatenate((pixels, pix), 0)
        cv2.imwrite(os.path.join("C:\IIMG", IMGs[i]), edges)
    np.savetxt("C:\IIMG\pixel_data.csv", pixels, delimiter=",")
