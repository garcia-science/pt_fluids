import cv2
import os
import numpy as np

def auto_canny(image, sigma):
    v = np.median(image)
    lower = int(max(0, (1.0 - sigma) * v))
    upper = int(min(255, (1.0 + sigma) * v))
    edged = cv2.Canny(image, lower, upper)
    return edged


def deteccion(main_path, file_i, file_o, REC, sigma):
    IMGs = os.listdir(main_path + file_i)  # lista de nombres de archivos en la carpeta indicada
    im = cv2.imread(main_path + file_i + '\cam0001.jpg')
    rec = list(REC)
    imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
    imBlur = cv2.GaussianBlur(imCrop, (9, 9), 6, 6, cv2.BORDER_DEFAULT)
    # edges = cv2.Canny(imBlur,10,200)
    edges = auto_canny(imBlur, sigma)
    cv2.imwrite(os.path.join(main_path + file_o, IMGs[0]), edges)
    for i in range(1, len(IMGs)):
        im = cv2.imread(main_path + file_i + '\\' + IMGs[i])
        imCrop = im[rec[1]:(rec[1] + rec[3]), rec[0]:(rec[0] + rec[2])]
        imBlur = cv2.GaussianBlur(imCrop, (5, 5), 6, 6, cv2.BORDER_DEFAULT)
        # edges = cv2.Canny(imBlur,10,200)
        edges = auto_canny(imBlur, sigma)
        cv2.imwrite(os.path.join(main_path + file_o, IMGs[i]), edges)
    return IMGs


def ROI_select(path):
    im = cv2.imread(path + '\cam0001.jpg')
    fromCenter = False
    RECs = cv2.selectROI(im)
    return RECs
