import os
cd= os.getcwd()


from skimage.measure import compare_ssim
#from skimage.metrics import structural_similarity
#import argparse
import imutils
import cv2
import matplotlib.pyplot as plt
import copy
import matplotlib.animation as animation
import time

'''
Parte para correrlo desde la terminal:
    
construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-f", "--first", required=True,
	help="first input image")
ap.add_argument("-s", "--second", required=True,
	help="second")
args = vars(ap.parse_args())
'''


#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=10,metadata=dict(artist='Me'), bitrate=1800)


# load the two input images
#img_name_1 = cd + '/sss/ss6.png'
#img_name_2 = cd + '/sss/ss7.png'
#img_name_3 = cd + '/sss/ss8.png'
ims1 = []
ims2 = []
ims3 = []
ims4 = []
fig1 = plt.figure()
ax = fig1.add_subplot(223)

cont=0

names=cd + '\\diff_ani4\\frame'
scores = []
scores_X = []

def my_plot(ax,listay,listax):
    p0, = ax.plot(listax)
    p1, = ax.plot(listay)
    return [p1, p0]   # return a list of the new plots

start = time.time()
while len(ims2)<len(os.listdir(cd+'\\diff_ani4'))-1:
    t0=time.time()
    imageA = cv2.imread(names + str(cont) + '.jpg')
    imageB = cv2.imread(names + str(cont+1) + '.jpg')

    # convert the images to grayscale
    grayA = cv2.cvtColor(imageA, cv2.COLOR_BGR2GRAY)
    grayB = cv2.cvtColor(imageB, cv2.COLOR_BGR2GRAY)

    # compute the Structural Similarity Index (SSIM) between the two
    # images, ensuring that the difference image is returned

    (score, diff) = compare_ssim(grayA, grayB, full=True)
    #(score, diff) = structural_similarity(grayA, grayB, full=True)

    diff = (diff * 255).astype("uint8")
    scores.append(score)
    scores_X.append(len(ims2)+1)
    print("SSIM: {}".format(score))

    # threshold the difference image, followed by finding contours to
    # obtain the regions of the two input images that differ
    thresh = cv2.threshold(diff, 0, 255,
    	cv2.THRESH_BINARY_INV | cv2.THRESH_OTSU)[1]
    cnts = cv2.findContours(thresh.copy(), cv2.RETR_EXTERNAL,
    	cv2.CHAIN_APPROX_SIMPLE)
    cnts = imutils.grab_contours(cnts)


    plt.subplot(224)
    im3 = plt.imshow(copy.copy(imageB),animated=True,interpolation='nearest',clim=[0,255])
    ims3.append([im3])


    # loop over the contours
    for c in cnts:
    	# compute the bounding box of the contour and then draw the
    	# bounding box on both input images to represent where the two
    	# images differ
    	(x, y, w, h) = cv2.boundingRect(c)
    	cv2.rectangle(imageA, (x, y), (x + w, y + h), (0, 0, 255), 2)
    	cv2.rectangle(imageB, (x, y), (x + w, y + h), (0, 0, 255), 2)

    plt.subplot(221)
    im2 = plt.imshow(copy.copy(diff),animated=True,interpolation='nearest',clim=[0,255])
    ims2.append([im2])
    plt.subplot(222)
    im1 = plt.imshow(copy.copy(imageB),animated=True,interpolation='nearest',clim=[0,255])
    ims1.append([im1])



    p0, = ax.plot(scores,color='blue')

    #ps = my_plot(ax,scores,scores_X)
    ims4.append([p0])


    #im4 = plt.plot(copy.copy(scores))
    #ims4.append([im4])


    cont+=1
    tf=time.time()
    tt=tf-t0
    print(f'tiempo del loop = {tt}')


#fig1 = plt.figure()
ani1 = animation.ArtistAnimation(fig1, ims1, interval=350, blit=True,repeat_delay=1)#interval ~ tiempo entre frames
ani2 = animation.ArtistAnimation(fig1, ims2, interval=350, blit=True,repeat_delay=1)#interval ~ tiempo entre frames
ani3 = animation.ArtistAnimation(fig1, ims3, interval=350, blit=True,repeat_delay=1)#interval ~ tiempo entre frames
ani4 = animation.ArtistAnimation(fig1, ims4, interval=350, blit=True,repeat_delay=1)#interval ~ tiempo entre frames

end = time.time()
print(f'tiempo total = {end-start}')

#ani1.save('box.mp4' , writer=writer)
#ani2.save('diff.mp4' , writer=writer)

#plt.figure()
#plt.plot(scores)

## show the output images
##cv2.imshow("Original", imageA)
##cv2.imshow("Modified", imageB)
#plt.figure()
#plt.imshow(diff)
##cv2.imshow("Diff", diff)
#plt.figure()
#plt.imshow(thresh)
##cv2.imshow("Thresh", thresh)
##cv2.waitKey(0)
