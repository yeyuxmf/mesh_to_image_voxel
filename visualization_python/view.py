# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import os
import numpy as np
import cv2

def read_txt_view_render_image():

    img = np.zeros((1024, 1024), np.float)
    data_txt = open("./render_image.txt", "r")
    for line in data_txt.readlines():
        line_  = line.strip().split(" ")
        img[int(line_[0]), int(line_[1])] = float(line_[2])
    #img = img / np.max(img)
    cv2.namedWindow("img", cv2.WINDOW_NORMAL)
    cv2.imshow("img", img)
    cv2.waitKey(0)

def read_txt_view_depth_image():

    img = np.zeros((1024, 1024), np.float)
    data_txt = open("./depth_image.txt", "r")
    for line in data_txt.readlines():
        line_ = line.strip().split(" ")
        img[int(line_[0]), int(line_[1])] = float(line_[2])

    nonv = np.nonzero(img)[1]  #normalization

    img = (img  - np.min(nonv))/ (np.max(nonv) - np.min(nonv))
    cv2.namedWindow("img", cv2.WINDOW_NORMAL)
    cv2.imshow("img", img)
    cv2.waitKey(0)



if __name__ == '__main__':

    read_txt_view_render_image()
    read_txt_view_depth_image()


