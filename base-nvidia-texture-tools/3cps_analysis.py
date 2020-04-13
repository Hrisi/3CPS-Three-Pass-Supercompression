import warnings
from subprocess import Popen, PIPE
import argparse
import os
warnings.filterwarnings('ignore', module='skimage')
import shutil
import cv2
import numpy as np
import math

BASE_PATH = "/nvidia-texture-tools/data"

# ======================================

def init_folder(nameFolder, erase):
   """ cerate the folder

   Args:
      nameFolder (str): it is the folder that should be created

      erase (bool): erase everything in the folder if it exists

   Returns:
      Nan
   """
   exist = os.path.isdir(nameFolder)

   if exist and erase:
      shutil.rmtree(nameFolder)

   if not exist or erase:  
      os.mkdir(nameFolder)

   # +++++++++++++++++++++++++++++++++++++++++++++++++++++

def extend_w_suffix(img, suffix):
      """ extend the corename with the suffix

      Args:
         img (str): the image
         suffix (str): the suffix

      Returns:
         filename (str): the full name
      """
      
      coreImage = ".".join(img.split(".")[:-1])

      return ".".join([coreImage, suffix])

   # +++++++++++++++++++++++++++++++++++++++++++++++++++++

def add_path(folder, filename):
      """ create the path

      Args:
         folder (str): the right folder

         filename (str): the file name

      Returns:
         fullpath (str): complete path
      """

      return os.path.join(folder, filename)

   # +++++++++++++++++++++++++++++++++++++++++++++++++++++
 
def psnr(img_name1, img_name2):
      """ compute psnr

      Args:
         img_name1 (str): the right folder

         img_name2 (str): the file name

      Returns:
         psnr (float)
      """

      img1 = cv2.imread(img_name1)
      img2 = cv2.imread(img_name2)

      mse = np.mean( (img1 - img2) ** 2 )
      if mse == 0:
          return 100
      PIXEL_MAX = 255.0
      return 20 * math.log10(PIXEL_MAX / math.sqrt(mse))

   # +++++++++++++++++++++++++++++++++++++++++++++++++++++
 
def compress_dxt(img, is_first_time):
      """ compress the image into DXT format

      Args:
         img (string): the name of the image to compress
         is_first_time (bool): True if first, False if second

      Returns:
         tCompress (float): the time it took to compress the image
      """

      original_image = add_path(BASE_PATH+"/webp", extend_w_suffix(img, "png"))
      if (is_first_time):
         destination_image = add_path(BASE_PATH+"/compressed_textures/dxt1", extend_w_suffix(img, "dds"))
      else:
         destination_image = add_path(BASE_PATH+"/compressed_textures/3cps", extend_w_suffix(img, "dds"))

      cmd = ["nvcompress", "-bc1"]
      cmd.extend([original_image])
      cmd.extend([destination_image])

      proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
      log, _ = proc.communicate()
      log_utf = log.decode('UTF-8')
      timeCompress = float(log_utf.split(" ")[-2])

      return timeCompress

   # +++++++++++++++++++++++++++++++++++++++++++++++++++++

def decompress_dxt(img, is_first_time):
      """decompress the DXT image

      Args:
         img (string): the name of the dds compressed image
         is_first_time (bool): True if first, False if second

      Return:
         psnr (float): the quality of the image with regards to the texture
      """

      ori = add_path(BASE_PATH+"/textures", extend_w_suffix(img, "jpg"))
      if is_first_time:
         dxt_image = add_path(BASE_PATH+"/compressed_textures/dxt1", extend_w_suffix(img, "dds"))
         tga_image = add_path(BASE_PATH+"/compressed_textures/dxt1/png", extend_w_suffix(img, "tga"))
         png_image = add_path(BASE_PATH+"/compressed_textures/dxt1/png", extend_w_suffix(img, "png"))

      else:
         dxt_image = add_path(BASE_PATH+"/compressed_textures/3cps", extend_w_suffix(img, "dds"))
         tga_image = add_path(BASE_PATH+"/compressed_textures/3cps/png", extend_w_suffix(img, "tga"))
         png_image = add_path(BASE_PATH+"/compressed_textures/3cps/png", extend_w_suffix(img, "png"))

      cmd = ["nvdecompress", dxt_image, tga_image]
      proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
      proc.communicate()

      cmd = ["convert", "-quality", "100"]
      cmd.extend([tga_image])
      cmd.extend([png_image])
      res = Popen(cmd)
      res.communicate()

      return psnr(ori, png_image)

# ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++

if __name__ == "__main__":

   parser = argparse.ArgumentParser(description='compress to DXT files')

   parser.add_argument('-i', dest='img',
         type = str,
         help = "the webp compressed image name")
   parser.add_argument('-f', dest='is_first',
         type = int,
         help = "indicator of whether it is the first compression or not")
   
   arg = parser.parse_args()

   init_folder(BASE_PATH+"/compressed_textures", False)

   init_folder(BASE_PATH+"/compressed_textures/dxt1", False)
   init_folder(BASE_PATH+"/compressed_textures/3cps", False)

   init_folder(BASE_PATH+"/compressed_textures/dxt1/png", False)
   init_folder(BASE_PATH+"/compressed_textures/3cps/png", False)

   # compression
   tReCompress = compress_dxt(arg.img, arg.is_first)
   psnr = decompress_dxt(arg.img, arg.is_first)
   print ("Texture ", arg.img)
   print ("Time for compression: ", tReCompress, "seconds")
   print("PSNR quality of compressed texture: ", psnr, "dB")
