from PIL import Image
import numpy as np
from matplotlib import pyplot as plt

#-----read test.png-----
filename = 'test_tilt.png'
img = Image.open(filename)
gray_img = img.convert('L')
f_xy = np.asarray(gray_img)
pil_img = Image.fromarray(f_xy)
#pil_img.show()
#print(img.mode) # RGBA
#rgb = img.convert('RGB')

#-----FFT-----
f_uv = np.fft.fft2(f_xy)
shifted_f_uv = np.fft.fftshift(f_uv)
ps = 20 * np.log(np.abs(shifted_f_uv))
plt.imshow(ps)
#plt.show()

#-----angle-----
print("size:", ps.shape)
print("row:", ps.shape[0])
print("column:", ps.shape[1])
nrow = int(ps.shape[0]) # x axis
ncolumn = int(ps.shape[1]) # y axis
#print (ps)
#
rcut = int(min(nrow,ncolumn)/2)
print ("FFT rcut:" ,rcut)
xo = int(nrow/2)
yo = int(ncolumn/2)
print ("FFT center: ",xo,yo)
#
max_evalue = 0.0
aim_angle = 0.0
for i in range(4500):
    angle = float(i)/100.0
    evalue = 0.0
    for r in range(int(rcut/2.0),rcut):
        x = int(r*np.cos(np.radians(angle)))
        y = int(r*np.sin(np.radians(angle)))
        evalue = evalue + ps[x+xo][y+yo] + ps[-y+xo][x+yo]
    if (evalue > max_evalue):
        max_evalue = evalue
        aim_angle = angle
    #print (angle, evalue)
#print (aim_angle,max_evalue)
#
for r in range(rcut):
   x = int(r*np.cos(np.radians(aim_angle)))
   y = int(r*np.sin(np.radians(aim_angle)))
   ps[x+xo][y+yo] = 0
   ps[-x+xo][-y+yo] = 0
   ps[-y+xo][x+yo] = 0
   ps[y+xo][-x+yo] = 0
plt.imshow(ps)
plt.show()

#-----rotate-----
gray_img = img.convert('L').rotate(-aim_angle, expand = 1, fillcolor = 255)
f_xy = np.asarray(gray_img)
pil_img = Image.fromarray(f_xy)
pil_img.show()
pil_img.save("case.png")

#-----resoluve-----
def resolve(f_xy,ncolumn_start,ncolumn_end,nrow_start,nrow_end,thresold):
   blunk_row_sum = 0
   for xr in range(ncolumn_start,ncolumn_end):
      blunk_row_sum = blunk_row_sum + f_xy[0][xr]
   blunk_column_sum = 0
   for yr in range(nrow_start,nrow_end):
      blunk_column_sum = blunk_column_sum + f_xy[yr][0]
   #
   nf_xy = np.copy(f_xy)
   nf_xy[:][:] = 255
   #
   for xr in range(ncolumn_start,ncolumn_end):
      column_sum = 0
      for yr in range(nrow_start,nrow_end):
         column_sum = column_sum + f_xy[yr][xr]
      if (column_sum < blunk_column_sum*thresold):
         for yrc in range(nrow_start,nrow_end):
            nf_xy[yrc][xr] = 1
   #
   for yr in range(nrow_start,nrow_end):
      row_sum = 0
      for xr in range(ncolumn_start,ncolumn_end):
         row_sum = row_sum + f_xy[yr][xr]
      if (row_sum < blunk_row_sum*thresold):
         for xrc in range(ncolumn_start,ncolumn_end):
            nf_xy[yr][xrc] = nf_xy[yr][xrc] + 1
   #
   for xr in range(ncolumn):
      for yr in range(nrow):
         if (nf_xy[yr][xr] != 2):
            nf_xy[yr][xr] = 255
   #
   return nf_xy

rnf_xy = np.copy(f_xy)
rnf_xy[:][:] = 255
# thresold = 0.5: main graph
nf_xy = resolve(f_xy,0,ncolumn,0,nrow,0.5)
xmg = np.full(100,-ncolumn)
ymg = np.full(100,-nrow)
i = 0
for xr in range(ncolumn):
   for yr in range(nrow):
      if (nf_xy[yr][xr] == 2):
         xmg[i] = xr
         ymg[i] = yr
         i = i + 1
         rnf_xy[yr][xr] = 125
#print (xmg, ymg)
x_start = min(np.abs(xmg))
y_last = min(np.abs(ymg))
x_last = max(xmg)
y_start = max(ymg)
x_center = int((x_start + x_last)/2.0)
y_center = int((y_start + y_last)/2.0)
print ("X, start, center, end:",x_start,x_center,x_last)
print ("Y, start, center, end:",y_start,y_center,y_last)
#
# thresold = 1.0: units and ytics
nf_xy = resolve(f_xy,0,x_start,0,y_start,1.0)
#for xr in range(ncolumn):
#   for yr in range(nrow):
#      if (nf_xy[yr][xr] == 2):
#         rnf_xy[yr][xr] = 125
#
# y unit and tics
yt = np.empty(x_start)
yu = np.full(x_start,0)
i = 0
for xr in range(1,x_start):
   yt[xr] = np.abs(nf_xy[y_center][xr] - nf_xy[y_center][xr-1])
   if (yt[xr] != 0):
      yu[i] = xr
      i = i + 1
yu[i] = x_start
y_unit_start = yu[0]
y_unit_last  = yu[1]
y_unit_center_x = int((yu[0]+yu[1])/2.0)
y_tics_start = yu[2]
y_tics_last  = yu[3]
y_tics_center_x = int((yu[2]+yu[3])/2.0)
print ("y, unit", y_unit_start, y_unit_center_x, y_unit_last)
y_unit_size_hight = int(y_unit_last - y_unit_start)
print ("y, unit size hight", y_unit_size_hight)
print ("y, tics", y_tics_start, y_tics_center_x, y_tics_last)
print ("y, tics size", (y_tics_last - y_tics_start))
#
nf_xy = resolve(f_xy,y_unit_start,y_unit_last,0,nrow,1.0)
for xr in range(y_unit_start,y_unit_last):
   for yr in range(nrow):
      if (nf_xy[yr][xr] == 2):
         rnf_xy[yr][xr] = 125
#
yty = np.empty(nrow)
yuy = np.full(nrow,-nrow)
i = 0
for yr in range(1,nrow):
   yty[xr] = np.abs(nf_xy[yr][y_unit_center_x] - nf_xy[yr-1][y_unit_center_x])
   if (yty[xr] != 0):
      yuy[i] = yr
      i = i + 1
y_unit_start_y = max(yuy)
y_unit_last_y  = min(np.abs(yuy))
print (y_unit_start_y, y_unit_last_y)
y_unit_size_width = int(y_unit_start_y - y_unit_last_y)
print ("unit size width",y_unit_size_width)
predit_nchar = int(y_unit_size_width/(y_unit_size_hight*0.5))
print ("predict number of characters",predit_nchar)
#
nf_xy = resolve(f_xy,y_tics_start,y_tics_last-5,0,nrow,1.0)
for xr in range(y_tics_start,y_tics_last-5):
   for yr in range(nrow):
      if (nf_xy[yr][xr] == 2):
         rnf_xy[yr][xr] = 125
#
# thresold = 1.0: units and xtics
nf_xy = resolve(f_xy,x_start,ncolumn,y_start,nrow,1.0)
#for xr in range(ncolumn):
#   for yr in range(nrow):
#      if (nf_xy[yr][xr] == 2):
#         rnf_xy[yr][xr] = 125
#
# x unit and tics
yutr = nrow - y_start
xt = np.empty(yutr+1)
xu = np.full(yutr+1,0)
i = 0
for yr in range(yutr,1,-1):
   xt[yr] = np.abs(nf_xy[y_start+yr][x_center] - nf_xy[y_start+yr-1][x_center])
   if (xt[yr] != 0):
      xu[i] = y_start+yr
      i = i + 1
xu[i] = y_start
x_tics_last  = xu[2]
x_tics_start = xu[3]
x_tics_center = int((xu[2]+xu[3])/2.0)
x_unit_last  = xu[0]
x_unit_start = xu[1]
x_unit_center_y = int((xu[0]+xu[1])/2.0)
print ("x, tics", x_tics_start, x_tics_center, x_tics_last)
print ("x, tics size", (x_tics_last - x_tics_start))
print ("x, unit", x_unit_start, x_unit_center_y, x_unit_last)
x_unit_size_hight = int(x_unit_last - x_unit_start)
print ("x, unit size", x_unit_size_hight)
#
nf_xy = resolve(f_xy,0,ncolumn,x_tics_start,x_tics_last,1.0)
for xr in range(ncolumn):
   for yr in range(x_tics_start,x_tics_last):
      if (nf_xy[yr][xr] == 2):
         rnf_xy[yr][xr] = 125
#
nf_xy = resolve(f_xy,0,ncolumn,x_unit_start,x_unit_last,1.0)
for xr in range(ncolumn):
   for yr in range(x_unit_start,x_unit_last):
      if (nf_xy[yr][xr] == 2):
         rnf_xy[yr][xr] = 125
#
xtx = np.empty(ncolumn)
xux = np.full(ncolumn,-ncolumn)
i = 0
for xr in range(1,ncolumn):
   xtx[xr] = np.abs(nf_xy[x_unit_center_y][xr] - nf_xy[x_unit_center_y][xr-1])
   if (xtx[xr] != 0):
      xux[i] = xr
      i = i + 1
x_unit_start_x = min(np.abs(xux))
x_unit_last_x  = max(xux)
print (x_unit_start_x, x_unit_last_x)
x_unit_size_width = int(x_unit_last_x - x_unit_start_x)
print ("unit size width",x_unit_size_width)
predit_nchar = int(x_unit_size_width/(x_unit_size_hight*0.5))
print ("predict number of characters",predit_nchar)
#
for xr in range(ncolumn):
   for yr in range(nrow):
      if (f_xy[yr][xr] != 255):
         rnf_xy[yr][xr] = f_xy[yr][xr]
#
pil_img = Image.fromarray(rnf_xy)
pil_img.show()
pil_img.save("case.png")

#-----IFFT-----
#unshifted_f_uv = np.fft.fftshift(shifted_f_uv)
#i_f_xy = np.fft.ifft2(unshifted_f_uv).real
#plt.imshow(i_f_xy)
#plt.show()

#-----output-----
#new_data = (((data - np.min(data)) / (np.max(data) - np.min(data))) * 255).astype(np.uint8)
#new_pil_img = pil_img.convert("L")
#new_pil_img.show()
#new_pil_img.save("case.png")
