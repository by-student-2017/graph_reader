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
def resolve(f_xy,ncolumn,nrow,thresold):
   blunk_row_sum = 0
   for xr in range(ncolumn):
      blunk_row_sum = blunk_row_sum + f_xy[0][xr]
   blunk_column_sum = 0
   for yr in range(nrow):
      blunk_column_sum = blunk_column_sum + f_xy[yr][0]
   #
   nf_xy = np.copy(f_xy)
   nf_xy[:][:] = 255
   #
   for xr in range(ncolumn):
      column_sum = 0
      for yr in range(nrow):
         column_sum = column_sum + f_xy[yr][xr]
      if (column_sum < blunk_column_sum*thresold):
         for yrc in range(nrow):
            nf_xy[yrc][xr] = 1
   #
   for yr in range(nrow):
      row_sum = 0
      for xr in range(ncolumn):
         row_sum = row_sum + f_xy[yr][xr]
      if (row_sum < blunk_row_sum*thresold):
         for xrc in range(ncolumn):
            nf_xy[yr][xrc] = nf_xy[yr][xrc] + 1
   #
   return nf_xy

# thresold = 0.5: main graph
nf_xy = resolve(f_xy,ncolumn,nrow,0.5)
xmg = np.full(100,-ncolumn)
ymg = np.full(100,-nrow)
i = 0
for xr in range(ncolumn):
   for yr in range(nrow):
      if (nf_xy[yr][xr] == 2):
         xmg[i] = xr
         ymg[i] = yr
         i = i + 1
#print (xmg, ymg)
x_start = min(np.abs(xmg))
y_last = min(np.abs(ymg))
x_last = max(xmg)
y_start = max(ymg)
print ("X:",x_start,x_last)
print ("Y:",y_start,y_last)
#
# thresold = 1.0: units
#nf_xy = resolve(f_xy,ncolumn,nrow,0.98)
pil_img = Image.fromarray(nf_xy)
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
