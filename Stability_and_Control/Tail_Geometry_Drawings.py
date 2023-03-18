import matplotlib.pyplot as plt
import numpy as np

def Geometry_Drawing(crh,cth,bh,Sweep_LE_HT,Sweep_half_HT,real):
    x1 = 0
    y1 = crh / 2
    x4 = 0
    y4 = -crh / 2
    y2 = crh / 2 - (bh / 2) * np.tan(Sweep_LE_HT)
    if Sweep_half_HT < 0:
        y3 = -0.5 * cth + (bh / 2) * np.tan(Sweep_half_HT)
    else:
        y3 = -0.5 * cth - (bh / 2) * np.tan(Sweep_half_HT)
    if real ==0:
        x2 = -bh / 2
        x3 = -bh / 2
    else:
        x2 = bh/2
        x3 = bh/2
    if real ==0:
        xupper = np.linspace(x2, x1)
        yupper = -((y1 - y2) / x2) * xupper + y1
        xlower = np.linspace(x3, x4)
        ylower = ((y3 - y4) / x3) * xlower + y4
    else:
        xupper = np.linspace(x1, x2)
        yupper = -((y1 - y2) / x2) * xupper + y1
        xlower = np.linspace(x4, x3)
        ylower = ((y3 - y4) / x3) * xlower + y4
    return x1,x2,y1,y2,y3,y4,xupper,yupper,xlower,ylower

def Geometry_elevator():
    return x1_e, x2_e, y1_e, y2_e, y3_e, y4_e, xupper_e, yupper_e, xlower_e, ylower_e

def Geometry_trimtab():
    return x1_t, x2_t, y1_t, y2_t, y3_t, y4_t, xupper_t, yupper_t, xlower_t, ylower_t

def Plot_Graph(crh, cth,bh,Sweep_LE_HT,Sweep_half_HT,real, colour1,colour2, trimtab, elevator):
    x1, x2, y1, y2, y3, y4, xupper, yupper, xlower, ylower = Geometry_Drawing(crh,cth,bh,Sweep_LE_HT,Sweep_half_HT,real)

    if elevator == 0:
        x1_e, x2_e, y1_e, y2_e, y3_e, y4_e, xupper_e, yupper_e, xlower_e, ylower_e = Geometry_elevator()
        plt.plot(xupper_e, yupper_e, color=colour1)
        plt.plot(xlower_e, ylower_e, color=colour1)
        plt.vlines(x=x2_e, ymin=y3_e, ymax=y2_e, color=colour1)
        plt.vlines(x=x1_e, ymin=y4_e, ymax=y1_e, color=colour1)

    if trimtab ==0:
        x1_t, x2_t, y1_t, y2_t, y3_t, y4_t, xupper_t, yupper_t, xlower_t, ylower_t = Geometry_trimtab()
        plt.plot(xupper_t, yupper_t, color=colour1)
        plt.plot(xlower_t, ylower_t, color=colour1)
        plt.vlines(x=x2_t, ymin=y3_t, ymax=y2_t, color=colour1)
        plt.vlines(x=x1_t, ymin=y4_t, ymax=y1_t, color=colour1)

    plt.plot(xupper, yupper, color=colour1)
    plt.plot(xlower, ylower, color=colour1)
    plt.vlines(x=x2, ymin=y3, ymax=y2, color=colour1)
    plt.vlines(x=x1, ymin=y4, ymax=y1, color=colour1)
    plt.fill_between(xlower, ylower, yupper, color=colour2)
    return

Plot_Graph(crh=4, cth=3,bh=12,Sweep_LE_HT=10*np.pi/180,Sweep_half_HT=0,real=0,colour1="#11aa00",colour2="#11aa0055", trimtab=0,elevator=0)
Plot_Graph(crh=5, cth=3,bh=14,Sweep_LE_HT=10*np.pi/180,Sweep_half_HT=0,real=1, colour1="#728FCE", colour2="#98AFC7", trimtab=1,elevator=1)
Plot_Graph(crh=7, cth=4,bh=14,Sweep_LE_HT=10*np.pi/180,Sweep_half_HT=0,real=1, colour1="#AFDCEC", colour2="#DCE9FA",trimtab=1,elevator=1)
plt.show()