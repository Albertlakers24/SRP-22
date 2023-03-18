import matplotlib.pyplot as plt
import numpy as np

def Geometry_Drawing(crh,cth,bh,Sweep_LE_HT,Sweep_half_HT,real):
    x1 = 0
    y1 = crh / 2
    x4 = 0
    y4 = -crh / 2
    y2 = cth / 2 - (bh / 2) * np.tan(Sweep_half_HT)
    if Sweep_half_HT < 0:
        y3 = -0.5 * cth + (bh / 2) * np.tan(Sweep_half_HT)
    else:
        y3 = -0.5 * cth - (bh / 2) * np.tan(Sweep_half_HT)
    if real ==0:
        x2 = -bh / 2
    else:
        x2 = bh/2
    if real ==0:
        xupper = np.linspace(x2, x1)
        yupper = -((y1 - y2) / x2) * xupper + y1
        xlower = np.linspace(x2, x4)
        ylower = ((y3 - y4) / x2) * xlower + y4
    else:
        xupper = np.linspace(x1, x2)
        yupper = -((y1 - y2) / x2) * xupper + y1
        xlower = np.linspace(x4, x2)
        ylower = ((y3 - y4) / x2) * xlower + y4
    return x1,x2,y1,y2,y3,y4,xupper,yupper,xlower,ylower

def Geometry_elevator(ylower,x_inboard,x_outboard,real,y3,y4,y1,y2,x2,cf_c,crh,cth,bh,w_horn):
    if real ==0:        #todo add horn
        x4_e = x_inboard
        x3_e = x_outboard
        x1_e = x4_e
        x2_e = x3_e
        cprime_i = ((cth-crh)/(-bh/2))*x_inboard+crh
        cprime_o = ((cth-crh)/(-bh/2))*x_outboard+crh
        c_i_e =cprime_i*cf_c
        c_o_e =cprime_o*cf_c
        cprime_horn = ((cth-crh)/(-bh/2))*(x2_e+w_horn)+crh
        c_o_horn=cprime_horn*cf_c
        ylower_i = ((y3 - y4) / x2) * x_inboard + y4
        ylower_o =((y3 - y4) / x2) * x_outboard + y4
        if ylower_i < 0:
            y1_e = ylower_i + c_i_e
        else:
            y1_e = ylower_i - c_i_e
        if ylower_o < 0:
            y2_e = ylower_o + c_o_e
        else:
            y2_e = ylower_o - c_o_e
        if ylower_o < 0:
            y2_horn = ylower_o + c_o_horn
        else:
            y2_horn= ylower_o - c_o_horn
        y3_e = ((y3 - y4) / x2) * x3_e + y4
        y4_e = ((y3 - y4) / x2) * x4_e + y4
        xupper_e = [x1_e,x2_e]
        yupper_e = [y1_e,y2_e]
        xlower_e = np.linspace(x3_e, x4_e)
        ylower_e = ((y3 - y4) / x2) * xlower_e + y4
    else:
        x4_e = x_inboard_e
        x3_e = x_outboard_e
        x1_e = x4_e
        x2_e = x3_e
        y1_e = -1.3
        y2_e = -0.8
        y3_e = ((y3 - y4) / x2) * x3_e + y4
        y4_e = ((y3 - y4) / x2) * x4_e + y4
        xupper_e = [x1_e, x2_e]
        yupper_e = [y1_e, y2_e]
        xlower_e = np.linspace(x3_e, x4_e)
        ylower_e = ((y3 - y4) / x2) * xlower_e + y4
    return x1_e, x2_e, y1_e, y2_e, y3_e, y4_e, xupper_e, yupper_e, xlower_e, ylower_e,y2_horn

def Plot_Graph(crh, cth,bh,Sweep_LE_HT,Sweep_half_HT,real, colour1,colour2, trimtab, elevator,x_inboard_e,x_outboard_e,x_inboard_t,x_outboard_t,ce_c,cf_c,horn,w_horn,h_horn):
    x1, x2, y1, y2, y3, y4, xupper, yupper, xlower, ylower = Geometry_Drawing(crh,cth,bh,Sweep_LE_HT,Sweep_half_HT,real)

    if elevator == 0:
        x1_e, x2_e, y1_e, y2_e, y3_e, y4_e, xupper_e, yupper_e, xlower_e, ylower_e,y2_horn = Geometry_elevator(ylower,x_inboard_e,x_outboard_e,real,y3,y4,y1,y2,x2,ce_c,crh,cth,bh,w_horn)
        if horn == 0:
            plt.plot(xlower_e, ylower_e, color="r")
            plt.vlines(x=x1_e, ymin=y4_e, ymax=y1_e, color="r")
            plt.vlines(x=x2_e, ymin=y3_e, ymax=y2_e+h_horn, color="r")
            plt.hlines(y=y2_e+h_horn,xmin=x2_e,xmax=x2_e+w_horn,color="r")
            plt.vlines(x=x2_e+w_horn, ymin=y2_horn, ymax=y2_e+h_horn, color="r")
            x_values = [x2_e+w_horn,x1_e]
            y_values = [y2_horn,y1_e]
            plt.plot(x_values,y_values,linestyle ="-", color="r")
        else:
            plt.plot(xupper_e, yupper_e, color="r")
            plt.plot(xlower_e, ylower_e, color="r")
            plt.vlines(x=x2_e, ymin=y3_e, ymax=y2_e, color="r")
            plt.vlines(x=x1_e, ymin=y4_e, ymax=y1_e, color="r")

    if trimtab ==0:
        x1_e, x2_e, y1_e, y2_e, y3_e, y4_e, xupper_e, yupper_e, xlower_e, ylower_e,y2_horn = Geometry_elevator(ylower,x_inboard_t,x_outboard_t,real,y3,y4,y1,y2,x2,cf_c,crh,cth,bh,w_horn)
        plt.plot(xupper_e, yupper_e, color="b")
        plt.plot(xlower_e, ylower_e, color="b")
        plt.vlines(x=x2_e, ymin=y3_e, ymax=y2_e, color="b")
        plt.vlines(x=x1_e, ymin=y4_e, ymax=y1_e, color="b")

    plt.plot(xupper, yupper, color=colour1, alpha=0.4)
    plt.plot(xlower, ylower, color=colour1, alpha=0.4)
    plt.vlines(x=x2, ymin=y3, ymax=y2, color=colour1, alpha=0.4)
    plt.vlines(x=x1, ymin=y4, ymax=y1, color=colour1, alpha=0.4)
    plt.fill_between(xlower, ylower, yupper, color=colour2, alpha=0.1)
    return

Plot_Graph(crh=4, cth=3,bh=12,Sweep_LE_HT=10*np.pi/180,Sweep_half_HT=10*np.pi/180,real=0,colour1="#11aa00",colour2="#11aa0055", trimtab=0,elevator=0,x_inboard_e=-1,x_outboard_e=-6,x_inboard_t=-2,x_outboard_t=-5,ce_c=0.3,cf_c=0.1, horn=0,w_horn=0.4,h_horn=0.8)
# Plot_Graph(crh=5, cth=3,bh=14,Sweep_LE_HT=10*np.pi/180,Sweep_half_HT=0,real=1, colour1="#728FCE", colour2="#98AFC7", trimtab=1,elevator=1,x_inboard_e=0,x_outboard_e=0)
# Plot_Graph(crh=7, cth=4,bh=14,Sweep_LE_HT=10*np.pi/180,Sweep_half_HT=0,real=1, colour1="#AFDCEC", colour2="#DCE9FA",trimtab=1,elevator=1,x_inboard_e=0,x_outboard_e=0)
plt.show()