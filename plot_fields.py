#!/usr/bin/python

"""plot_fields.py
Danielle McDermott
June 22, 2018

Collection of subroutines designed for plotting 
2D position data and making figures with Gridspec. 
The data files were generated by code Vortex2D.
Takes as a 'main' function Fig1_400.py, 
so by default it will create Figure 1 of 'Stripes on Stripes' paper
written in May 2014

Some notes:
Currently working to make a more general version of 
ascii_multiplot.py

ascii plot and heat plot have a lot in common,
currently combined into a single subroutine

TODO: add a plot class
"""



#######################################################
#Import Python Library Modules
#######################################################
import matplotlib
matplotlib.use('Agg')   #for batch jobs!!!

import os, math, sys, csv, glob, numpy

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#contour plots
from matplotlib import cm
#3D plots
#from mpl_toolkits.mplot3d import axes3d

import data_importerDM as diDM

from matplotlib import ticker

from mpl_toolkits.axes_grid1 import AxesGrid

def letter_range(start, stop):
    for c in range(ord(start), ord(stop)):
        yield chr(c)
    return

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    SOURCE:
    http://stackoverflow.com/questions/
    /7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = numpy.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = numpy.hstack([
        numpy.linspace(0.0, midpoint, 128, endpoint=False), 
        numpy.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


##################################################################
#plot the velocity field of a system of format
#x y vx vy
##################################################################
def plot_velocity_field(ax,filename="coarse_grain_field.txt",
                        Sx=[0,60.0],Sy=[0,60.0],
                        plot_curl=0,
                        plot_coarse_field=0,
                        plot_cbar=0,
                        vmin=None,vmax=None,
                        shift=None):
    
    '''Read in a formatted file containing a velocity field
    and make a vector plot
    '''

    #ax is declared in main function
    ax.set_aspect('equal')
    
    #get the formatted data from c code calculation
    (x,y,vx,vy,curl,enstrophy) = diDM.get_data(filename,6,sep='\t',path1=os.getcwd())

    if plot_curl:
        #create regular x-y grid
        xi = numpy.linspace(min(x), max(x))
        yi = numpy.linspace(min(y), max(y))

        #create a mesh x-y grid from the numpy.linspace
        #these are necessary for the plot
        X, Y = numpy.meshgrid(xi, yi)
    
        #take the Heat (x,y,z) data and grid it on the xi, yi (not the mesh)
        Z = matplotlib.mlab.griddata(x, y, curl, xi, yi, interp='linear')

        #############################################
        # color map of curl of coarse grain field
        #############################################
        
        orig_cmap=cm.RdBu
        kwargs_cset={'offset':0}
        if shift:
            shifted_cmap = shiftedColorMap(orig_cmap, 
                                           midpoint=shift, 
                                           name='shifted')
            kwargs_cset['cmap']=shifted_cmap
        else:
            kwargs_cset['cmap']=orig_cmap

        if vmin:
            kwargs_cset['vmin']=vmin
        if vmax:
            kwargs_cset['vmax']=vmax
            
        cset = ax.contourf(X, Y, Z, zdir='curl',**kwargs_cset)
        #cbar = fig.colorbar(cset, ax=ax)

        '''
        if plot_cbar:
            cbar_ax=add_subplot_axes(ax,rect=[-0.5,0.9,2.0,0.05])
            cb = plt.colorbar(cset,ax=ax,cax=cbar_ax,orientation='horizontal') 
        cbar_ax.set_title(r"$\vec{\omega}=\nabla \times \vec{v}$",fontsize=18)
        '''

    if plot_coarse_field:
        ax.quiver(x,y,vx,vy,pivot='middle',headwidth=20,headlength=20)

    ax.set_xlim(*Sx)        
    ax.set_ylim(*Sy)



    return #cset

##################################################################
#plot the velocity field of particles, format
#id x y vx vy |v|
##################################################################
def particle_velocity_field(ax, filename="velocity_frame22000", 
                            Sx=[0,60.0], Sy=[0,60.0],
                            BW=0, plot_log=0, scale=1.0,
                            xshift=0.0, yshift=0.0, 
                            arrow_color='black',shrink=0,
                            abs_SX = 60.0):

    '''Read in a formatted file containing a velocity field
    and make a vector plot
    '''
    
    size_vectors=12 #currently hardcoded
    #ax is declared in main function
    ax.set_aspect('equal')

    #get the formatted data from c code calculation
    print("The filename is:", filename)
    print("")
    (i,type,x,y,vx,vy,radius,speed) = diDM.get_data(filename,8,path1=os.getcwd())    

    #############################################################
    #Here we are typically centering the particles in the figure
    #############################################################
    #beware, intentionally hardcoded, meant as an absolute size
    
    if xshift:
        for j in range(len(x)):
            x[j] = x[j]+xshift
            if x[j] > abs_SX:
                x[j] -= abs_SX
            elif x[j] < 0.0:
                x[j] += abs_SX

    if yshift:
        for j in range(len(y)):
            y[j] = y[j]+yshift

            if y[j] > abs_SX:
                y[j] -= abs_SX
            elif y[j] < 0.0:
                y[j] += abs_SX

    
    if shrink:
        center=0.5*Sx[-1]
        ax.set_xlim(center-5.0,center+7.0)
    else:
        ax.set_xlim(*Sx)        

    ax.set_ylim(*Sy)


    #use a list to track the arguments we will call with quiver
    quiver_args=[x,y,vx,vy]


    #use a dictionary to track the optional arguments
    quiver_kargs={'pivot':'middle',
                  'headwidth':size_vectors,'headlength':size_vectors}

    if scale != 1.0:
        scale_vectors=1
        quiver_kargs['scale']=scale

    if BW == 0:
        #speend the list to color the arrows
        quiver_args.append(speed)

        #use a logscale to color the arrows by speed
        if plot_log:
            vmin=1e-7
            vmax=1e-3
            norm1=matplotlib.colors.LogNorm(vmin,vmax)
            #append our optional arguments with color map information
            quiver_kargs['cmap']=cm.jet
            quiver_kargs['norm']=norm1

    #Black and white arrows
    else:

        if arrow_color!='black':
            quiver_kargs['color']=arrow_color
            
    #now call quiver with the args and optional args we have defined!
    ax.quiver(*quiver_args,**quiver_kargs) 

    ax.set_xlim(*Sx)        
    ax.set_ylim(*Sy)

    return

##################################################################
#plot the velocity field of particles, format
#id x y vx vy |v|
##################################################################
def enstrophy_field(ax,filename="coarse_grain_field.txt",
                    Sx=[0,60.0],Sy=[0,60.0],BW=0,plot_log=0,plot_cbar=1):
    
    '''Read in a formatted file containing a velocity field
    and make a vector plot
    you could totally combine this with the curl field
    and possibly even the heat map
    as a general subroutine
    '''

    #ax is declared in main function
    ax.set_aspect('equal')
    
    #get the formatted data from c code calculation
    (x,y,vx,vy,speed,enstrophy) = diDM.get_data(filename,6,sep='\t',path1=os.getcwd())

    ax.set_xlim(*Sx)        
    ax.set_ylim(*Sy)

    #create regular x-y grid
    xi = numpy.linspace(min(x), max(x))
    yi = numpy.linspace(min(y), max(y))

    #create a mesh x-y grid from the numpy.linspace
    #these are necessary for the plot
    X, Y = numpy.meshgrid(xi, yi)
    
    #take the Heat (x,y,z) data and grid it on the xi, yi (not the mesh)
    Z = matplotlib.mlab.griddata(x, y, enstrophy, xi, yi, interp='linear')

    #############################################
    # color map of curl of coarse grain field
    #############################################
    cset = ax.contourf(X, Y, Z, zdir='enstrophy', offset=0, cmap=cm.Purples)
    #cbar = fig.colorbar(cset, ax=ax)
    '''
    if plot_cbar:
        cbar_ax=add_subplot_axes(ax,rect=[-1.0,0.0,2.0,0.05])
        cb = plt.colorbar(cset,ax=ax,cax=cbar_ax,orientation='horizontal')
        cbar_ax.set_title(r"$\varepsilon(\vec{\omega})$",fontsize=18)

        #,ticks=tick_loc )
    '''


    return #cset


#---------------------------------------------------------------------------
#nominal main.c-------------------------------------------------------------
#---------------------------------------------------------------------------
if __name__ == "__main__":

    ###################################################
    #Get Parameters 
    #(1) single value of Nv, Fp, Np
    #2x2 figure
    ###################################################
    
    N = 1075
    Sy = [0,60.0]
    Sx = [0,60.0]
    verbose = 1

    if(verbose == 1):
        print("Verbose settings, change as you see fit")
        print("Size hardcoded to: ", Sx, Sy)
        print("Number of particles (not used)", N)
        
    #######################################################
    #define number of row/columns to make a gridded figure
    #######################################################
    rows = 2
    columns = 2

    #########################################################
    ####################parameters got!######################
    #########################################################

    #turn on/off the "X" "Y" and axis ticks
    labels = 1

    #########################################################
    #identify data files
    #########################################################
        
    field_file = "velocity_data/field_00009950"
    particle_file = "velocity_data/XV_data_t=00009950"

    ###########################
    #Define standard figure
    ###########################
    fig = plt.figure( figsize=(columns*3,rows*3) )
        
    ############################################################
    #created grid for figure (subplots connect them, IMPORTANT)
    ############################################################
    G = gridspec.GridSpec(rows, columns, wspace=0.1, hspace=0.1)
        
    #######################################
    #populate a list for annotations
    #######################################
    letter = []
    for x in letter_range('a', 'z'):
        letter.append( x )  
    
    #########################################################
    #add letter annotation at (xt,yt)
    #########################################################
    xt = Sx[0]+1.5
    #y -position of text depends on whether or not label exists
    
    if labels:
        yt = Sy[1]-5.0        
    else:
        yt = Sy[1]-4.0

    #######################################
    #Iterate through grid adding subplots
    #the function called is different every time
    #######################################
    
    i = 0
    
    for a in range(rows):
        for b in range(columns):
            
            if verbose:
                print("Plotting grid position: a=",a,"b=",b)

            #grab a grid pointer
            plot_position = G[a,b]    

            #make a subplot at that location
            ax = fig.add_subplot(plot_position)        
            
            #grab file from the declared list
            file = field_file

            ##########################################
            #Case/Switch (ish) to plot (a-d)
            ##########################################
            if i == 0:    #(a) GR plot
                particle_velocity_field(ax,filename=particle_file,Sx=Sx,Sy=Sy)
                
                if verbose:
                    print("plotting velocity vectors in panel (a) - colored by velocity")

            elif i == 1:  #(b) plot the curl of the velocity field
                plot_velocity_field(ax,filename=file,plot_curl=1,Sx=Sx,Sy=Sy)

                if verbose:
                    print("plotting curl of field in panel (b)")

            elif i == 2:  #(c) plot the arrows, sized by magnitude
                plot_velocity_field(ax,filename=file,plot_coarse_field=1,Sx=Sx,Sy=Sy)
                
                if verbose:
                    print("plotting goofy vectors in (c) - sized by magnitude")
                

            elif i == 3:  #(d) plot the enstrophy
                enstrophy_field(ax,filename=file,plot_cbar=1,Sx=Sx,Sy=Sy)
                
                if verbose:
                    print("plotting enstrophy (vortical kinetic energy) (d)")
                

            else:
                print("some error in grid, plots, etc")
                exit(-1)

            ###########################################
            ##Add subplot letter (a), etc
            ##text position is hardwired, needs work
            ###########################################

            xt0 = xt
            yt0 = yt

            ax.text(xt0,yt0,"("+letter[i]+")",backgroundcolor='white',size=14)
            print(i)
            #increase i for every new a,b instance
            i += 1

            ###########################################
            #for loop ends here
            ###########################################

    ########################################
    #Tweak the figure
    ########################################
        
    if labels:
        G.tight_layout(fig, rect=[0, 0, 1, 1], h_pad=0.2, w_pad=0.2,pad=0.2)     
    else:
        p_val = -1.4
        G.tight_layout(fig, rect=[0, 0, 1, 1], h_pad=p_val,w_pad=p_val,pad=0.0)
    ###################################################
    #Save the Figure
    ###################################################
    out_file="fields.png"
    plt.savefig(out_file) #another cmd line opportunity
        
    sys.exit()
