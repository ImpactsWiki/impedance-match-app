# import python libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import panel as pn

# Required Impedance Match Calculation classes and functions for this app
import IM_module as IM 

def IM_app(webappbool=False):
    """ Shock Impedance Matching Tool and code. 
        Usage: IM_app(matdata,imat,webappbool=False)
        Inputs: materials parameters DataFrame and database indices objects.
                Optional boolean to use local copy of IHED in subdirectory database-ihed.
                If material not local, will fetch from IHED web site and save a local copy.
        Output: panel app (https://panel.holoviz.org/index.html)
        Requires IM_module loaded as IM. Variables in MKS.
        v1.0.0 - November 12, 2022 - S.T.Stewart
        v1.0.1 - 11/13/22 STS code cleanup and documentation
        v1.0.2 - 11/14/22 STS framing the app
        v1.0.3 - 11/15/22 STS auto-resize the plot; added panel to add a new material
    """
    #========================================================
    # WIDGET SETTINGS & CUSTOMIZATIONS
    # set default plot limits
    upmaxfactor = 1.2 # maximum x value is this factor * impact velocity
    pmaxfactor = 1.5 # maximum y value is this factor * first impact pressure
    # IM Hugoniot/release/reshock arrays
    # big factor in case mix of low and high impedance materials
    impvel_factor = 2 # make up arrays up to impvel_factor*vel 
    # widget splash image and default save file name
    default_image_filename = 'Impact-solution.pdf'
    #
    plt.rcParams["figure.figsize"] = (9,6)
    plt.rcParams["figure.dpi"] = 100
    plt.rc('font', size=10)
    menuwidth = 300 # pixel width for the left menu column
    #========================================================
    # read in materials database
    # needed by the add material function
    global matdata
    matdata, imat = IM.ReadMaterials(matfilename='materials-data.csv') # materials-data.csv is the default file name; can replace with your own file name

    # w* variables are widgets

    ab_materials = list(matdata.loc[:,'Material'].values)
    ab_materials.insert(0, 'Choose material')   # empty material at the top of the list
    bc_materials = list(matdata.loc[:,'Material'].values)
    bc_materials.insert(0, 'Choose material / No material')   # empty material at the top of the list

    wmat1=pn.widgets.Select(
        name='Material 1',
        options=ab_materials,
    )
    
    wmat2=pn.widgets.Select(
        name='Material 2',
        options=ab_materials,
    )
    
    wmat3=pn.widgets.Select(
        name='Material 3',
        options=bc_materials,
    )
    
    wmat4=pn.widgets.Select(
        name='Material 4',
        options=bc_materials,
    )

    wpmax = pn.widgets.FloatInput(name='Set plot max P (GPa), 0 to autoscale', value=0,step=5, start=0,width=int(menuwidth/2))

    wvel = pn.widgets.FloatInput(name='Impact Velocity (km/s)', value=0, step=1e-1, start=0, end=100)

    wshowdata = pn.widgets.Checkbox(
        value=False,
        name='Show IHED data (takes a minute)',
    )
    wuselocaldata = pn.widgets.Checkbox(
        value=False,
        name='Use/save local IHED data',
    )
    wusehugoniot = pn.widgets.Checkbox(
        value=True,
        name='Use Hugoniot for release and reshock',
    )
    wusehugoniotnote = pn.widgets.StaticText(value='(Mie-Grueneisen model can be unstable)')

    winfo = pn.widgets.StaticText(value='')

    # group to save PDF file when running in a Jupyter notebook ; does not work in Heroku web app
    winstruct = pn.widgets.StaticText(value='Enter file name with extension (.pdf, .png)')
    winstruct2 = pn.widgets.StaticText(value='Saves to local dir (online in CoLab)')
    wbutton = pn.widgets.Button(name='Save Plot', button_type='primary',width=100)
    wfilename = pn.widgets.TextInput(value=default_image_filename)
    wsaveimage = pn.Column(winstruct,winstruct2,wfilename,wbutton,width=menuwidth)
    #end group

    if webappbool:
        wcolumn1 = pn.Column('## Shock Impedance Match Tool\nSelect 2 or more materials and impact velocity.<br>Changing impact velocity updates the plot.', wmat1, wmat2, wmat3, wmat4, wvel, wshowdata, wuselocaldata, wusehugoniot, wusehugoniotnote, wpmax, winfo, width=menuwidth)
    else:
        wcolumn1 = pn.Column('## Shock Impedance Match Tool\nSelect 2 or more materials and impact velocity<br>Changing impact velocity updates the plot.', wmat1, wmat2, wmat3, wmat4, wvel, wshowdata, wuselocaldata, wusehugoniot, wusehugoniotnote, wpmax, wsaveimage, winfo, width=menuwidth)

    @pn.depends(vel=wvel,usehugoniot=wusehugoniot,showdata=wshowdata)
    def plot(vel,usehugoniot,showdata): 
        # usehugoniot and showdata are function parameters to trigger redraw of plot
        # the code below accesses the widget values directly
        vel = vel*1.e3 # put impact velocity in m/s
        userinfostr=''
        if vel > 0:
            #print('INCLUDES Mie-Gruneisen Release or Reshock for Material 2.')
            #print('Requires g0 and q in materials database.')
            # initialize strings for output information
            string1=''
            string2=''
            string3=''
            string4=''
            # initialize index for material layer
            id1=[]
            id2=[]
            id3=[]
            id4=[]

            # INITIALIZE MATERIALS FROM DROPDOWN MENUS
            # initialize common up array for all the curves used in IM solution
            up = np.arange(0,2001)/2000.*vel*impvel_factor # m/s

            # find index and initialize Hugoniot for each material layer (2 to 3 materials needed)
            id1 = np.where(matdata.loc[:,'Material'].values == wmat1.value)[0]
            if len(id1)>0: # material 1 has been set
                mat1 = IM.Material()
                mat1.DefineParamsID(wmat1.value,matdata,imat)
                mat1.MakeHugoniot(up)
                if mat1.s2 != 0 and mat1.d == 0:
                    userinfostr=userinfostr+'<br> WARNING: '+mat1.name+' Hugoniot is quadratic.'
                if (wshowdata.value) and (mat1.ihed.id > -1):
                    mat1.GetIHED(uselocalbool=wuselocaldata.value)

            id2 = np.where(matdata.loc[:,'Material'].values == wmat2.value)[0]
            if len(id2)>0:
                mat2 = IM.Material()
                mat2.DefineParamsID(wmat2.value,matdata,imat)
                mat2.MakeHugoniot(up)
                if mat2.s2 != 0 and mat2.d == 0:
                    userinfostr=userinfostr+'<br> WARNING: '+mat2.name+' Hugoniot is quadratic.'
                if wshowdata.value and (mat2.ihed.id > -1):
                    mat2.GetIHED(uselocalbool=wuselocaldata.value)

            id3 = np.where(matdata.loc[:,'Material'].values == wmat3.value)[0]
            if len(id3)>0:
                mat3 = IM.Material()
                mat3.DefineParamsID(wmat3.value,matdata,imat)
                mat3.MakeHugoniot(up)
                if mat3.s2 != 0 and mat3.d == 0:
                    userinfostr=userinfostr+'<br> WARNING: '+mat3.name+' Hugoniot is quadratic.'
                if wshowdata.value and (mat3.ihed.id > -1):
                    mat3.GetIHED(uselocalbool=wuselocaldata.value)

            id4 = np.where(matdata.loc[:,'Material'].values == wmat4.value)[0]
            if len(id4)>0:
                mat4 = IM.Material()
                mat4.DefineParamsID(wmat4.value,matdata,imat)
                mat4.MakeHugoniot(up)
                if mat4.s2 != 0 and mat4.d == 0:
                    userinfostr=userinfostr+'<br> WARNING: '+mat4.name+' Hugoniot is quadratic.'
                if wshowdata.value and (mat4.ihed.id > -1):
                    mat4.GetIHED(uselocalbool=wuselocaldata.value)

            # check that at least 2 materials have been defined
            if len(id1)==0 or len(id2)==0:
                winfo.value = 'No plot. Need materials 1 and 2.'
                return

            fig = plt.figure()
            def on_button_clicked(b):
                fig.savefig(wfilename.value,bbox_inches='tight',dpi=300)
            wbutton.on_click(on_button_clicked)

            # SOLVE FOR FIRST IM STATE: mat1--> mat2 at vel
            if len(id1) >0 and len(id2)>0:
                # solve for first impedance match state mat1--> mat2 at vel
                # use general function for maximum flexibility in curve construction: Intersection(x1,y1,x2,y2):
                res = IM.Intersection(up,mat2.hug.parr,vel-up,mat1.hug.parr)
                #print(res)
                mat1.im1.up=res[0][0] # m/s
                mat1.im1.p=res[1][0] # Pa
                mat1.im1.v=1./(np.interp(res[0],up,1./mat1.hug.varr)[0]) # assumes f is monotonic and increasing
                mat1.im1.e=np.interp(res[0],up,mat1.hug.earr)[0] # assumes f is monotonic and increasing
                mat2.im1.up=res[0][0]
                mat2.im1.p=res[1][0]
                mat2.im1.v=1./(np.interp(res[0],up,1./mat2.hug.varr)[0]) # assumes f is monotonic and increasing
                mat2.im1.e=np.interp(res[0],up,mat2.hug.earr)[0] # assumes f is monotonic and increasing
                #print(mat2.im1.p,mat1.im1.up)
                string1 = wmat1.value+' impacts '+wmat2.value+' at '+str(vel/1.e3)+' km/s\nImp. Match Mat2: Up='+IM.ClStr(mat2.im1.up/1.e3)+' (km/s) P='+IM.ClStr(mat2.im1.p/1.e9)+' (GPa)'
                plt.plot((vel-up)/1.e3,mat1.hug.parr/1.e9,label='Mat1 '+wmat1.value+' Hug.')
                plt.plot(up/1.e3,mat2.hug.parr/1.e9,label='Mat2 '+wmat2.value+' Hug.')
                plt.plot(mat2.im1.up/1.e3,mat2.im1.p/1.e9,'o',color='black',label='Mat2 Shock')


                #plt.plot(2*up_match_fix-material2.uparr/1.e3,material2.psarr/1.e9,label=mat2.value+' MG-release')

                if (wshowdata.value) and (mat1.ihed.id > -1):
                    indm1=np.where(mat1.ihed.marr == 1.)
                    if mat1.name == 'Ice':
                        indm1 = np.where(mat1.ihed.marr == 1.093)[0] # IHED used liquid water density
                    plt.scatter(vel/1.e3-mat1.ihed.uparr[indm1]/1.e3,mat1.ihed.parr[indm1]/1.e9,label=mat1.ihed.matname)
                    indm1=np.where(mat2.ihed.marr == 1.)
                if (wshowdata.value) and (mat2.ihed.id > -1):               
                    if mat2.name == 'Ice':
                        indm1 = np.where(mat2.ihed.marr == 1.093)[0] # IHED used liquid water density
                    plt.scatter(mat2.ihed.uparr[indm1]/1.e3,mat2.ihed.parr[indm1]/1.e9,label=mat2.ihed.matname)

                if len(id3)>0:
                    # optional third material included
                    # plot mat3 principal Hugoniot
                    plt.plot(up/1.e3,mat3.hug.parr/1.e9,label='Mat3 '+wmat3.value+' Hug.')
                    if wshowdata.value and (mat3.ihed.id > -1):
                        indm1=np.where(mat3.ihed.marr == 1.)
                        if mat3.name == 'Ice':
                            indm1 = np.where(mat3.ihed.marr == 1.093)[0] # IHED used liquid water density
                        plt.scatter(mat3.ihed.uparr[indm1]/1.e3,mat3.ihed.parr[indm1]/1.e9,label=mat3.ihed.matname)

                    # Calculate driver release
                    mat2.MakeMGIsentrope(mat2.im1,useHugoniotbool=wusehugoniot.value)
                    #print('Mat2 Isentrope from P=',mat2.im1.p/1.e9,' GPa up=',mat2.im1.up/1.e3,' km/s')
                    up_offset = np.interp(mat2.im1.p,mat2.isen.parr,mat2.isen.uparr)
                    #plt.plot(mat2.isen.uparr/1.e3,mat2.isen.parr/1.e9,label='Mat2 Isentrope 2')
                    res = IM.Intersection(up,mat3.hug.parr,2*up_offset-mat2.isen.uparr,mat2.isen.parr)
                    #print('release or reshock: ',res[0],up_offset)
                    if (up_offset > 0) and (len(res) > 0):
                        # Mie-Grueneisen model mat3 should be OK                        
                        if res[0][0] > up_offset:
                            # release
                            res2 = IM.Intersection(up,mat3.hug.parr,(2*up_offset-mat2.isen.uparr),mat2.isen.parr)
                            #print('Release res2=',res2)
                            mat2.im2.up=res2[0][0] # m/s
                            mat2.im2.p=res2[1][0] # Pa
                            mat2.im2.v=1./(np.interp(res2[0],up,1./mat2.isen.varr)[0]) # assumes f is monotonic and increasing
                            mat2.im2.e=np.interp(res2[0],up,mat2.isen.earr)[0] # assumes f is monotonic and increasing
                            mat3.im1.up=res2[0][0]
                            mat3.im1.p=res2[1][0]
                            mat3.im1.v=1./(np.interp(res2[0],up,1./mat3.hug.varr)[0]) # assumes f is monotonic and increasing
                            mat3.im1.e=np.interp(res2[0],up,mat3.hug.earr)[0] # assumes f is monotonic and increasing
                            plt.plot((2*up_offset-mat2.isen.uparr)/1.e3,mat2.isen.parr/1.e9,label='Mat2 Isentrope')
                            plt.scatter(mat3.im1.up/1.e3,mat3.im1.p/1.e9,label='IM mat3')
                            string3B = ' RELEASE'
                        else:
                            # reshock
                            mat2.MakeReshockHug(mat2.im1,useHugoniotbool=wusehugoniot.value)
                            res2 = IM.Intersection(up,mat3.hug.parr,mat2.reshock.uparr,mat2.reshock.parr)
                            #print('Reshock res2=',res2[0])
                            mat2.im2.up=res2[0][0] # m/s
                            mat2.im2.p=res2[1][0] # Pa
                            mat2.im2.v=1./(np.interp(res2[0],up,1./mat2.reshock.varr)[0]) # assumes f is monotonic and increasing
                            mat2.im2.e=np.interp(res2[0],up,mat2.reshock.earr)[0] # assumes f is monotonic and increasing
                            mat3.im1.up=res2[0][0]
                            mat3.im1.p=res2[1][0]
                            mat3.im1.v=1./(np.interp(res2[0],up,1./mat3.hug.varr)[0]) # assumes f is monotonic and increasing
                            mat3.im1.e=np.interp(res2[0],up,mat3.hug.earr)[0] # assumes f is monotonic and increasing
                            plt.plot(mat2.reshock.uparr/1.e3,mat2.reshock.parr/1.e9,label='Mat2 reshock')
                            plt.scatter(mat3.im1.up/1.e3,mat3.im1.p/1.e9,label='IM mat3')
                            string3B = ' RESHOCK'
                        if wusehugoniot.value:
                            string3 = '\nImp. Match Mat3 (HUG): Up='+IM.ClStr(mat3.im1.up/1.e3)+' (km/s) P='+IM.ClStr(mat3.im1.p/1.e9)+' (GPa)'+string3B
                        else:
                            string3 = '\nImp. Match Mat3 (MG): Up='+IM.ClStr(mat3.im1.up/1.e3)+' (km/s) P='+IM.ClStr(mat3.im1.p/1.e9)+' (GPa)'+string3B
                    else:
                        # something wrong with mat3 Mie Grueneisn model
                        userinfostr = userinfostr + '<br> ERROR WITH MIE-GRUENEISEN MODEL FOR MAT 3. Use Hugoniot for reshock/release.'
                        id4=[] # skip 4th material if there is a problem with mat3
                        winfo.value = 'ERROR WITH MIE-GRUENEISEN MODEL FOR MAT 3. Use Hugoniot for reshock/release.'
                        plt.close(fig)
                        return fig
                if len(id4)>0:
                    # optional fourth material included
                    # plot mat4 principal Hugoniot
                    plt.plot(up/1.e3,mat4.hug.parr/1.e9,label='Mat4 '+wmat4.value+' Hug.')
                    if (wshowdata.value) and (mat4.ihed.id > -1):               
                        indm1=np.where(mat4.ihed.marr == 1.)
                        if mat4.name == 'Ice':
                            indm1 = np.where(mat4.ihed.marr == 1.093)[0] # IHED used liquid water density
                        plt.scatter(mat4.ihed.uparr[indm1]/1.e3,mat4.ihed.parr[indm1]/1.e9,label=mat4.ihed.matname)

                    # Calculate driver release
                    mat3.MakeMGIsentrope(mat3.im1,useHugoniotbool=wusehugoniot.value)
                    #print('Mat2 Isentrope from P=',mat2.im1.p/1.e9,' GPa up=',mat2.im1.up/1.e3,' km/s')
                    up_offset = np.interp(mat3.im1.p,mat3.isen.parr,mat3.isen.uparr)
                    #plt.plot(mat2.isen.uparr/1.e3,mat2.isen.parr/1.e9,label='Mat2 Isentrope 2')

                    res = IM.Intersection(up,mat4.hug.parr,2*up_offset-mat3.isen.uparr,mat3.isen.parr)
                    #print('release or reshock: ',res[0],up_offset)
                    if (up_offset > 0) and (len(res) > 0):
                        if res[0][0] > up_offset:
                            # release
                            res3 = IM.Intersection(up,mat4.hug.parr,(2*up_offset-mat3.isen.uparr),mat3.isen.parr)
                            #print('Release res2=',res2)
                            mat3.im2.up=res3[0][0] # m/s
                            mat3.im2.p=res3[1][0] # Pa
                            mat3.im2.v=1./(np.interp(res3[0],up,1./mat3.isen.varr)[0]) # assumes f is monotonic and increasing
                            mat3.im2.e=np.interp(res3[0],up,mat3.isen.earr)[0] # assumes f is monotonic and increasing
                            mat4.im1.up=res3[0][0]
                            mat4.im1.p=res3[1][0]
                            mat4.im1.v=1./(np.interp(res3[0],up,1./mat4.hug.varr)[0]) # assumes f is monotonic and increasing
                            mat4.im1.e=np.interp(res3[0],up,mat4.hug.earr)[0] # assumes f is monotonic and increasing
                            plt.plot((2*up_offset-mat3.isen.uparr)/1.e3,mat3.isen.parr/1.e9,label='Mat3 Isentrope')
                            plt.scatter(mat4.im1.up/1.e3,mat4.im1.p/1.e9,label='IM mat4')
                            string4B = ' RELEASE'
                        else:
                            # reshock
                            mat3.MakeReshockHug(mat3.im1,useHugoniotbool=wusehugoniot.value)
                            res3 = IM.Intersection(up,mat4.hug.parr,mat3.reshock.uparr,mat3.reshock.parr)
                            #print('Reshock res3=',res3[0])
                            mat3.im2.up=res3[0][0] # m/s
                            mat3.im2.p=res3[1][0] # Pa
                            mat3.im2.v=1./(np.interp(res3[0],up,1./mat3.reshock.varr)[0]) # assumes f is monotonic and increasing
                            mat3.im2.e=np.interp(res3[0],up,mat3.reshock.earr)[0] # assumes f is monotonic and increasing
                            mat4.im1.up=res3[0][0]
                            mat4.im1.p=res3[1][0]
                            mat4.im1.v=1./(np.interp(res3[0],up,1./mat4.hug.varr)[0]) # assumes f is monotonic and increasing
                            mat4.im1.e=np.interp(res3[0],up,mat4.hug.earr)[0] # assumes f is monotonic and increasing
                            plt.plot(mat3.reshock.uparr/1.e3,mat3.reshock.parr/1.e9,label='Mat3 reshock')
                            plt.scatter(mat4.im1.up/1.e3,mat4.im1.p/1.e9,label='IM mat4')
                            string4B = ' RESHOCK'
                        if wusehugoniot.value:
                            string4 = '\nImp. Match Mat4 (HUG): Up='+IM.ClStr(mat4.im1.up/1.e3)+' (km/s) P='+IM.ClStr(mat4.im1.p/1.e9)+' (GPa)'+string4B
                        else:
                            string4 = '\nImp. Match Mat4 (MG): Up='+IM.ClStr(mat4.im1.up/1.e3)+' (km/s) P='+IM.ClStr(mat4.im1.p/1.e9)+' (GPa)'+string4B
                    else:
                        # something wrong with mat4 Mie Grueneisn model
                        userinfostr = userinfostr + '<br> ERROR WITH MIE-GRUENEISEN MODEL FOR MAT 4. Use Hugoniot for reshock/release.'

                plt.title(string1+string2+string3+string4)            
                plt.legend(bbox_to_anchor=(1,1), loc="upper left") #bbox_to_anchor=(1.05, 1)
                plt.xlabel('Particle Velocity (km/s)\nhttps://impactswiki.net/impact-tools-book/')
                plt.ylabel('Pressure (GPa)')
                if wpmax.value > 0: # don't use negative values
                    plt.ylim(0,wpmax.value)
                else:
                    plt.ylim(0,pmaxfactor*mat1.im1.p/1.e9)
                plt.xlim(0,upmaxfactor*vel/1.e3)
                plt.tight_layout()
                if wusehugoniot.value:
                    userinfostr = userinfostr + '<br> Using Hugoniot for reshock and release.'
                else:
                    userinfostr = userinfostr + '<br> Using Mie-Grueneisen model for reshock and release.'
                winfo.value = 'Updated plot, impact vel (km/s)='+IM.ClStr(vel/1.e3)+userinfostr
                plt.close(fig)
            return fig
        else:
            winfo.value = 'No plot. Impact velocity <= 0'

    wplot=pn.panel(plot, sizing_mode='scale_width') # panel to display the impedance match plot

    wbottomtext = pn.widgets.StaticText(value='<b>Manual</b> <a href="https://impactswiki.net/impact-tools-book/">https://impactswiki.net/impact-tools-book/</a><br><b>Repo</b> <a href="https://github.com/ImpactsWiki/impedance-match-app">https://github.com/ImpactsWiki/impedance-match-app</a><br>If crashing or unresponsive, use Hugoniot for release and reshock. Common problems with large impedance mismatch materials and MG models.')
    
    # Widgets below the main IM Tool
    # display current matdata DataFrame
    wdf_widget = pn.widgets.Tabulator(matdata)
    # entry boxes for new material parameters
    waddmattext = pn.widgets.StaticText(value='Enter values for new material in MKS. <a href="http://www.ihed.ras.ru/rusbank/substsearch.php" target="_blank">Search the IHED Database to find the material number</a>; enter -1 if not available. The substance ID number is given in the web page address when plotting the Hugoniot or displaying the data as plain text.<p> 2, 3, or 4 non-zero parameters define the form of the Hugoniot:<br>2=Linear: Us = c0 + s1*up<br>3=Quadratic: Us = c0 + s1*up + s2*up^2<br>4=Mod. Universal Liquid: Us = c0 + s1*up - c*up*exp(-d*up)<p>Mie-Grueneisen parameter: g(v) = g0*(v/v0)^q')
    if webappbool:
        waddmattext.value = 'Changes to materials database will be lost upon refreshing this web app. <a href="https://impactswiki.net/impact-tools-book/" target="_blank">Use the Jupyter Notebook version</a> to access a personal database.<p>'+waddmattext.value
    wnewname = pn.widgets.TextInput(name='New Material Name',value='Enter Name')
    wrho0 = pn.widgets.FloatInput(name='Density [kg/m^3]', page_step_multiplier=1)
    wc0 = pn.widgets.FloatInput(name='c0 [m/s]', page_step_multiplier=1)
    ws1 = pn.widgets.FloatInput(name='s1 [-]', page_step_multiplier=1)
    ws2 = pn.widgets.FloatInput(name='s2 [s/m] or c [-]', page_step_multiplier=1)
    wd  = pn.widgets.FloatInput(name='d [s/m]', page_step_multiplier=1)
    wg0 = pn.widgets.FloatInput(name='g0 [-]', value=1., page_step_multiplier=1)
    wq = pn.widgets.FloatInput(name='q [-]', value=1., page_step_multiplier=1)
    wihed = pn.widgets.IntInput(name='IHED substance id number [int]', value=-1, page_step_multiplier=1)
    wdate = pn.widgets.TextInput(name='Date', value='YYYY-MM-DD')
    
    wnote = pn.widgets.TextInput(name='Note', value='Enter reference for material parameters.')
    waddmatbutton = pn.widgets.Button(name='Add material to database', button_type='primary')
    wnewparams = pn.Column(waddmattext,wnewname,wrho0,wc0,ws1,ws2,wd,wg0,wq,wihed,wdate,wnote,waddmatbutton,width=500)
    
    @pn.depends(wmat1,wmat2,wmat3,wmat4)
    def on_addmatbutton_clicked(event):
        print('CLICKED BUTTON')
        # somebody please tell me the better way to access these variables in a button function....
        global matdata
        # load an empty new material parameter DataFrame
        # Columns must match the original materials-data.csv file
        newmatdata, imat = IM.ReadMaterials(matfilename='materials-new.csv')
        newmatdata.iloc[:,imat.name] = wnewname.value
        newmatdata.iloc[:,imat.rho0] = wrho0.value
        newmatdata.iloc[:,imat.c0] = wc0.value
        newmatdata.iloc[:,imat.s1] = ws1.value
        newmatdata.iloc[:,imat.s2] = ws2.value
        newmatdata.iloc[:,imat.d] = wd.value
        newmatdata.iloc[:,imat.g0] = wg0.value
        newmatdata.iloc[:,imat.q] = wq.value
        newmatdata.iloc[:,imat.ihed] = wihed.value
        newmatdata.iloc[:,imat.date] = wdate.value
        newmatdata.iloc[:,imat.note] = wnote.value
        matdata = pd.concat([newmatdata,matdata],ignore_index=True)
        wdf_widget.value = matdata
        ab_materials = list(matdata.loc[:,'Material'].values)
        ab_materials.insert(0, 'Choose material')   # empty material at the top of the list
        wmat1.options=ab_materials
        wmat2.options=ab_materials
        bc_materials = list(matdata.loc[:,'Material'].values)
        bc_materials.insert(0, 'Choose material/nomaterial')   # empty material at the top of the list
        wmat3.options=bc_materials
        wmat4.options=bc_materials
        
    waddmatbutton.on_click(on_addmatbutton_clicked)#,matdata=matdata,ab=ab_materials,bc=bc_materials)
    
    # author into at bottom of app
    wauthortext = pn.widgets.StaticText(value='&#169; 2022 S. T. Stewart, Planetary Impacts Community Wiki')

    # collect the various parts of the web app
    wtop_pane = pn.pane.PNG('PetaviusLangrenus_Poupeau_3000.png',link_url="https://impacts.wiki",sizing_mode="scale_width")
    wmain_pane = pn.Row(wcolumn1,wplot,sizing_mode="scale_width")
    wmatdata_pane = pn.Card(wdf_widget, title="Materials Database", sizing_mode='scale_width', collapsed=True)
    waddmat_pane = pn.Card(wnewparams, title="Add Material", sizing_mode='scale_width', collapsed=True)
    
    wcombo_pane = pn.Column(wtop_pane,wmain_pane,wbottomtext,wmatdata_pane,waddmat_pane,pn.layout.Divider(),wauthortext,width=1200,sizing_mode="scale_width")
    
    return wcombo_pane

## this function creates a panel app that can be run as a standalone app in a web browser using
## IM_app(matdata,imat)
## use IM_app.servable() to run the widget alone in a web browser using command line: bokeh serve --show filename.ipynb
## see https://panel.holoviz.org/index.html
### END of IM_app.py ###
    