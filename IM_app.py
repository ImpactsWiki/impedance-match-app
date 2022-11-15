# import python libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import panel as pn

# Required Impedance Match Calculation classes and functions for this app
import IM_module as IM 

def IM_app(matdata,imat,webappbool=False):
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

    ab_materials = list(matdata.loc[:,'Material'].values)
    ab_materials.insert(0, 'Choose material')   # empty material at the top of the list
    bc_materials = list(matdata.loc[:,'Material'].values)
    bc_materials.insert(0, 'Choose material/nomaterial')   # empty material at the top of the list

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

    wpmax = pn.widgets.FloatInput(name='Set plot Pmax (GPa), 0 to autoscale', value=0,step=5, start=0)

    wvel = pn.widgets.FloatInput(name='Impact Velocity (km/s) - change to update plot', value=0, step=1e-1, start=0, end=100)

    wshowdata = pn.widgets.Checkbox(
        value=False,
        name='Show IHED data (takes a minute)',
    )
    wuselocaldata = pn.widgets.Checkbox(
        value=False,
        name='Use/save local IHED data',
    )
    wusehugoniot = pn.widgets.Checkbox(
        value=False,
        name='Use Hugoniot for release and reshock',
    )

    winfo = pn.widgets.StaticText(value='')

    # group to save PDF file when running in a Jupyter notebook ; does not work in Heroku web app
    winstruct = pn.widgets.StaticText(value='Enter image file name with extension (.pdf, .png)')
    wbutton = pn.widgets.Button(name='Save Plot', button_type='primary',width=100)
    wfilename = pn.widgets.TextInput(value=default_image_filename)
    wsaveimage = pn.Column(winstruct,wfilename,wbutton,width=menuwidth)
    #end group

    if webappbool:
        column1 = pn.Column('## Shock Impedance Match Tool\nSelect 2 or more materials and impact velocity', wmat1, wmat2, wmat3, wmat4, wshowdata, wuselocaldata, wusehugoniot, wvel, wpmax, winfo, width=menuwidth)#, background='WhiteSmoke')
    else:
        column1 = pn.Column('## Shock Impedance Match Tool\nSelect 2 or more materials and impact velocity', wmat1, wmat2, wmat3, wmat4, wshowdata, wuselocaldata, wusehugoniot, wvel, wpmax, wsaveimage, winfo, width=menuwidth)#, background='WhiteSmoke')


    def plot(wvel):
        vel = wvel*1.e3 # put impact velocity in m/s
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
                if wshowdata.value:
                    mat1.GetIHED(uselocalbool=wuselocaldata.value)

            id2 = np.where(matdata.loc[:,'Material'].values == wmat2.value)[0]
            if len(id2)>0:
                mat2 = IM.Material()
                mat2.DefineParamsID(wmat2.value,matdata,imat)
                mat2.MakeHugoniot(up)
                if mat2.s2 != 0 and mat2.d == 0:
                    userinfostr=userinfostr+'<br> WARNING: '+mat2.name+' Hugoniot is quadratic.'
                if wshowdata.value:
                    mat2.GetIHED(uselocalbool=wuselocaldata.value)

            id3 = np.where(matdata.loc[:,'Material'].values == wmat3.value)[0]
            if len(id3)>0:
                mat3 = IM.Material()
                mat3.DefineParamsID(wmat3.value,matdata,imat)
                mat3.MakeHugoniot(up)
                if mat3.s2 != 0 and mat3.d == 0:
                    userinfostr=userinfostr+'<br> WARNING: '+mat3.name+' Hugoniot is quadratic.'
                if wshowdata.value:
                    mat3.GetIHED(uselocalbool=wuselocaldata.value)

            id4 = np.where(matdata.loc[:,'Material'].values == wmat4.value)[0]
            if len(id4)>0:
                mat4 = IM.Material()
                mat4.DefineParamsID(wmat4.value,matdata,imat)
                mat4.MakeHugoniot(up)
                if mat4.s2 != 0 and mat4.d == 0:
                    userinfostr=userinfostr+'<br> WARNING: '+mat4.name+' Hugoniot is quadratic.'
                if wshowdata.value:
                    mat4.GetIHED(uselocalbool=wuselocaldata.value)

            # check that at least 2 materials have been defined
            if len(id1)==0 or len(id2)==0:
                winfo.value = 'No plot. Need materials 1 and 2.'
                return

            fig = plt.figure()
            #def on_button_clicked(b):
            #    fig.savefig(wfilename.value,bbox_inches='tight',dpi=300)
            #wbutton.on_click(on_button_clicked)

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

                if wshowdata.value:
                    indm1=np.where(mat1.ihed.marr == 1.)
                    if mat1.name == 'Ice':
                        indm1 = np.where(mat1.ihed.marr == 1.093)[0] # IHED used liquid water density
                    plt.scatter(vel/1.e3-mat1.ihed.uparr[indm1]/1.e3,mat1.ihed.parr[indm1]/1.e9,label=mat1.ihed.matname)
                    indm1=np.where(mat2.ihed.marr == 1.)
                    if mat2.name == 'Ice':
                        indm1 = np.where(mat2.ihed.marr == 1.093)[0] # IHED used liquid water density
                    plt.scatter(mat2.ihed.uparr[indm1]/1.e3,mat2.ihed.parr[indm1]/1.e9,label=mat2.ihed.matname)

                if len(id3)>0:
                    # optional third material included
                    # plot mat3 principal Hugoniot
                    plt.plot(up/1.e3,mat3.hug.parr/1.e9,label='Mat3 '+wmat3.value+' Hug.')
                    if wshowdata.value:
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

                if len(id4)>0:
                    # optional fourth material included
                    # plot mat4 principal Hugoniot
                    plt.plot(up/1.e3,mat4.hug.parr/1.e9,label='Mat4 '+wmat4.value+' Hug.')
                    if wshowdata.value:
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
                winfo.value = 'Updated plot, impact vel (km/s)='+IM.ClStr(wvel)+userinfostr
                plt.close(fig)
            return fig
        else:
            winfo.value = 'No plot. Impact velocity <= 0'


    out=pn.bind(plot,wvel=wvel) # link the velocity widget to the plot function

    wbottomtext = pn.widgets.StaticText(value='<b>Manual</b> <a href="https://impactswiki.net/impact-tools-book/">https://impactswiki.net/impact-tools-book/</a><br><b>Repo</b> <a href="https://github.com/ImpactsWiki/impedance-match-app">https://github.com/ImpactsWiki/impedance-match-app</a><br>If crashing or unresponsive, use Hugoniot for release and reshock.')
    df_widget = pn.widgets.Tabulator(matdata)
    wtemptext = pn.widgets.StaticText(value='in the queue....')
    wauthortext = pn.widgets.StaticText(value='S. T. Stewart, 2022')

    top_pane = pn.pane.PNG('PetaviusLangrenus_Poupeau_3000.png',link_url="https://impacts.wiki",sizing_mode="scale_width")
    main_pane = pn.Row(column1,out,sizing_mode="scale_width")
    matdata_pane = pn.Card(df_widget, title="Materials Database", sizing_mode='scale_width', collapsed=True)
    addmat_pane = pn.Card(wtemptext, title="Add Material", sizing_mode='scale_width', collapsed=True)
    
    combo_pane = pn.Column(top_pane,main_pane,wbottomtext,matdata_pane,addmat_pane,pn.layout.Divider(),wauthortext,width=1200,sizing_mode="scale_width")
    
    return combo_pane

## this function creates a panel app that can be run as a standalone app in a web browser using
## IM_app(matdata,imat)
## use IM_app.servable() to run the widget alone in a web browser using command line: bokeh serve --show filename.ipynb
## see https://panel.holoviz.org/index.html
### END of IM_app.py ###
    