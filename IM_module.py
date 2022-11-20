# IM_module.py
# Define Material classes and functions for Impedance Match interactive
# Everything in MKS units
# S.T.Stewart
#
import numpy as np
import urllib.request   # part of python standard library; used to extract data from IHED database
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
#
__version__ = '1.1.1' # November 19, 2022 Debug MG release E,V calc
#__version__ = '1.1.0' # November 19, 2022 Added IM_match function; lots of debugging and error messaging.
#__version__ = '1.0.4' # November 14, 2022 Debugging.
#__version__ = '1.0.3' # November 13, 2022 Fixed symmetric impact bug; code cleanup and documentation.
#__version__ = '1.0.2' # November 12, 2022 Added Universal Liquid Hugoniot
#__version__ = '1.0.1' # November 12, 2022 Added more graceful failure mores
#
def UniversalHugoniot(up,a,b,c,d):
    """ Function returns modified universal liquid Hugoniot.
        Usage: UniversalHugoniot(up,a,b,c,d):
        Us = (a + b * up - c * up *np.exp(-d*up))
    """
    return (a + b * up - c * up *np.exp(-d*up)) # shock vel in m/s
    
class MaterialIndices:
    """Class for book-keeping the column indices in the material database file."""
    def __init__(self):
        self.name=''
        self.rho0=0
        self.c0=0
        self.s1=0
        self.s2=0
        self.d =0
        self.g0=0
        self.q=0
        self.uplow=0
        self.uphigh=0
        self.ihed=0
        self.date=0
        self.note=0
class EOS_Line:
    """Class for EOS line (Hugoniot, isentrope, reshock Hugoniot)"""
    def __init__(self):
        self.v0 = 0.
        self.p0 = 0.
        self.e0 = 0.
        self.t0 = 0.
        self.up0 = 0.
        self.parr = 0. # pressure Pa
        self.varr = 0. # specific volume m3/kg
        self.earr = 0. # specific energy J/kg
        self.uparr = 0. # particle velocity m/s
        self.uparr2 = 0. # particle velocity m/s FOR DEBUGGING
        self.usarr = 0. # shock velocity m/s
        self.tarr = 0. # temperature K
        self.sarr = 0. # specific entropy J/K/kg
        self.garr = 0. # Mie Grueneisen parameter [-]
class EOS_Point:
    """Class for EOS point."""
    def __init__(self):
        self.p   = 0. # Pa
        self.v   = 0. # m3/kg
        self.e   = 0. # J/kg
        self.up  = 0. # m/s
    def __str__(self):
        # print friendly output for EOS point to display to the user
        output_str = 'P = '+ClStr(self.p/1.e9)+' (GPa), u<sub>p</sub> = '+ClStr(self.up/1.e3)+' (km/s)'
        output_str = output_str + ', E = '+ClStr(self.e/1.e6)+' (MJ/kg)'
        if self.v > 0:
            output_str = output_str +', &#961; = '+ClStr(1/self.v/1.e3)+ ' (g/cm<sup>3</sup>)'
        else:
            output_str = output_str +', V = '+str(self.v)+ ' (m<sup>3</sup>/kg)'
        return output_str

class IHED:
    """Class for holding IHED or user-supplied shock wave data."""
    def __init__(self):
        self.id       = -1 # integer index for IHED material number
        self.rho0     = 0.
        self.uparr    = 0.
        self.usarr    = 0.
        self.parr     = 0.
        self.earr     = 0.
        self.varr     = 0.
        self.marr     = 0. # ratio to full density
        self.c0       = 0.
        self.s1       = 0.
        self.s2       = 0.
        self.d        = 0.
        self.matname  = ''
class Material:
    """Material class for impedance matching calculations."""
    def __init__(self):
        """Initialize material class for impedance matching calculations."""
        self.name = ''
        self.note = ''
        self.rho0 = 0. # initial density kg/m3
        self.v0   = 0. # initial specific volume = 1/rho0
        self.c0   = 0. # linear eos intercept m/s
        self.s1   = 0. # linear eos slope
        self.s2   = 0. # quadratic eos slope
        self.g0   = 0. # Mie Grueneisen parameter [-]
        self.q    = 0. # Mie Grueneisen parameter exponent [-] g0/v0=(g/v)^q
        self.p0   = 0. # initial pressure [Pa]
        self.e0   = 0. # initial specific internal energy [J/kg]
        self.t0   = 0. # initial temperature [K] 
        self.hug     = EOS_Line() # Principal Hugoniot
        self.im1     = EOS_Point() # impedance match state 1 
        self.im2     = EOS_Point() # impedance match state 2 
        self.isen    = EOS_Line() # Release isentrope
        self.reshock = EOS_Line() # Reshock Hugoniot
        self.ihed    = IHED() # structure to store IHED data for this material
        self.ihed2   = IHED() # structure to store a second group of IHED or user data for this material
    def DefineParams(self,name,rho0,c0,s1,s2,d,g0,q,ihednum,note):
        """ Initialize material parameters with manual parameter input.
            Usage: DefineParams(self,name,rho0,c0,s1,s2,d,g0,q,ihednum,note):
        """
        self.name = name # material name string
        self.rho0 = rho0 # kg/m3
        self.v0   = 1/rho0 # m3/kg
        self.c0   = c0 # m/s
        self.s1   = s1 # [-]
        self.s2   = s2 # [1/(m/s)] if quadratic or [-] if universal liquid
        self.d    = d  # [1/(m/s)]
        self.g0   = g0 # [-]
        self.q    = q # [-]        
        self.ihed.id  = int(ihednum) # [integer] -1 indicates no data in IHED
        self.note  = note # string with user comment on sources for parameters
    def DefineParamsID(self,matdatstr,matdata,imat):
        """ Initialize material parameters from using values from material database file.
            Usage: DefineParamsID(self,matdatstr,matdata,imat):
            Inputs: material name string, matdata DataFrame, matdata index structure
        """
        idx = np.where(matdata.loc[:,'Material'].values == matdatstr)[0]
        #print(idx,matdatstr)
        if len(idx) > 0:
            #print('Material index=',idx,matdatstr)
            self.name     = matdatstr
            self.rho0     = matdata.iloc[idx[0],imat.rho0] # mks
            self.v0       = 1/self.rho0 # m3/kg
            self.c0       = matdata.iloc[idx[0],imat.c0] # m/s
            self.s1       = matdata.iloc[idx[0],imat.s1] # [-]
            self.s2       = matdata.iloc[idx[0],imat.s2] # [-]
            self.d        = matdata.iloc[idx[0],imat.d] # [-]
            self.g0       = matdata.iloc[idx[0],imat.g0] # [-]
            self.q        = matdata.iloc[idx[0],imat.q] # [-]        
            self.ihed.id  = int(matdata.iloc[idx[0],imat.ihed]) # [integer] -1 indicates no data in IHED
            self.note     = matdata.iloc[idx[0],imat.note] # string with parameter source information
        else:
            print('WARNING: Cannot find this material in the database: ',matdatstr)
    def GetIHED(self,formflag=1,upmin=0.,upmax=1E99,id2=-1,moredata=[-1],uselocalbool=False):
        """Gather material data points from IHED online database.
           Usage: GetIHED(self,formflag=1,upmin=0.,upmax=1E99,id2=-1,moredata=[-1],uselocalbool=False):
           Inputs: Optional formflag 1 for linear fit.
                                     2 for quadratic fit. 
                                     3 for mod. universal liquid Hugoniot. Us = A + B Up − C Up exp(−D Up)
                   Optional upmin & upmax for fit is in m/s.
                   Optional second IHED material ID number to add for Hugoniot fit (id2) 
                   OR user-supplied data with format moredata = [rho0, [uparr], [usarr]]
                   WARNING: Cannot use both id2 and moredata simultaneously.
                   Optional uselocalbool boolean flag to read and save IHED data in local directory.
        """
        if formflag not in [1,2,3]:
            # must be 1 or 2 or 3. 
            # 1 is line Us = c0 + s1 * up
            # 2 is quadratic Us = c0 + s1 * up + s2 * up^2
            # 3 is modified universal liquid equation Us = A + B Up − C Up exp(−D Up)
            #      A = c0 B=s1 C=s2, D=D material parameters
            print('Hugoniot flag=',formflag)            
            print('Hugoniot forms can only be line (1) or quadratic (2) or universal Hugoniot (3).')
            return
        if self.ihed.id > -1:
            if uselocalbool: # check if local copy of IHED data exists; if not load it
                ihedfname = 'database-ihed/IHED-'+str(self.ihed.id)+'.txt'
                if not os.path.isdir('database-ihed/'):
                    os.mkdir('database-ihed')
                    print('made local directory: database-ihed')
                if not os.path.exists(ihedfname):
                    print('fetching IHED table from web server')
                    url = 'http://www.ihed.ras.ru/rusbank/plaintext.php?substid='+str(self.ihed.id)+'&type=0'
                    hds = {'user-agent':'Mozilla/5.0'}
                    request = urllib.request.Request(url,headers=hds)
                    response = urllib.request.urlopen(request).read()
                    li=str(response).split("\\n")
                    #print(li)
                    with open(ihedfname, 'w') as fp:
                        # save a local copy of the IHED data file
                        print('writing local copy of IHED table')
                        # a little cleanup
                        tmp = str(response) # converts each newline into characters \n
                        tmp = tmp[2::] # get rid of extra b'
                        clean_response = tmp[0:len(tmp)-1] # get rid of ending '
                        # split by \n and then join with newlines
                        fp.write('\n'.join(clean_response.split("\\n")))  
                        tmp=''
                        clean_response=''
                else:
                    # read in the local copy
                    #print('reading in local copy of IHED table #',self.ihed.id)
                    with open(ihedfname) as fp:
                        li = fp.readlines() # this is a list for each line
                        #print(li)
                        #return # for debugging
            else:
                # read IHED data from web server
                url = 'http://www.ihed.ras.ru/rusbank/plaintext.php?substid='+str(self.ihed.id)+'&type=0'
                hds = {'user-agent':'Mozilla/5.0'}
                request = urllib.request.Request(url,headers=hds)
                response = urllib.request.urlopen(request).read()
                li=str(response).split("\\n")
                #print(li)
                #return # for debugging
            df={} # initialize variable for table data frame 
            # IHED table should be split into lines in list variable li
            value=[]
            material=li[1] # second line has material name and initial density in g/cm3
            materialname=material.split(",")[0] # extract name
            for k in li:
                if k==li[0] or k==li[1] or k==li[2] or k==li[3]:
                    continue
                y=k.split()
                if y[0]=="References:":
                    break # stop gathering data 
                value.append(y)
                df[material]= DataFrame(value)
            for key in df: # some data have remarks, some don't
                if df[key].shape[1]==9:
                    df[key].columns=['m','Up','Us','Pressure','R/R0_ratio','Density','E-E0','Rem','Reference']
                elif df[key].shape[1]==8:
                    df[key].columns=['m','Up','Us','Pressure','R/R0_ratio','Density','E-E0','Reference']
            dataset = []
            for key in df:
                x = np.array(df[key]['Up'].astype('float64'))
                y = np.array(df[key]['Us'].astype('float64'))    
                p = np.array(df[key]['Pressure'].astype('float64'))
                d = np.array(df[key]['Density'].astype('float64'))    
                m = np.array(df[key]['m'].astype('float64'))
                densitytemp=key.split()
                density=float(densitytemp[-2])
            self.ihed.rho0    = density * 1000. # kg/m3
            self.ihed.uparr   = x * 1000. # m/s
            self.ihed.usarr   = y * 1000. # m/s
            self.ihed.varr    = 1/(d*1000.) # m3/kg
            self.ihed.parr    = p * 1E9 # Pa
            self.ihed.marr    = m # [-]
            self.ihed.matname = materialname
            #print('Got IHED data: ',materialname,self.ihed.id,key)
            # fit the nonporous data
            ind = np.where( (self.ihed.uparr > upmin) & (self.ihed.uparr < upmax) & (self.ihed.marr == 1.0))[0]
            if self.name == 'Ice':
                ind = np.where( (self.ihed.uparr > upmin) & (self.ihed.uparr < upmax) & (self.ihed.marr == 1.093))[0]
            if len(ind)>0:
                if formflag in [1,2]:
                    # line or quadratic
                    fitparams = np.polyfit(self.ihed.uparr[ind],self.ihed.usarr[ind],formflag) # mks
                if formflag == 3:
                    # curve fit to universal liquid Hugoniot
                    p0=[min(self.ihed.usarr),1.2,2.,.3e-3]
                    fitparams, con = curve_fit(UniversalHugoniot, self.ihed.uparr[ind], self.ihed.usarr[ind],p0=p0)
                    #print('p0=',p0)
                    #print('Curvefit=',fitparams,con)
            elif moredata[0] == -1:
                print('ERROR: No points in up range to fit!')
                return
            #print(fitparams)
            if len(ind)>0:
                if formflag == 1:
                    self.ihed.c0 = fitparams[1] # m/s
                    self.ihed.s1 = fitparams[0] # [-]
                    self.ihed.s2 = 0. # set to zero to keep variables clean
                    self.ihed.d  = 0. # set to zero to keep variables clean
                elif formflag == 2:
                    self.ihed.c0 = fitparams[2] # m/s
                    self.ihed.s1 = fitparams[1] # [-]
                    self.ihed.s2 = fitparams[0] # (m/s)^-1
                    self.ihed.d  = 0. # set to zero to keep variables clean
                elif formflag == 3: # universal liquid Hugoniot
                    self.ihed.c0 = fitparams[0] # m/s
                    self.ihed.s1 = fitparams[1] # [-]
                    self.ihed.s2 = fitparams[2] # [-]
                    self.ihed.d  = fitparams[3] # (m/s)^-1              
            if id2>-1:
                # get a second set of data from IHED
                if uselocalbool: # check if local copy of IHED data exists; if not load it
                    ihedfname = 'database-ihed/IHED-'+str(id2)+'.txt'
                    if not os.path.isdir('database-ihed/'):
                        os.mkdir('database-ihed')
                        print('made local directory: database-ihed')
                    if not os.path.exists(ihedfname):
                        print('fetching IHED table from web server')
                        url = 'http://www.ihed.ras.ru/rusbank/plaintext.php?substid='+str(id2)+'&type=0'
                        hds = {'user-agent':'Mozilla/5.0'}
                        request = urllib.request.Request(url,headers=hds)
                        response = urllib.request.urlopen(request).read()
                        li=str(response).split("\\n")
                        #print(li)
                        with open(ihedfname, 'w') as fp:
                            # save a local copy of the IHED data file
                            print('writing local copy of IHED table')
                            # a little cleanup
                            tmp = str(response) # converts each newline into characters \n
                            tmp = tmp[2::] # get rid of extra b'
                            clean_response = tmp[0:len(tmp)-1] # get rid of ending '
                            # split by \n and then join with newlines
                            fp.write('\n'.join(clean_response.split("\\n")))  
                            tmp=''
                            clean_response=''
                    else:
                        # read in the local copy
                        #print('reading in local copy of IHED table #',id2)
                        with open(ihedfname) as fp:
                            li = fp.readlines() # this is a list for each line
                            #print(li)
                            #return # for debugging
                else:
                    # read IHED data from web server
                    url = 'http://www.ihed.ras.ru/rusbank/plaintext.php?substid='+str(id2)+'&type=0'
                    hds = {'user-agent':'Mozilla/5.0'}
                    request = urllib.request.Request(url,headers=hds)
                    response = urllib.request.urlopen(request).read()
                    li=str(response).split("\\n")
                    #print(li)
                    #return # for debugging
                df={}
                value=[]
                material=li[1] # first line has material name and initial density in g/cm3
                materialname=material.split(",")[0] # extract name
                for k in li:
                    if k==li[0] or k==li[1] or k==li[2] or k==li[3]:
                        continue
                    y=k.split()
                    if y[0]=="References:":
                        break # stop gathering data 
                    value.append(y)
                    df[material]= DataFrame(value)
                for key in df: # some data have remarks, some don't
                    if df[key].shape[1]==9:
                        df[key].columns=['m','Up','Us','Pressure','R/R0_ratio','Density','E-E0','Rem','Reference']
                    elif df[key].shape[1]==8:
                        df[key].columns=['m','Up','Us','Pressure','R/R0_ratio','Density','E-E0','Reference']
                dataset = []
                for key in df:
                    x = np.array(df[key]['Up'].astype('float64'))
                    y = np.array(df[key]['Us'].astype('float64'))    
                    p = np.array(df[key]['Pressure'].astype('float64'))
                    d = np.array(df[key]['Density'].astype('float64'))    
                    m = np.array(df[key]['m'].astype('float64'))
                    densitytemp=key.split()
                    density=float(densitytemp[-2])
                self.ihed2.rho0    = density * 1000. # kg/m3
                self.ihed2.uparr   = x * 1000. # m/s
                self.ihed2.usarr   = y * 1000. # m/s
                self.ihed2.varr    = 1/(d*1000.) # m3/kg
                self.ihed2.parr    = p * 1E9 # Pa
                self.ihed2.marr    = m # [-]
                self.ihed2.matname = materialname
                self.ihed2.id = id2
                print('Got IHED data: ',materialname,self.ihed2.id,key)
                # fit the nonporous data from both datasets in desired particle velocity range
                ind2 = np.where( (self.ihed2.uparr > upmin) & (self.ihed2.uparr < upmax) & (self.ihed2.marr == 1.0))[0]
                ind1 = np.where( (self.ihed.uparr > upmin) & (self.ihed.uparr < upmax) & (self.ihed.marr == 1.0))[0]
                fitup = np.append(self.ihed.uparr[ind1],self.ihed2.uparr[ind2])
                fitus = np.append(self.ihed.usarr[ind1],self.ihed2.usarr[ind2])
                if formflag in [1,2]:
                    fitparams = np.polyfit(fitup,fitus,formflag) # mks
                if formflag == 3:
                    ## CHANGE TO CURVEFIT
                    p0=[min(fitus),1.2,2.,.3e-3]
                    fitparams, con = curve_fit(UniversalHugoniot, fitup, fitus, p0=p0)
                    #print('p0=',p0)
                    #print('Curvefit=',fitparams,con)
               #print(fitparams)
                if formflag == 1:
                    self.ihed2.c0 = fitparams[1] # m/s
                    self.ihed2.s1 = fitparams[0] # [-]
                    self.ihed2.s2 = 0. # set to zero to keep variables clean
                    self.ihed2.d  = 0. # set to zero to keep variables clean
                elif formflag == 2:
                    self.ihed2.c0 = fitparams[2] # m/s
                    self.ihed2.s1 = fitparams[1] # [-]
                    self.ihed2.s2 = fitparams[0] # (m/s)^-1
                    self.ihed2.d  = 0 # set to zero to keep variables clean
                else:
                    # formflag == 3
                    self.ihed2.c0 = fitparams[0] # m/s a
                    self.ihed2.s1 = fitparams[1] # [-] b
                    self.ihed2.s2 = fitparams[2] # [-] c
                    self.ihed2.d  = fitparams[3] # (m/s)^-1 d
            if moredata[0] != -1:
                # get a second set of data from [rho0,uparr,usarr] array
                print('Including user data')
                nmoredata=len(moredata)
                tmp=int((len(moredata)-1)/2)
                x=moredata[1:tmp+1] # up
                y=moredata[tmp+1::] # us
                self.ihed2.rho0    = moredata[0] # kg/m3
                self.ihed2.uparr   = x # m/s
                self.ihed2.usarr   = y # m/s
                self.ihed2.varr    = (1/self.ihed2.rho0)*(1-self.ihed2.uparr/self.ihed2.usarr) # m3/kg
                self.ihed2.parr    = self.ihed2.rho0*self.ihed2.uparr*self.ihed2.usarr # Pa
                self.ihed2.marr    = np.ones(tmp)
                self.ihed2.matname = self.name
                self.ihed2.id      = 9999
                # fit the nonporous data from both datasets in desired particle velocity range
                ind2 = np.where( (self.ihed2.uparr > upmin) & (self.ihed2.uparr < upmax) & (self.ihed2.marr == 1.0))[0]
                ind1 = np.where( (self.ihed.uparr > upmin) & (self.ihed.uparr < upmax) & (self.ihed.marr == 1.0))[0]
                fitup = np.append(self.ihed.uparr[ind1],self.ihed2.uparr[ind2])
                fitus = np.append(self.ihed.usarr[ind1],self.ihed2.usarr[ind2])
                #fitparams = np.polyfit(fitup,fitus,formflag) # mks
                if formflag in [1,2]:
                    fitparams = np.polyfit(fitup,fitus,formflag) # mks
                if formflag == 3:
                    p0=[min(fitus),1.2,2.,.3e-3]
                    fitparams, con = curve_fit(UniversalHugoniot, fitup, fitus, p0=p0)
                    #print('p0=',p0)
                    #print('Curvefit=',fitparams,con)
                #print(fitparams)
                if formflag == 1:
                    self.ihed2.c0 = fitparams[1] # m/s
                    self.ihed2.s1 = fitparams[0] # [-]
                    self.ihed2.s2 = 0. # set to zero to keep variables clean
                    self.ihed2.d  = 0. # set to zero to keep variables clean
                elif formflag == 2:
                    self.ihed2.c0 = fitparams[2] # m/s
                    self.ihed2.s1 = fitparams[1] # [-]
                    self.ihed2.s2 = fitparams[0] # (m/s)^-1
                    self.ihed2.d  = 0. # set to zero to keep variables clean
                elif formflag == 3:
                    self.ihed2.c0 = fitparams[0] # m/s
                    self.ihed2.s1 = fitparams[1] # [-]
                    self.ihed2.s2 = fitparams[2] # (m/s)^-1
                    self.ihed2.d  = fitparams[3] # (m/s)^-1                 
        else:
            print('No IHED material ID number entered for this material.')
    def PlotIHED(self,savebool=False,fname=''):
        """ Plot IHED data and fit.
            Usage: PlotIHED(self,savebool=False,fname=''):
            Inputs: Optional 0 or 1 save PDF figure boolean.
                    Optional figure name string. Default is 'IHED-plot-'+self.name+'-v'+__version__+'.pdf'
        """
        fig = plt.figure() # initialize the figure object
        paramstring = 'r$_0$='+ClStr(self.rho0/1.e3)+' (g/cm$^3$), Us (km/s)='+ClStr(self.c0/1.e3)+'+'+ClStr(self.s1)+'up'
        if self.s2 != 0 and self.d == 0: # s2 is s/m -> s/km
            if self.s2<0:
                paramstring += str(self.s2*1.e3)+'up$^2$'
            if self.s2>0:
                paramstring += '+'+str(self.s2*1.e3)+'up$^2$'
        if self.d != 0: # s2 is dimless and d is s/m -> s/km
            if self.s2<0:
                paramstring += ClStr(self.s2)+'up exp(-'+str(self.d*1.e3)
            if self.s2>0:
                paramstring += '+'+ClStr(self.s2)+'up exp(-'+str(self.d*1.e3)+'up)'

        if self.ihed.d != 0:
            lab_txt = 'IHED univeral liquid fit'
        elif self.ihed.s2 != 0:
            lab_txt = 'IHED quadratic fit'
        else:
            lab_txt = 'IHED linear fit'
        if self.ihed2.id > 0:
            # combine 2 datasets
            ind = np.where(self.ihed.marr == 1.0)[0] # don't plot porous data 
            up1=self.ihed.uparr[ind]
            us1=self.ihed.usarr[ind]
            ind = np.where(self.ihed2.marr == 1.0)[0] # don't plot porous data 
            up2=self.ihed2.uparr[ind]
            us2=self.ihed2.usarr[ind]
            up = np.append(up1,up2)
            us = np.append(us1,us2)
            c0=self.ihed2.c0 # ihed2 contains fit to both datasets
            s1=self.ihed2.s1
            s2=self.ihed2.s2
            d=self.ihed2.d
        else:
            ind = np.where(self.ihed.marr == 1.0)[0] # don't plot porous data 
            if self.name == 'Ice':
                ind = np.where(self.ihed.marr == 1.093)[0] # for some reason IHED did not use ice's actual density
            up1=self.ihed.uparr[ind]
            us1=self.ihed.usarr[ind]
            c0=self.ihed.c0
            s1=self.ihed.s1
            s2=self.ihed.s2
            d=self.ihed.d
            up=up1
            us=us1
        plt.scatter(up1/1000.,us1/1000.,label='IHED '+self.ihed.matname)
        if self.ihed2.id > 0:
            if self.ihed2.id == 9999:
                plt.scatter(up2/1000.,us2/1000.,label='User '+self.ihed2.matname)
            else:
                plt.scatter(up2/1000.,us2/1000.,label='IHED '+self.ihed2.matname)
        uptmp = np.arange(max(up))
        if d==0:
            plt.plot(uptmp/1000.,(c0+s1*uptmp+s2*uptmp*uptmp)/1000.,label=lab_txt)
            diff = us-(c0+s1*up+s2*up*up)
        else:
            plt.plot(uptmp/1000.,(c0+s1*uptmp-s2*uptmp*np.exp(-d*uptmp))/1000.,label=lab_txt)            
            diff = us-(c0+s1*up-s2*up*np.exp(-d*up))
#        diff = self.ihed.usarr[ind]-(self.ihed.c0+self.ihed.s1*self.ihed.uparr[ind]+self.ihed.s2*self.ihed.uparr[ind]*self.ihed.uparr[ind])
#        diff = us-(c0+s1*up+s2*up*up)
        if self.ihed2.id > -1:
            print(self.ihed.matname+' '+self.ihed2.matname+' fit (mks) c0=',self.ihed2.c0,' s1=',self.ihed2.s1,' s2=',self.ihed2.s2,' d=',self.ihed2.d)
        else:
            print(self.ihed.matname+' fit (mks) c0=',self.ihed.c0,' s1=',self.ihed.s1,' s2=',self.ihed.s2,' d=',self.ihed.d)
        print('N='+str(len(up))+'. Fit stdev (km/s)='+str(np.std(diff/1000.)))
        # plot main Hugoniot fit from the materials database csv file
        #uptmp = np.arange(max(up))
        if self.d == 0:
            if self.s2 ==0:
                labelstr = 'Primary Linear Hugoniot fit'
            else:
                labelstr = 'Primary Quadratic Hugoniot fit'
            plt.plot(uptmp/1000.,(self.c0+self.s1*uptmp+self.s2*uptmp*uptmp)/1000.,'--',label=labelstr)
        else:
            plt.plot(uptmp/1000.,(self.c0+self.s1*uptmp-self.s2*uptmp*np.exp(-self.d*uptmp))/1000.,'--',label='Primary UL Hugoniot fit')            
        print('Primary Hugoniot parameters c0,s1,s2,d=',self.c0,self.s1,self.s2,self.d)
        plt.xlabel('Particle Velocity (km/s)')
        plt.ylabel('Shock Velocity (km/s)')
        plt.title('Material = '+self.name+', IHED data '+self.ihed.matname+'\n'+paramstring)
        plt.legend()
        if savebool:
            if fname == '':
                fname='IHED-plot-'+self.name+'-v'+__version__+'.pdf'
            plt.savefig(fname,dpi=300)
        #plt.show()
        #fig.close()
        uptmp=[]
    def MakeHugoniot(self,uparr,ihedbool=False):
        """ Make Hugoniot from material data parameters
            Usage: MakeHugoniot(self,uparr,ihedbool=False):
            Form depends on number of non-zero material parameters: c0, s1, s2, D
            c0 s1 !=0 linear
            c0 s1 s2 !=0 quadratic
            c0 s1 s2 d !=0 modified universal liquid Hugoniot Us = (a + b * up - c * up *np.exp(-d*up))
            Isentrope and reshock Hugoniot use the same volume array as defined by the input particle velocity array.
            Inputs: particle velocity numpy array in m/s.
            Optional ihedbool boolean to use material database IHED fitted Hugoniot parameters instead of matdata.
        """
        if ihedbool:
            if self.ihed2.id>-1:
                c0=self.ihed2.c0 # use fit to combined datasets
                s1=self.ihed2.s1
                s2=self.ihed2.s2
                d =self.ihed2.d
            else:
                c0=self.ihed.c0 # use fit to 1st (only) IHED dataset
                s1=self.ihed.s1
                s2=self.ihed.s2
                d =self.ihed.d
        else:
            c0=self.c0 # use fit in materials database csv file
            s1=self.s1
            s2=self.s2
            d =self.d
        self.hug.uparr = uparr # array in m/s
        if d == 0:
            # line or quadratic
            self.hug.usarr = c0+s1*self.hug.uparr+s2*self.hug.uparr*self.hug.uparr # m/s
        else:
            # universal liquid Hugoniot a=c0, b=s1, c=s2, d=d
            self.hug.usarr = c0+s1*self.hug.uparr-s2*self.hug.uparr*np.exp(-d*self.hug.uparr) # m/s            
        self.hug.parr  = self.rho0*self.hug.uparr*self.hug.usarr # Pa
        self.hug.varr  = self.v0*(1-self.hug.uparr/self.hug.usarr) # m3/kg
        self.hug.earr  = self.e0+0.5*self.hug.parr*(self.v0-self.hug.varr) # J/kg
        if self.g0 <= 0.:
            print('Gamma value not set. Using defaults: g0=1, q=1')
            self.g0=1.
            self.q=1.
        self.hug.garr  = self.g0*np.power((self.rho0*self.hug.varr),self.q) # [-]
    def MakeMGIsentrope(self,pstart,usemgmodelbool=False,impvel=0.):
        """ Calculate the isentrope from pstart on the principal Hugoniot.
            MakeMGIsentrope(self,pstart,usemgmodelbool=False,impvel=0.)
            Default is to use the Hugoniot for isentrope. Flag uses the Mie-Grueneisen model (can be unstable).
            Usage: MakeMGIsentrope(self,pstart,usemgmodelbool=False):
            Uses the same volume and gamma arrays as in the principal Hugoniot.
            Input: EOS_Point object
            Optional boolean to use Hugoniot for the isentrope.
            Impact velocity to check for symmetric impact cases required for MG model.
            Returns: successful calculation boolean
        """
        MGIsuccess=False
        #Must create Hugoniot first
        if len(self.hug.parr) <= 10:
            MGIsuccess=False
            print('MakeMGIsentrope ERROR: need to generate a Hugoniot first.')
            return MGIsuccess
        if usemgmodelbool and (impvel==0):
            MGIsuccess=False
            print('MakeMGIsentrope ERROR: need to provide the impact velocity.')
            return MGIsuccess
        if not usemgmodelbool:
            self.isen.parr = np.copy(self.hug.parr)
            self.isen.varr = np.copy(self.hug.varr)
            self.isen.earr = np.copy(self.hug.earr)
            self.isen.uparr = np.copy(self.hug.uparr)
            self.isen.uparr2 = np.copy(self.hug.uparr)
            self.isen.garr = np.copy(self.hug.garr)
            MGIsuccess = True
            return MGIsuccess
        #print(pstart.p, max(self.hug.parr))
        if pstart.p >= max(self.hug.parr):
            MGIsuccess=False
            print('MakeMGIsentrope ERROR: pstart > max Hugoniot pressure')
            return MGIsuccess
        istart  = np.where(self.hug.parr > pstart.p)[0][0]  # ASSUMES THAT THE VOLUME ARRAY DECREASES WITH INCREASING INDEX; Pressure increases with increasing index
        #print('ISTART=',istart,' Pstart=',pstart.p)
        self.isen.parr = np.zeros(len(self.hug.varr))
        self.isen.varr = np.copy(self.hug.varr)
        self.isen.earr = np.zeros(len(self.hug.varr))
        self.isen.uparr = np.zeros(len(self.hug.varr))
        self.isen.uparr2 = np.zeros(len(self.hug.varr))
        self.isen.garr = np.copy(self.hug.garr)
        # calc for volumes on release
        #print('istart,v,p,up,e=',istart,pstart.v,pstart.p/1.e9,pstart.up/1.e3,pstart.e)
        #with open('log.txt', 'a') as fp: # debugging
        #    fp.write('pstart='+str(pstart.p)+' = '+str(self.hug.parr[istart])+'\n')
        #    fp.write('vstart='+str(pstart.v)+' = '+str(self.hug.varr[istart])+'\n')
        for i in range(istart-1,-1,-1):
            if i == istart-1:
                dv = (self.isen.varr[i]-pstart.v)
                print('istart down dv, startup = ',dv,pstart.up,self.isen.uparr[istart])
                # check for roundoff errors for symmetric impact
                print('impvel, up/impvel', impvel, ((2*pstart.up/impvel)))
                if (dv <= 0) or ((abs(2*pstart.up/impvel)-1)<1.e-5): 
                    # sometimes dv is a very tiny number; not catching positive tiny right now
                    # up vs impvel check is because there are floating precision problems
                    # case of symmetric impact messing up the math
                    print('symmetric i=',i)
                    self.isen.parr[i] = np.copy(pstart.p)
                    self.isen.earr[i] = np.copy(pstart.e)
                    self.isen.uparr[i]= np.copy(pstart.up)
                    self.isen.uparr2[i]= np.copy(pstart.up)
                    print('symm up = ',self.isen.uparr[i],self.isen.uparr2[i])
                else:
                    self.isen.parr[i] = (self.hug.parr[i]-(self.hug.earr[i]-pstart.e+pstart.p*dv/2.)*(self.hug.garr[i]/self.hug.varr[i])) / (1.-dv*self.hug.garr[i]/self.hug.varr[i]/2.)
                    self.isen.earr[i] = pstart.e-(pstart.p+self.isen.parr[i])*dv/2.
                    if (-(pstart.p-self.isen.parr[i])/(pstart.v-self.isen.varr[i])) > 0:
                        # no sqrt of negative numbers
                        self.isen.uparr[i]= pstart.up-(pstart.p-self.isen.parr[i])/np.sqrt(-(pstart.p-self.isen.parr[i])/(pstart.v-self.isen.varr[i]))
                        self.isen.uparr2[i]= pstart.up-(pstart.p-self.isen.parr[i])/np.sqrt(-(pstart.p-self.isen.parr[i])/(pstart.v-self.isen.varr[i]))
                        print('what is up with up1, up2 = ',self.isen.uparr[i],self.isen.uparr2[i])
            else: # not the first step on the release
                dv = (self.isen.varr[i]-self.isen.varr[i+1])
                #print('istart down rest dv = ',dv)
                gov = self.isen.garr[i]/self.isen.varr[i]
                self.isen.parr[i] = (self.hug.parr[i]-(self.hug.earr[i]-self.isen.earr[i+1]+self.isen.parr[i+1]*dv/2.)*(self.hug.garr[i]/self.hug.varr[i])) / (1.-dv*self.hug.garr[i]/self.hug.varr[i]/2.)
                self.isen.earr[i] = self.isen.earr[i+1]-(self.isen.parr[i+1]+self.isen.parr[i])*dv/2.
                if (-(self.isen.parr[i+1]-self.isen.parr[i])/(self.isen.varr[i+1]-self.isen.varr[i])) > 0:
                    # no sqrt of negative numbers
                    # these two equations for Uparr are identical along an isentrope
                    self.isen.uparr[i]= self.isen.uparr[i+1]-(self.isen.parr[i+1]-self.isen.parr[i])/np.sqrt(-(self.isen.parr[i+1]-self.isen.parr[i])/(self.isen.varr[i+1]-self.isen.varr[i]))
                    self.isen.uparr2[i]= self.isen.uparr2[i+1]-np.sqrt(-(self.isen.parr[i+1]-self.isen.parr[i])*(self.isen.varr[i+1]-self.isen.varr[i]))
                    #print('down side i, up1, up2, dv = ',i,self.isen.uparr[i],self.isen.uparr2[i],dv)
                else:
                    self.isen.uparr[i]= self.isen.uparr[i+1]
                    self.isen.uparr2[i]= self.isen.uparr2[i+1]
                    print('check negative i up1 up2 dv = ',i,self.isen.uparr[i],self.isen.uparr2[i],dv)

            #print(i,dv,self.isen.parr[i]/1.e9,self.isen.uparr[i]/1.e3,self.isen.earr[i])
        # calc for greater volumes for plotting purposes
        for i in range(istart,len(self.hug.parr)):
            if i == istart:
                dv = (pstart.v-self.isen.varr[i]) # check this!
                print('istart up dv = ',i,dv)
                print('upside impvel, up/impvel', impvel, ((2*pstart.up/impvel)))
#                if dv<=0: # sometimes a very tiny number; not catching positive tiny right now
                if (dv <= 0) or ((abs(2*pstart.up/impvel)-1)<1.e-5): 
                    # sometimes dv is a very tiny number; not catching positive tiny right now
                    # up vs impvel check is because there are floating precision problems
                    # case of symmetric impact messing up the math
                    print('symmetric i=',i)
                    # symmetric case check
                    self.isen.parr[i] = pstart.p
                    self.isen.earr[i] = pstart.e
                    self.isen.uparr[i]= pstart.up
                    self.isen.uparr2[i]= pstart.up
                else:  
                    self.isen.parr[i] = (self.hug.parr[i]-(self.hug.earr[i]-pstart.e+pstart.p*dv/2.)*(self.hug.garr[i]/self.hug.varr[i])) / (1.-dv*self.hug.garr[i]/self.hug.varr[i]/2.)
                    self.isen.earr[i] = pstart.e-(pstart.p+self.isen.parr[i])*dv/2.
                    if (-(pstart.p-self.isen.parr[i])/(pstart.v-self.isen.varr[i])) > 0:
                        # no sqrt of negative numbers
                        self.isen.uparr[i]= pstart.up-(pstart.p-self.isen.parr[i])/np.sqrt(-(pstart.p-self.isen.parr[i])/(pstart.v-self.isen.varr[i]))
                        self.isen.uparr2[i]= pstart.up-(pstart.p-self.isen.parr[i])/np.sqrt(-(pstart.p-self.isen.parr[i])/(pstart.v-self.isen.varr[i]))
            else:
                dv = -(self.isen.varr[i]-self.isen.varr[i-1])
                #print('istart up rest dv = ',i,dv)
                self.isen.parr[i] = (self.hug.parr[i]-(self.hug.earr[i]-self.isen.earr[i-1]+self.isen.parr[i-1]*dv/2.)*(self.hug.garr[i]/self.hug.varr[i])) / (1.-dv*self.hug.garr[i]/self.hug.varr[i]/2.)
                self.isen.earr[i] = self.isen.earr[i-1]+(self.isen.parr[i-1]+self.isen.parr[i])*dv/2.
#                if ((self.isen.parr[i-1]-self.isen.parr[i])/(self.isen.varr[i-1]-self.isen.varr[i])) > 0:
#                    # no sqrt of negative numbers
#                    # these two equations for Uparr are identical along an isentrope
#                    self.isen.uparr[i]= self.isen.uparr[i-1]-(self.isen.parr[i-1]-self.isen.parr[i])/np.sqrt((self.isen.parr[i-1]-self.isen.parr[i])/(self.isen.varr[i-1]-self.isen.varr[i]))
#                    self.isen.uparr2[i]= self.isen.uparr2[i-1]+np.sqrt((self.isen.parr[i-1]-self.isen.parr[i])*(self.isen.varr[i-1]-self.isen.varr[i]))
                if (-(self.isen.parr[i-1]-self.isen.parr[i])/(self.isen.varr[i-1]-self.isen.varr[i])) > 0:
                    # no sqrt of negative numbers
                    # these two equations for Uparr are identical along an isentrope
                    self.isen.uparr[i]= self.isen.uparr[i-1]-(self.isen.parr[i-1]-self.isen.parr[i])/np.sqrt(-(self.isen.parr[i-1]-self.isen.parr[i])/(self.isen.varr[i-1]-self.isen.varr[i]))
                    self.isen.uparr2[i]= self.isen.uparr2[i-1]+np.sqrt(-(self.isen.parr[i-1]-self.isen.parr[i])*(self.isen.varr[i-1]-self.isen.varr[i]))
                else:
                    # MG model is failing at higher pressures; this happens all the time so don't call it a fail
                    # but instead fill the upper parts of the isentrope with NaNs
                    print('up side negative sqrt nans i, up(i-1)=',i, self.isen.uparr[i-1])
                    self.isen.parr[i::] = np.nan
                    self.isen.earr[i::] = np.nan
                    self.isen.uparr[i::]= np.nan
                    self.isen.uparr2[i::]= np.nan
                    MIGsuccess = True
                    return MGIsuccess
            #print(i,dv,self.isen.parr[i]/1.e9,self.isen.uparr[i]/1.e3,self.isen.earr[i])
            # made it here without crashing; hopefully the arrays are all OK
            MGIsuccess=True
        return MGIsuccess
    #------------------------------------------------------------------------------------------
    def MakeReshockHug(self,pstart,usemgmodelbool=False):
        """ Calculate reshock Hugoniot from pstart on the principal Hugoniot.
            Default is to use the Hugoniot for reshock. Flag uses the Mie-Grueneisen model (can be unstable).
            Usage: MakeReshockHug(self,pstart,usemgmodelbool=False):
            Uses the same volume and gamma arrays as in the principal Hugoniot.
            Input: Pressure in Pa.
            Returns: successful calculation boolean
        """
        MGRsuccess=False
        if not usemgmodelbool:
            self.reshock.parr = np.copy(self.hug.parr)
            self.reshock.varr = np.copy(self.hug.varr)
            self.reshock.earr = np.copy(self.hug.earr)
            self.reshock.uparr = 2*pstart.up-np.copy(self.hug.uparr) # contains reshock to lower up
            self.reshock.uparr2 = np.copy(self.hug.uparr) # contains original Hugoniot up
            self.reshock.garr = np.copy(self.hug.garr)
            ind  = np.where(self.reshock.parr < pstart.p)[0]  # ASSUMES THAT THE VOLUME ARRAY DECREASES WITH INCREASING INDEX; Pressure increases with increasing index
            self.reshock.uparr[ind]=np.nan
            self.reshock.parr[ind]=np.nan
            self.reshock.earr[ind]=np.nan
            MGRsuccess = True
            return MGRsuccess
        #print('Making reshock Hugoniot from P=',pstart.p/1.e9,' GPa')
        # ASSUMES THAT THE VOLUME ARRAY DECREASES WITH INCREASING INDEX; Pressure increases with increasing index
        ind  = np.where((self.hug.parr > pstart.p))[0]
        #Must create Hugoniot first
        self.reshock.parr = np.zeros(len(self.hug.parr))
        self.reshock.varr = np.copy(self.hug.varr )
        self.reshock.garr = np.copy(self.hug.garr)
        self.reshock.earr = np.zeros(len(self.hug.parr))
        self.reshock.uparr = np.zeros(len(self.hug.parr))
        self.reshock.uparr2 = np.zeros(len(self.hug.parr))
        # calc reshock Hugoniot
        dv = pstart.v-self.reshock.varr[ind]
        self.reshock.parr[ind] = (self.hug.parr[ind]+(pstart.p-self.hug.parr[ind])*(self.v0-self.reshock.varr[ind])*(self.reshock.garr[ind]/self.reshock.varr[ind]/2))/(1-self.reshock.garr[ind]/self.reshock.varr[ind]/2*(pstart.v-self.reshock.varr[ind]))
        self.reshock.earr[ind] = pstart.e+0.5*(self.reshock.parr[ind]+pstart.p)*dv
        # check for bad points
        ind2 = np.where((self.reshock.parr[ind]-pstart.p)*(pstart.v-self.reshock.varr[ind]) > 0)[0]
        self.reshock.uparr[ind[ind2]] = pstart.up-np.sqrt((self.reshock.parr[ind[ind2]]-pstart.p)*(pstart.v-self.reshock.varr[ind[ind2]]))
        self.reshock.uparr2[ind[ind2]] = pstart.up+np.sqrt((self.reshock.parr[ind[ind2]]-pstart.p)*(pstart.v-self.reshock.varr[ind[ind2]]))
        ind2 = np.where((self.reshock.parr[ind]-pstart.p)*(pstart.v-self.reshock.varr[ind]) <= 0)[0]
        if len(ind2) > 0:
            # MG model is failing
            self.reshock.parr[ind[ind2]] = np.nan
            self.reshock.earr[ind[ind2]] = np.nan
            self.reshock.uparr[ind[ind2]]= np.nan
            self.reshock.uparr2[ind[ind2]]= np.nan 
        ind  = np.where(self.reshock.parr < pstart.p)[0]  # ASSUMES THAT THE VOLUME ARRAY DECREASES WITH INCREASING INDEX; Pressure increases with increasing index
        if len(ind) > 0:
            self.reshock.uparr[ind]=np.nan
            self.reshock.uparr2[ind]=np.nan
            self.reshock.parr[ind]=np.nan
            self.reshock.earr[ind]=np.nan
        # check for reshock Hugoniot rollover
        hugshift = np.roll(self.reshock.parr,1)
        ind = np.where((self.reshock.parr-hugshift)<0)[0]
        if len(ind)>0:
            self.reshock.uparr[ind]=np.nan
            self.reshock.uparr2[ind]=np.nan
            self.reshock.parr[ind]=np.nan
            self.reshock.earr[ind]=np.nan
        MGRsuccess = True # made it to end without crashing
        return MGRsuccess
    #------------------------------------------------------------------------------------------
    def PlotCurves(self,pstart,savebool=False,fname=''):
        """ Plot principal Hugoniot and isentropes and reshock Hugoniot from pstart.
            Usage: PlotCurves(self,pstart,savebool=False,fname=''):
            Inputs: initial pressure for isentrope & reshock in Pa.
                    Optional savebool boolean to save PDF. 
                    Optional figure filename, default is 'EOS-plots-'+self.name+'-v'+__version__+'.pdf'
        """
        paramstring = 'r$_0$='+ClStr(self.rho0/1.e3)+' (g/cm$^3$), Us (km/s)='+ClStr(self.c0/1.e3)+'+'+ClStr(self.s1)+'up'
        if self.s2 != 0 and self.d == 0: # s2 is s/m -> s/km
            if self.s2<0:
                paramstring += str(self.s2*1.e3)+'up$^2$'
            if self.s2>0:
                paramstring += '+'+str(self.s2*1.e3)+'up$^2$'
        if self.d != 0: # s2 is dimless and d is s/m -> s/km
            if self.s2<0:
                paramstring += ClStr(self.s2)+'up exp(-'+str(self.d*1.e3)+'up)'
            if self.s2>0:
                paramstring += '+'+ClStr(self.s2)+'up exp(-'+str(self.d*1.e3)+'up)'
            
        plt.rcParams["figure.figsize"] = (12,10)
        plt.rc('font', size=15)
        self.im1.p=pstart
        self.im1.v=1/np.interp(pstart,self.hug.parr,1/self.hug.varr)
        self.im1.e=np.interp(pstart,self.hug.parr,self.hug.earr)
        MGIsuccess = self.MakeMGIsentrope(self.im1)
        #print('MGI success flag = ',MGIsuccess)
        self.MakeReshockHug(self.im1)
        rhomax = np.interp(pstart*3,self.hug.parr,1/self.hug.varr)/1.e3
        rhomin = 0.8*(self.rho0/1.e3)
        #print(rhomax)
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
        st = fig.suptitle(r' '+self.name+'\n'+paramstring)
        ax1.plot(1/self.hug.varr/1.e3,self.hug.parr/1.e9,label='Hugoniot')
        if MGIsuccess: 
            ax1.plot(1/self.isen.varr/1.e3,self.isen.parr/1.e9,label='Isentrope')
        ax1.plot(1/self.reshock.varr/1.e3,self.reshock.parr/1.e9,'--',label='Reshock Hug.')
        ax1.set_xlabel('Density (g/cm$^3$)')
        ax1.set_ylabel('Pressure (GPa)')
        ax1.set_ylim([0,pstart*3/1.e9])
        ax1.set_xlim([rhomin,rhomax])
        ax1.legend()
        
        ind=np.where(self.hug.parr < 2*pstart)[0]
        ax2.plot(1/self.hug.varr[ind]/1.e3,self.hug.earr[ind]/1.e6,label='Hugoniot')
        if MGIsuccess: 
            ax2.plot(1/self.isen.varr[ind]/1.e3,self.isen.earr[ind]/1.e6,label='Isentrope')
        ax2.plot(1/self.reshock.varr[ind]/1.e3,self.reshock.earr[ind]/1.e6,'--',label='Reshock Hug.')
        ax2.set_xlabel('Density (g/cm$^3$)')
        ax2.set_xlim([rhomin,rhomax])
        ax2.set_ylabel('Energy (MJ/kg)')
        ax2.legend()

        ax3.plot(1/self.hug.varr/1.e3,self.hug.garr,label='$\gamma(v)$')
        ax3.plot(1/self.hug.varr/1.e3,self.g0*self.hug.varr/self.v0,'--',label='$\gamma_0$/$v_0$=const')
        ax3.set_xlabel('Density (g/cm$^3$)')
        ax3.legend()
        ax3.set_ylabel('Mie Grueneisen Parameter')
        ax3.set_xlabel('Density (g/cm$^3$)')
        ax3.set_xlim([rhomin,rhomax])

        ind=np.where(self.hug.parr < 3*pstart)[0]
        up0 = np.interp(pstart,self.hug.parr,self.hug.uparr)
        ax4.plot(self.hug.uparr/1.e3,self.hug.parr/1.e9,label='Hugoniot')
        if MGIsuccess: 
            ax4.plot(up0/1.e3+self.isen.uparr/1.e3,self.isen.parr/1.e9,label='Isentrope')
        #ax4.plot(up0/1.e3+self.isen.uparr2/1.e3,self.isen.parr/1.e9,'--',label='isen2')
        ax4.plot(up0/1.e3+self.reshock.uparr/1.e3,self.reshock.parr/1.e9,label='Reshock A')
        ax4.plot(up0/1.e3+self.reshock.uparr2/1.e3,self.reshock.parr/1.e9,'--',label='Reshock B')
        ax4.set_ylim([0,pstart*3/1.e9])
        ax4.set_xlim([0,up0*2/1.e3])
        ax4.legend()
        ax4.set_xlabel('Particle Velocity (km/s)')
        ax4.set_ylabel('Pressure (GPa)')
        #print(mat3.isen.uparr2)

        fig.tight_layout()
        if savebool:
            # now an input variable: fname='EOS-plots-'+self.name+'-v'+__version__+'.pdf'
            if fname == '':
                fname='EOS-plots-'+self.name+'-v'+__version__+'.pdf'
            fig.savefig(fname,dpi=300)
        #plt.show()
        plt.close(fig)
        return fig
    #------------------------------------------------------------------------------------------
    def IM_match(self,matb,vel=0.,pstart=0,usemgmodelbool=False):
        """ Plot principal Hugoniot and isentropes and reshock Hugoniot from pstart.
            Usage: PlotCurves(self,pstart,savebool=False,fname=''):
            Inputs: initial pressure for isentrope & reshock in Pa.
                    Optional savebool boolean to save PDF. 
                    Optional figure filename, default is 'EOS-plots-'+self.name+'-v'+__version__+'.pdf'
        """
        IM_match_success = False
        # mat1.IM_match(mat2,up=up,vel=vel) calculates mat1 -> mat2 at vel with Hugoniot up array
        if (pstart==0) and (vel>0):
            # This is the first pair of materials: mata -> matb at vel (mks)
            # find index and initialize Hugoniot for each material layer
            res = Intersection(matb.hug.uparr,matb.hug.parr,vel-self.hug.uparr,self.hug.parr)
            print('1st IM: ',vel,res[0],res[1]/1.e9)
            self.im1.up=res[0][0] # m/s
            self.im1.p=res[1][0] # Pa
            self.im1.v=1./(np.interp(res[0],self.hug.uparr,1./self.hug.varr)[0]) # assumes f is monotonic and increasing
            self.im1.e=np.interp(res[0],self.hug.uparr,self.hug.earr)[0] # assumes f is monotonic and increasing
            matb.im1.up=res[0][0]
            matb.im1.p=res[1][0]
            matb.im1.v=1./(np.interp(res[0],matb.hug.uparr,1./matb.hug.varr)[0]) # assumes f is monotonic and increasing
            matb.im1.e=np.interp(res[0],matb.hug.uparr,matb.hug.earr)[0] # assumes f is monotonic and increasing
            IM_match_success = True
            return IM_match_success
        
        # mat2.IM_match(mat3,pstart=mat2.im1) calculates reshock/release from mat2 into mat
        if (pstart!=0):
           # This is a target pair of materials mata.im1 release/reshock into matb initially at rest & 0 pressure
            matbhug_patup = np.interp(self.im1.up,matb.hug.uparr,matb.hug.parr) # monotonic increasing density
            print('matbhug_patup (GPa) = ',matbhug_patup/1.e9)
            if matbhug_patup > self.im1.p:
                # reshock mata into matb
                reshock_success = self.MakeReshockHug(self.im1,usemgmodelbool=usemgmodelbool)
                print('RESHOCK SUCCESS = ',reshock_success)
                #self.reshock.uparr is decreasing
                res_reshock = Intersection(matb.hug.uparr,matb.hug.parr,np.flip(self.reshock.uparr),np.flip(self.reshock.parr))
                #plt.plot(matb.hug.uparr,matb.hug.parr)
                #plt.plot(self.reshock.uparr,self.reshock.parr))
                print('res_reshock = ',res_reshock[0],res_reshock[1]/1.e9,len(res_reshock))
                if len(res_reshock[0])>0:
                    #with open('log.txt', 'a') as fp: # debugging
                    #    fp.write('CHECK res_reshock vel,up,P= '+str(vel)+' '+str(res_reshock[0][0])+' '+str(res_reshock[1][0]/1.e9)+'\n')
                    #    print('CHECK res_reshock vel,up,P= '+str(vel)+' '+str(res_reshock[0][0])+' '+str(res_reshock[1][0]/1.e9)+'\n')
                    self.im2.up=res_reshock[0][0] # m/s
                    self.im2.p=res_reshock[1][0] # Pa
                    # reshock up array is decreasing, so flip
                    self.im2.v=1./(np.interp(res_reshock[0],np.flip(self.reshock.uparr),1./np.flip(self.reshock.varr))[0]) # assumes f is monotonic and increasing
                    self.im2.e=np.interp(res_reshock[0],np.flip(self.reshock.uparr),np.flip(self.reshock.earr))[0] # assumes f is monotonic and increasing
                    matb.im1.up=res_reshock[0][0]
                    matb.im1.p=res_reshock[1][0]
                    matb.im1.v=1./(np.interp(res_reshock[0],matb.hug.uparr,1./matb.hug.varr)[0]) # assumes f is monotonic and increasing
                    matb.im1.e=np.interp(res_reshock[0],matb.hug.uparr,matb.hug.earr)[0] # assumes f is monotonic and increasing
                    IM_match_success = True
                    return IM_match_success
                else:
                    # something wrong with mata Mie Grueneisen model
                    #with open('log.txt', 'a') as fp: # debugging
                    #    fp.write('ERROR res_reshock vel mat= '+str(vel)+self.name+'\n')
                    #    print('ERROR res_reshock vel mat= '+str(vel)+self.name+'\n')
                    self.im2.up=0 # m/s
                    self.im2.p=0 # Pa
                    self.im2.v=0
                    self.im2.e=0
                    matb.im1.up=0
                    matb.im1.p=0
                    matb.im1.v=0
                    matb.im1.e=0
                    IM_match_success = False
                    return IM_match_success
            else:
                # release mata into matb
                isen_success = self.MakeMGIsentrope(self.im1,usemgmodelbool=usemgmodelbool,impvel=vel)
                print('ISEN SUCCESS = ',isen_success)
                # find intersection between mata isentrope and matb Hugoniot
                ind = np.where((self.isen.parr > 0))[0]
                res_release = Intersection(matb.hug.uparr,matb.hug.parr,(2*self.im1.up-self.isen.uparr),self.isen.parr)
                print('res_release = ',res_release,len(res_release))
                if len(res_release[0])>0:
                    #with open('log.txt', 'a') as fp: # debugging
                    #    fp.write('CHECK res_release vel,up,P= '+str(vel)+' '+str(res_release[0][0])+' '+str(res_release[1][0]/1.e9)+'\n')
                    #    print('CHECK res_release vel,up,P= '+str(vel)+' '+str(res_release[0][0])+' '+str(res_release[1][0]/1.e9)+'\n')
                    self.im2.up=res_release[0][0] # m/s
                    self.im2.p=res_release[1][0] # Pa
                    ind = np.where(self.isen.varr > 0)[0]
                    self.im2.v=1./(np.interp(res_release[0],2*self.im1.up-self.isen.uparr[ind],1./self.isen.varr[ind])[0]) # assumes f is monotonic and increasing
                    self.im2.e=np.interp(res_release[0],2*self.im1.up-self.isen.uparr,self.isen.earr)[0] # assumes f is monotonic and increasing
                    matb.im1.up=res_release[0][0]
                    matb.im1.p=res_release[1][0]
                    matb.im1.v=1./(np.interp(res_release[0],matb.hug.uparr,1./matb.hug.varr)[0]) # assumes f is monotonic and increasing
                    matb.im1.e=np.interp(res_release[0],matb.hug.uparr,matb.hug.earr)[0] # assumes f is monotonic and increasing
                    IM_match_success = True
                    return IM_match_success
                else:
                    # something wrong with mata Mie Grueneisen model
                    #with open('log.txt', 'a') as fp: # debugging
                    #    fp.write('ERROR res_release vel mat= '+str(vel)+self.name+'\n')
                    #    print('ERROR res_release vel mat= '+str(vel)+self.name+'\n')
                    self.im2.up=0 # m/s
                    self.im2.p=0 # Pa
                    self.im2.v=0
                    self.im2.e=0
                    matb.im1.up=0
                    matb.im1.p=0
                    matb.im1.v=0
                    matb.im1.e=0
                    IM_match_success = False
                    return IM_match_success
                        
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------        
def ClStr(value):
    """ Return string with value rounded to 2 decimal places.
        Usage: ClStr(value): 
        Output: (str(round(value*100)/100.))
    """
    if np.isnan(value):
        return '0'
    else:
        return(str(round(value*100)/100.))

def ReadMaterials(matfilename='materials-data.csv'):
    """Reads in materials data CSV file and converts to mks. 
       Usage: matdata, imat = ReadMaterials(matfilename='materials-data.csv')
       Output: 2 objects: DataFrame and DF index structure.
       To view the output: display(matdata) and vars(imat)
    """
    #print('Reading materials data file and converting to mks: ',matfilename)
    matdata=pd.read_csv(matfilename) 

    # assign column indices for material property variables [so it is easier to change the materials file format]
    # these index numbers correspond to the appropriate column in the materials data CSV file
    imat = MaterialIndices()
    imat.name=0
    imat.rho0=1
    imat.c0=3
    imat.s1=5
    imat.s2=7
    imat.d=9
    imat.g0=11
    imat.q=13
    imat.uplow=15
    imat.uphigh=16
    imat.ihed=17
    imat.date=18
    imat.note=19

    # convert to MKS because everything is better in one unit system
    matdata.iloc[:,imat.rho0] *= 1000. # g/cm3 to kg/m3
    matdata.rename(columns = {'Density(g/cm3)':'Density(kg/m3)'}, inplace = True)
    matdata.iloc[:,imat.c0] *= 1000. # km/s to m/s
    matdata.rename(columns = {'c0(km/s)':'c0(m/s)'}, inplace = True)
    ind = np.where(matdata.iloc[:,imat.d] != 0)[0]
    if len(ind) > 0:
        # universal Hugoniot C=s2 [no units] D=s/km -> s/m
        matdata.iloc[ind,imat.d] /= 1000. # 1/(km/s)= s/km /1000 to s/m
        matdata.rename(columns = {'d(s/km)':'d(s/m)'}, inplace = True)
    ind = np.where(matdata.iloc[:,imat.d] == 0)[0]
    if len(ind) > 0:                        
        # line (s2=0) or quadratic Hugoniot s2=s/km -> s/m
        matdata.iloc[ind,imat.s2] /= 1000. # 1/(km/s)= s/km /1000 to s/m
        matdata.rename(columns = {'s2(s/km)':'s2(s/m)'}, inplace = True)
    matdata.rename(columns = {'s2(s/km)':'s2(s/m)'}, inplace = True) # change column name to mks
    matdata.rename(columns = {'d(s/km)':'d(s/m)'}, inplace = True) # change column name to mks
    for iii in np.where(matdata.iloc[:,imat.uplow] != -1):
        matdata.iloc[iii,imat.uplow] *= 1000 # km/s to m/s
    matdata.rename(columns = {'up_low(km/s)':'up_low(m/s)'}, inplace = True) # change column name to mks
    for iii in np.where(matdata.iloc[:,imat.uphigh] != -1):
        matdata.iloc[iii,imat.uphigh] *= 1000 # km/s to m/s
    matdata.rename(columns = {'up_high(km/s)':'up_high(m/s)'}, inplace = True) # change column name to mks

    # convert date column to string
    matdata['Date'] = matdata['Date'].astype(str)
    matdata['Notes'] = matdata['Notes'].astype(str)
    return matdata, imat # DataFrame of csv file and MaterialIndices object that defines the columns of the DF

#==================================================================================================
#### INTERSECTION FUNCTION FROM https://github.com/sukhbinder/intersection/
#### MIT LICENSE
#"""
#Give, two x,y curves this gives intersection points,
#autor: Sukhbinder
#5 April 2017
#Based on: http://uk.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
#"""
def _rect_inter_inner(x1, x2):
    n1 = x1.shape[0]-1
    n2 = x2.shape[0]-1
    X1 = np.c_[x1[:-1], x1[1:]]
    X2 = np.c_[x2[:-1], x2[1:]]
    S1 = np.tile(X1.min(axis=1), (n2, 1)).T
    S2 = np.tile(X2.max(axis=1), (n1, 1))
    S3 = np.tile(X1.max(axis=1), (n2, 1)).T
    S4 = np.tile(X2.min(axis=1), (n1, 1))
    return S1, S2, S3, S4
def _rectangle_intersection_(x1, y1, x2, y2):
    S1, S2, S3, S4 = _rect_inter_inner(x1, x2)
    S5, S6, S7, S8 = _rect_inter_inner(y1, y2)

    C1 = np.less_equal(S1, S2)
    C2 = np.greater_equal(S3, S4)
    C3 = np.less_equal(S5, S6)
    C4 = np.greater_equal(S7, S8)

    ii, jj = np.nonzero(C1 & C2 & C3 & C4)
    return ii, jj
def Intersection(x1, y1, x2, y2):
    """
    INTERSECTIONS Intersections of curves.
       Computes the (x,y) locations where two curves intersect.  The curves
       can be broken with NaNs or have vertical segments.
    usage:
    x,y=Intersection(x1,y1,x2,y2)
        Example:
        a, b = 1, 2
        phi = np.linspace(3, 10, 100)
        x1 = a*phi - b*np.sin(phi)
        y1 = a - b*np.cos(phi)
        x2=phi
        y2=np.sin(phi)+2
        x,y=intersection(x1,y1,x2,y2)
        plt.plot(x1,y1,c='r')
        plt.plot(x2,y2,c='g')
        plt.plot(x,y,'*k')
        plt.show()
    """
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)

    ii, jj = _rectangle_intersection_(x1, y1, x2, y2)
    n = len(ii)

    dxy1 = np.diff(np.c_[x1, y1], axis=0)
    dxy2 = np.diff(np.c_[x2, y2], axis=0)

    T = np.zeros((4, n))
    AA = np.zeros((4, 4, n))
    AA[0:2, 2, :] = -1
    AA[2:4, 3, :] = -1
    AA[0::2, 0, :] = dxy1[ii, :].T
    AA[1::2, 1, :] = dxy2[jj, :].T

    BB = np.zeros((4, n))
    BB[0, :] = -x1[ii].ravel()
    BB[1, :] = -x2[jj].ravel()
    BB[2, :] = -y1[ii].ravel()
    BB[3, :] = -y2[jj].ravel()

    for i in range(n):
        try:
            T[:, i] = np.linalg.solve(AA[:, :, i], BB[:, i])
        except:
            T[:, i] = np.Inf

    in_range = (T[0, :] >= 0) & (T[1, :] >= 0) & (
        T[0, :] <= 1) & (T[1, :] <= 1)

    xy0 = T[2:, in_range]
    xy0 = xy0.T
    return xy0[:, 0], xy0[:, 1]
#
### END of IM_module.py ###
#==================================================================================================
