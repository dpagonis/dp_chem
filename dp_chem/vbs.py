import math
import numpy as np
import matplotlib.pyplot as plt

# VBS_May13 : May et al 2013 JGRA VBS
# VBS_FIREX : Pagonis et al 2023 VBS

OA_MINIMUM = 1e-6 #ug/m3 ; to avoid divzero

class VBS:
  #This class is **always** at equilibirum
  #initialize with Temp (K), Pressure (mbar), and Total mass concentration of VBS (ug sm-3)
  #outputs: equilibrium OA concentrations at self.OA and self.OA_vol (ug sm-3 and ug m-3, respectively)
  #if you know OA but not total VBS conc, use self.SetOAConc(OA_std) to automatically find the total VBS conc.

  def Clausius_Clapeyron(self, Cstar298 , T , dH_298):
    Cstar = ((Cstar298*298.15)/T)*math.e**((dH_298/8.3145)*((1/298.15)-(1/T)))
    return Cstar

  def SetTP(self,T=None,P=None):                # set temperature and pressure conditions for VBS
    self.Temp = self.Temp if T is None else T
    self.Pres = self.Pres if P is None else P
    self.CStar = [self.Clausius_Clapeyron(cs , self.Temp , dH) for cs,dH in zip(self.CStar_298,self.dHvap)]
    self.PartitionVBS()

  def CalcVolConcs(self):
    self.StdToVol = (273.15/self.Temp) * (self.Pres/1013.25) #Std * StdToVol = Vol ; ATTN: this is flipped from the Jimenez group convention
    self.VBSConc_vol = self.VBSConc * self.StdToVol
    self.VBSMass_vol = self.VBSConc_vol * self.Fi
    self.OA_vol = self.OA * self.StdToVol
    #standard concentration floats: VBSConc, OA
    #volumetric conctration floats: VBSConc_vol, OA_vol
    #volumetric concntration array: VBSMass_vol

  def SetVBSConc(self,VBSConc_std): #directly set the total VBS concentration (ug sm-3)
    self.VBSConc = VBSConc_std
    self.PartitionVBS()

  def SetOAConc(self,OA_std):   #set the total VBS concentration by inputting the equilibrium-partitioned OA standard concentration (ug sm-3)
    while abs(self.OA-OA_std)/OA_std > 1e-6:
      self.VBSConc = self.VBSConc * OA_std / self.OA
      self.PartitionVBS()

  def SetOAConc_vol(self,OA_vol):   #set the total VBS concentration by inputting the equilibrium-partitioned OA concentration (ug m-3)
    self.SetOAConc(OA_vol / self.StdToVol)

  def ScaleConc(self,ScaleFactor):  #scale the total VBS concentration. e.g. ScaleFactor = 0.9 represents 10% dilution
    self.VBSConc = self.VBSConc * ScaleFactor
    self.PartitionVBS()

  def PartitionVBS(self):
    self.CalcVolConcs()
    OACalc = max(self.OA_vol,OA_MINIMUM)
    dOA = 1
    while dOA>1e-6:
      Fp = [(1 + (x/OACalc))**-1 for x in self.CStar]
      Mass_p = self.VBSMass_vol * Fp
      LastCOA = OACalc
      OACalc = np.sum(Mass_p)
      dOA = (abs(LastCOA - OACalc)/OACalc)
    self.Fp = Fp
    self.OA_vol = OACalc
    self.OA = self.OA_vol / self.StdToVol

  def plot(self):
      plt.bar(range(len(self.CStar_298)),self.VBSMass_vol,tick_label=[str(cs) for cs in self.CStar_298])
      OAbars = self.Fp * self.VBSMass_vol
      plt.bar(range(len(self.CStar_298)),OAbars)
      plt.xlabel('C*298 (ug m-3)')
      plt.ylabel('Concentration (ug m-3)')
      plt.figtext(0.2,0.9,'VBS = %s\nVBS = %.02f ug sm-3\nOA = %.02f ug sm-3' % (self.Name,self.VBSConc, self.OA))
      plt.figtext(0.6,0.9,'T = %.01f K\nP = %.0f mbar' % (self.Temp,self.Pres))
      plt.show()

  def __init__(self,Temp=None,Pres=None,VBSConc=None):
    self.Temp = 298.15 if Temp is None else Temp
    self.Pres = 1013.25 if Pres is None else Pres
    self.VBSConc = 100 if VBSConc is None else VBSConc #ugsm3
    self.CStar = [self.Clausius_Clapeyron(cs , self.Temp , dH) for cs,dH in zip(self.CStar_298,self.dHvap)]
    self.OA = self.VBSConc / 2
    self.PartitionVBS()
#end of class VBS

class VBS_May13(VBS):
  #VBS defined in May et al JGRA 2013 (https://doi.org/10.1002/jgrd.50828)
  Name = "May_2013_JGRA"
  Fi = np.array([0.2,0,0.1,0.1,0.2,0.1,0.3])
  CStar_298 = np.array([10**(i-2) for i in range(7)])
  dHvap = np.array([(85-4*math.log10(x))*1000 for x in CStar_298])
#end of class VBS_May13

class VBS_FIREX(VBS):
    #VBS fitted to FIREX-AQ thermal denuder data
    Name = "FIREX-AQ"
    Fi = np.array([0.2,0.1,0.2,0.1,0,0.2,0.2])
    CStar_298 = np.array([0.01,0.1,1,10,100,1000,10000])
    dHvap = np.array([(85-4*math.log10(x))*1000 for x in CStar_298])
  #end of class VBS_FIREX


