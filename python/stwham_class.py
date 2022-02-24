import numpy as np
import pickle
#from walker import walkers

class st_wham:
    def __init__(self,hist,Ghist,hcounts,minval,maxval,lambdas,eta,Ho,bl):
        self.kb = 0.0019872041
        self.nbins, self.nwalkers = np.shape(hist)[0],np.shape(hist)[1]
        self.checklimit=1e-4
        self.lambdas = lambdas
        # PDF Calculations
        self.hist  = hist
        self.pdf = np.zeros((self.nbins,self.nwalkers))
        self.count = np.zeros(self.nwalkers)
        self.totpdf = np.zeros(self.nbins)
        # Gfactor Calculations
        self.hcounts = hcounts
        self.Ghist = Ghist
        self.gfactor= np.zeros((self.nbins,self.nwalkers))
        self.totgfactor=np.zeros(self.nbins)
        self.totcount=0
        self.binsize=(maxval-minval)/float(self.nbins)
        # Stuff for the calculation
        self.TH, self.S = np.zeros(self.nbins), np.zeros(self.nbins)
        self.betaH, self.betaW = np.zeros(self.nbins),np.zeros(self.nbins)
        self.A = np.zeros(self.nbins)
        self.hfrac = np.zeros((self.nbins,self.nwalkers))
        self.bH_alpha = np.zeros((self.nbins,self.nwalkers))
        self.bH = np.zeros(self.nbins)
        self.dln_A_alpha = np.zeros((self.nbins,self.nwalkers))
        self.dln_A = np.zeros(self.nbins)
        self.f_alpha = np.zeros((self.nbins,self.nwalkers))
        # Limits
        self.bstart, self.bstop = None,None
        self.minval, self.maxval = minval,maxval
        self.eta,self.Ho = eta,Ho
        self.bl = bl
        # Prints
        if self.bl == -1: self._PrintParams()
        self._runSTWHAM()
        return
    def _EffTemp(self,shift,H):
        Teff = shift + self.eta*(H-self.Ho)
        return Teff
    def _PrintParams(self):
        print("ST WHAM Parameters")
        print("binsize = %10.5f" % self.binsize)
        print("minval  = %10.5f" % self.minval)
        print("maxval  = %10.5f" % self.maxval)
        print("nbins   = %10.5f" % self.nbins)
    def _runSTWHAM(self):
        self._calcPDF()
        self._calcGfactor()
        self._setLimits()
        self._calcWHAM()
        if self.bl == -1: self._writeWHAM()
    def _calcPDF(self):
        """
        This function calculates the individual pdf and the total pdf
        as well as the counts
        """
        for rep in range(self.nwalkers):
            self.count[rep]=0
            for i in range(self.nbins):
                self.pdf[i][rep] = self.hist[i][rep]  #pdf of each replica
                self.count[rep] = self.count[rep] + self.pdf[i][rep] #count
                self.totpdf[i] = self.totpdf[i]+self.pdf[i][rep] # adds to total pdf
            self.pdf[:,rep] = self.pdf[:,rep]/self.count[rep] # Normalized replica pdf
            self.totcount = self.totcount + self.count[rep] # Total number of data points
        self.totpdf = self.totpdf/self.totcount # Normalized total PDF
        return
    def _calcGfactor(self):
        """
        This calculates the g_factor from my derivation for the derivative
        """
        for rep in range(self.nwalkers):
            for i in range(self.nbins):
                self.totgfactor[i] = self.totgfactor[i] + self.Ghist[i][rep]
        self.totgfactor = self.totgfactor/np.sum(self.hcounts)
        
        for rep in range(self.nwalkers):
            for i in range(self.nbins):
                if self.totgfactor[i] != 0:
                    self.gfactor[i][rep] = self.Ghist[i][rep]/(self.totgfactor[i]*np.sum(self.hcounts,axis=0)[rep])
                else:
                    self.gfactor[i][rep]=0
        return
    def _setLimits(self):
        for i in range(self.nbins):
            if (self.totpdf[i] > self.checklimit):
                self.bstart = i + 3
                break
        if (self.bstart == None): exit("Error: Lower end of logarithm will be undefined")
        for i in range(self.nbins-1,self.bstart, -1):
            if (self.totpdf[i] > self.checklimit):
                self.bstop = i - 3
                break
        if (self.bstop == None): exit("Error: Upper end of logarithm will be undefined")
    def _interpolate(self,val,i,rep):
        if rep == -1:
            return np.log(val[i+1]/val[i-1])/(2*self.binsize)*self.kb
        else:
            return np.log(val[i+1][rep]/val[i-1][rep])/(2*self.binsize)*self.kb
    def _checkLimit(self,quant,i,rep):
        if rep == -1:
            return quant[i+1] > self.checklimit and quant[i-1] > self.checklimit
        else:
            return quant[i+1][rep] > self.checklimit and quant[i-1][rep] > self.checklimit

    def _calcWHAM(self):
        for i in range(self.bstart,self.bstop):
            if self._checkLimit(self.totpdf,i,-1): # Check it is above the checklimit
                self.betaH[i] = self._interpolate(self.totpdf,i,-1)
            else: # Otherwise 0
                betaH[i]=0

            for rep in range(self.nwalkers):
                if self._checkLimit(self.pdf,i,rep):
                    self.bH_alpha[i][rep] = self._interpolate(self.pdf,i,rep)
                else:
                    self.bH_alpha[i][rep] = 0
                if self._checkLimit(self.Ghist,i,rep):
                    self.dln_A_alpha[i][rep] = self._interpolate(self.Ghist,i,rep)
                if self.totpdf[i] > 0:
                    self.f_alpha[i][rep]=(self.count[rep]*self.pdf[i][rep])/(self.totcount*self.totpdf[i])
                    H      = self.minval + (i*self.binsize+self.binsize/2)
                    weight = np.nan_to_num(1/self._EffTemp(self.lambdas[rep],H))
                    self.betaW[i] = self.betaW[i] + self.f_alpha[i][rep]*weight
                    self.bH[i] = self.bH[i] + self.f_alpha[i][rep]*self.bH_alpha[i][rep]
                    self.A[i] = self.A[i] +  self.Ghist[i][rep]/np.sum(self.hcounts,axis=1)[i]
                    self.dln_A[i] = self.dln_A[i] + self.f_alpha[i][rep]*(self.gfactor[i][rep]-1)*(self.bH_alpha[i][rep] + weight) + self.f_alpha[i][rep]*self.gfactor[i][rep]*self.dln_A_alpha[i][rep]
            self.TH[i] = 1 / (self.bH[i] + self.betaW[i])
    def _writeWHAM(self):
        #with open("TandS_A.out",'w') as tout:
        #    for i in range(self.bstart,self.bstop):
        #        tout.write("%f %f %f %f\n" % (self.minval+(i*self.binsize), self.TH[i],self.A[i],self.dln_A[i]))
        self.xvals = self.minval + np.arange(self.bstart,self.bstop)*self.binsize
        np.savetxt("TandS_A.out", np.c_[self.xvals, self.TH[self.bstart:self.bstop],self.A[self.bstart:self.bstop],self.dln_A[self.bstart:self.bstop]])

def calc_halffrac(H,T,A,minval,maxval):
    TH_spline = CubicSpline(H,T)
    A_spline  = CubicSpline(H,A)
    xs = np.arange(-1930,-1360,0.001)
    close = abs(A_spline(xs[0]) - 0.5)
    hval = xs[0]
    for x in xs:
        dval = abs(A_spline(x) - 0.5)
        if dval < close:
            close = dval
            hval = x
    return TH_spline(hval)



if __name__ == "__main__":
    import argparse
    from scipy import stats
    from scipy.interpolate import CubicSpline
    parser=argparse.ArgumentParser()
    parser.add_argument('-nbins', default=200, type=int, help="Number of bins")
    parser.add_argument('-nreps', default=76, type=int, help="Number of replicas")
    parser.add_argument('-nblocks', default=5, type=int, help='Number of blocks')
    parser.add_argument('-prefix', default="all_dumps/DBPC_xGel_walker-", type=str, help="File prefix for ensemble info")
    args = parser.parse_args()

    nbins = args.nbins
    nblocks = args.nblocks
    nreps = args.nreps
    prefix = args.prefix

    t_val = stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

    # Read in the allwalkers pckl file
    allwalkers=pickle.load(open('allwalkers.pckl','rb'),encoding='latin1')
    allwalkers.post_process(nbins,nblocks)

    # Store Hdata
    Hdata = allwalkers.repdata["PotEng"]
    xGel,xGel_hist = [],[]

    # Histogram
    hcounts = []
    for r in range(nreps):
        xGel.append(np.genfromtxt(prefix+"%d.out"%r,usecols=0)[:10000])
        test, bin2 = np.histogram(Hdata[r][::10],bins=nbins,range=(allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"]),weights=xGel[r])
        hc, bin2 = np.histogram(Hdata[r][::10],bins=nbins,range=(allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"]))
        xGel_hist.append(test)
        hcounts.append(hc)
    Ghist = np.array(xGel_hist).T
    hcounts = np.array(hcounts).T

    # Call STWHAM + Run It
    TandS = st_wham(allwalkers.hist["PotEng"],Ghist,hcounts,allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"],allwalkers.lambdas,allwalkers.eta,allwalkers.Ho,-1)
    
    bl_data = []
    bl_xgel = []
    for r in range(nreps):
        bl_data.append(np.split(Hdata[r][::10],nblocks))
        bl_xgel.append(np.split(xGel[r],nblocks))
    bl_data=np.transpose(np.array(bl_data),(1,0,2))
    bl_xgel=np.transpose(np.array(bl_xgel),(1,0,2))

    bl_Ahist = []
    bl_hcounts = []
    for b in range(nblocks):
        bl_Ahist.append([])
        bl_hcounts.append([])
        for r in range(nreps):
            test,bin2 = np.histogram(bl_data[b][r],bins=nbins,range=(allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"]),weights=bl_xgel[b][r])
            hc, bin2 = np.histogram(bl_data[b][r],bins=nbins,range=(allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"]))
            bl_Ahist[b].append(test)
            bl_hcounts[b].append(hc)
        bl_Ahist[b] = np.array(bl_Ahist[b]).T
        bl_hcounts[b] = np.array(bl_hcounts[b]).T
    bl_TandS = []
    bl_TH = []
    bl_A = []
    for b in range(nblocks):
        bl_TandS.append(st_wham(allwalkers.bl_hist["PotEng"][b],bl_Ahist[b],bl_hcounts[b],allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"],allwalkers.lambdas,allwalkers.eta,allwalkers.Ho,b))
        bl_TH.append(bl_TandS[b].TH)
        bl_A.append(bl_TandS[b].A)
    errTH=np.std(bl_TH,axis=0)*t_val
    errA =np.std(bl_A,axis=0)*t_val
    bstart,bstop = TandS.bstart, TandS.bstop
    np.savetxt("TandS_final_A.out",np.c_[TandS.xvals,TandS.TH[bstart:bstop],errTH[bstart:bstop],TandS.A[bstart:bstop],errA[bstart:bstop]])
    Tm_est=calc_halffrac(TandS.xvals,TandS.TH[bstart:bstop],TandS.A[bstart:bstop],allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"])
    bl_Tm_est = []
    for b in range(nblocks):
        bl_Tm_est.append(calc_halffrac(TandS.xvals,bl_TH[b][bstart:bstop],bl_A[b][bstart:bstop],allwalkers.minval["PotEng"],allwalkers.maxval["PotEng"]))
    print("ML Tm Etimate %10.5f +/- %10.5f" % (Tm_est,np.std(bl_Tm_est)*t_val))
        
