#!/usr/bin/env python2
# @Copyright 2007 Kristjan Haule
# @Modified by Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2025

import sys, re, os, glob, shutil, types
from scipy import *
from scipy import interpolate
from utils import DmftEnvironment
from cmpEimp2 import GiveImpFunctional
import brod
import subprocess
import numpy
import dmft_DC

nv = map(int,numpy.__version__.split('.'))
if (nv[0],nv[1]) < (1,6):
    loadtxt = io.read_array
    savetxt = io.write_array

def add_nf(self, params):
    first_10_lines = []  # Define the list before using it
    sum_value = None  # It's also a good practice to initialize variables like sum_value

    with open(self.dir + 'nmat.txt', 'r') as file:
        for i, line in enumerate(file):
            if i < 16:  # Only read the first 16 lines
                first_10_lines.append(line)
                if line.strip().startswith('sum'):
                    sum_value = float(line.split()[1].strip())

    # Step 2: Add nf to the beginning of Gf.out and Sig.out
    for filename in [self.dir + 'Gf.out', self.dir + 'Sig.out']:
        with open(filename, 'r') as file:
            data = file.read()
        with open(filename, 'w') as file:
            file.write("# nf=" + str(sum_value)  + '\n' + data)


def loadtxtE(filename):
    if (nv[0],nv[1]) < (1,6):
        return io.read_array(filename)[:,:-2]
    else:
        return loadtxt(filename,ndmin=2)

def mprint(Us):
    Qcomplex =  type(Us[0,0])==complex or type(Us[0,0])==complex128
    for i in range(shape(Us)[0]):
        for j in range(shape(Us)[1]):
            if Qcomplex:
                print "%11.8f %11.8f  " % (real(Us[i,j]), imag(Us[i,j])),
            else:
                print "%11.8f  " % Us[i,j],
        print
    
Default_CoulombF = 'Ising'
            
class IMP_NORG(object):
    """Wrapper for norg
    """
    def __init__(self, Znuc, dire, params, Sigind, CF):
        # Environment
        self.env = DmftEnvironment()
        self.mpi_prefix = self.env.MPI
        
        # impurity directory
        self.dir = dire+'/'

        # create impurity directory if it does not yet exist
        if len(glob.glob(dire))==0 : os.mkdir( dire )
        
        # log-file is opened
        self.fh_info = open(self.dir + 'info.out', 'w')

        if params.has_key('CoulombF'):
            CoulombF = params['CoulombF'][0].strip("'").strip('"')
        else:
            params['CoulombF']=["'"+Default_CoulombF+"'", "# Full Coulomb repulsion"]
            CoulombF = Default_CoulombF
        
        if params.has_key('icase'):
            self.icase = params['icase']
        if params.has_key('UlamJ'):
            self.UlamJ = params['UlamJ']
        
        self.Sigind = array(Sigind)
        self.CF = CF
        
        # find smallest sigind to shift siginds so they start from 1
        minsigind = min([s for s in self.Sigind.flat if s != 0])
        self.Sigind_shifted = where(self.Sigind != 0, self.Sigind-minsigind+1, 0)
        
        # Transformation between spheric harmonics and local coordinates needs to be
        # used in exact diagonalization of the atom.
        fh_T = open(self.dir + 'Trans.dat', 'w')
        print >> fh_T, len(Sigind), len(Sigind[0]), '#  size of Sigind and CF'
        print >> fh_T, '#---- Sigind follows'
        for row in self.Sigind_shifted:
            for elem in row:
                print >> fh_T, ("%3s " % elem),
            print >> fh_T
        print >> fh_T, '#---- CF follows'

        for row in CF:
            for elem in row:
                print >> fh_T, "%12.8f %12.8f " % (elem.real, elem.imag),
            print >> fh_T

        # parameters, which will be given to norg and will be written to PARAMS
        _nom_ = 100
        if params.has_key('beta'): _nom_ = int(3*params['beta'][0])
        self.sparams={'exe':'norg', 'Sig':'Sig.out', 'Delta':'Delta.inp', 'cix':'actqmc.cix', 'nom':_nom_, 'mode':'GH', 'aom':1, 'Ncout':1000000, 'GlobalFlip':1000000, 'Gf':'Gf.out', 'OffDiagonal':'real'}

        # If old distribution of kinks exists, it is copied to imp directory
        status = glob.glob("status.*")
        if self.dir!='./':
            #print 'self.dir=', self.dir
            for fil in status:
                shutil.copy2(fil, self.dir+'/'+fil)
        

        self.PARAMS = 'PARAMS.norg'
        
        # file from database which contains atomic information for this particular atom
        database_arg = self.env.ROOT + '/database/' + 'atom.'+str(Znuc)+'.py'

        if os.path.exists(database_arg):
            execfile(database_arg)
            database = locals() # it seems python does not know about variables obtained by execfile
        

        # The following variables will be set by the following procedure:
        #   if the variable is defined in params, it takes its value from params otherwise it checks its definition in database
        vars = ['l', 'para', 'qOCA', 'kOCA', 'mOCA', 'Ex', 'Ep', 'J', 'cx', 'Eoca', 'qatom', 'n', 'Ekeep', 'Ekeepc', 'Nmaxc', 'Ewindow', 'max_M_size', 'CoulombF', 'add_occupancy', 'mode']
        
        self.spr={}
        for var in vars:
            if params.has_key(var):
                self.spr[var] = params[var][0]
            elif database.has_key(var):
                self.spr[var] = database[var]
            if params.has_key('atom_'+var):
                self.spr[var] = params['atom_'+var][0]
        
        self.spr['OCA_G']=False # no need for oca diagrams here
        # Some variables have different name in database and input file
        # spr['n'] -- params['nc']  -- database['n']
        if params.has_key('nc'):
            self.spr['n'] = params['nc'][0]
        elif database.has_key('nc'):
            self.spr['n'] = database['nc']
        
        if params.has_key('Ncentral'):
            self.spr['Ncentral'] = params['Ncentral'][0]
        elif self.spr.has_key('n'):
            self.spr['Ncentral'] = self.spr['n'][1:len(self.spr['n'])-1]

        # UlamJ= {0: (U, lmbda, epsilon, [J2, J4])}
        if params.has_key('UlamJ') and self.UlamJ[1]>1e-5:  # DCs='exact'
            # If we use Yukawa screening, we need to change Jhunds to be compatible with Yukawa form of screening
            # For dielectric screening we do not change Jhunds from user specified values
            Jh = self.UlamJ[3]
            if self.spr['l']==3: Jh = sum(Jh)/len(Jh)  # Currently does not handle an array of J's. You should correct that.
            self.spr['J']=Jh

        # stores some of just defined variables
        self.J = self.spr['J']
        if self.spr.has_key('l'):
            self.l = self.spr['l']
        else:
            if (Znuc>=57 and Znuc<=71) or (Znuc>=89 and Znuc<=103): self.l=3  # lanthanides & actinied
            else: self.l=2
        # extra off-diagonal spin-orbit
        if self.spr.has_key('cx'):
            self.cx = self.spr['cx']
        else:
            self.cx=0.0
    

    def _exe(self, params, DCs, extn, UpdateAtom, gbroad=0.0, kbroad=0.0, maxAc=200.):
        """ Executes the CTQMC impurity solver
        """
        REMOVE_STATUS=False
        Qaverage = False
        # Reading impurity levels
        (Eimp,Olap,Edc) = loadtxtE(self.dir+'/Eimp.inp')
        # Reading Delta
        Delta_dat = loadtxt(self.dir+'/Delta.imp').transpose()
        om = Delta_dat[0]
        Delta = Delta_dat[1::2] + Delta_dat[2::2]*1j
        # subtracting Edc, because the ksum used (Hk+s_oo-Edc)
        self.Eimps = Eimp-Edc
        #
        ##! add for the norg needs. test for the t2g system.
        # norg_eimps = ones(3) * (Eimp-Edc)
        # norg_degs  = self.Sigind.diagonal()[0:3]
        # norg_eimps = ones(len(self.Eimps)) * (Eimp-Edc)        #! test for the t2g and eg system.
        # norg_degs  = self.Sigind.diagonal()[0:len(self.Eimps)] #! test for the t2g and eg system.
        norg_degs  = [self.Sigind[i, i] for i in range(len(self.Sigind.diagonal())) if self.Sigind[i, i] != 0]    #! test for the t2g, eg and full-d system.
        # norg_eimps = ones(len(norg_degs)) * (Eimp-Edc)                                                            #! test for the t2g, eg and full-d system.
        norg_eimps = ones(len(self.Eimps)) * (Eimp-Edc)                                                           #! test for deg-d system.
        params['Ed'] = [norg_eimps, "# Impurity levels"]
        params['deg_idx'] = [norg_degs, "# Impurity leve degs"]
        # params['restrain'] = [numpy.array([0, -2, -3, 0, 3, 2]), "# norg restraion"]
        # params['restrain'] = [numpy.array([0, -1, -1, 0, 1, 1]), "# norg restraion"] ##! for testing
        # params['distribute'] = [numpy.array([1, 3, 0, 1, 0, 3]), "# norg distribute"]
        #
        if DCs=='fixn0':
            self.Eimps = ones(len(self.Eimps))*sum(self.Eimps)/len(self.Eimps)  # All impurity levels should be the same
            print >> self.fh_info, 'Eimps(fixn0)=', Eimps
	    
        # params['Ed'] = [self.Eimps.tolist(), "# Impurity levels"]

        
        # Our choice for mu_QMC is the following:
        #   impurity level for 5/2 is E_{5/2}=E0-2*cx and for 7/2 is E_{7/2}=E0+3/2*cx
        #   Here we correct E_{5/2} for 2*cx, which is already counted in cix file
        #E0 = params['Ed'][0][0]+2*self.cx

        if (self.l == 3 and Qaverage):
            # 5/2 energy levels:
            Eimps_5_2 = [self.Eimps[self.Sigind[i,i]-1] for i in range(2*self.l)]
            aEimps_5_2 = sum(Eimps_5_2)/(2*self.l)
            print '<E_5/2>=', aEimps_5_2
            #E0 = Eimps[0]+2*self.cx
            E0 = aEimps_5_2 + 2*self.cx
        else:
            E0 = self.Eimps[0]+2*self.cx
            
        mu_QMC = -E0

        Ident = zeros(len(self.Eimps))
        for i in range(len(self.Sigind)):
            if self.Sigind_shifted[i,i]>0: Ident[self.Sigind_shifted[i,i]-1]=1.0
        #print Ident
        '''
        # Below exact diagonalization of the atom is run to produce 'cix.out' file
        if UpdateAtom or len(glob.glob(self.dir+'/'+self.sparams['cix']))==0 : # no cix file yet -> create it

            #if type(Eimps)!= list :  Eimps = Eimps.tolist()
            #print 'Eimps=', Eimps, 'E0=', E0
            Eimps = (array(self.Eimps)-E0*Ident).tolist()
            print 'Eimps_new=', Eimps, ' muQMC=', mu_QMC
            Impq = [self.Sigind_shifted[i,i]-1 for i in range(len(self.Sigind_shifted))]
            print 'Impq=', Impq
            
            # Variables which are updated from the current paramateres
            for var in ['CoulombF','mode']:
                if params.has_key(var):
                    self.spr[var] = params[var][0]
                    
            if (self.l == 3) :
                # executable for ED of the atom
                self.atom_exe = self.env.ROOT + '/atomh'
                # argument list is created to later run the executable for exact diagonalization of the atom
                atom_arg = '-ns %d -ne %d ' % (self.spr['n'][0], self.spr['n'][-1])
                if (not Qaverage):
                    atom_arg += ' -Eimp \"'+str(Eimps)+'\" '
                    atom_arg += ' -Impq \"'+str(Impq)+'\" '
                    
                for k,p in self.spr.items():
                    if k=='mode':
                        if p[0]=='S':
                            atom_arg += '-HB2 1 '
                        continue
                    if type(p)==types.ListType:
                        atom_arg += '-'+k+' "'+str(p)+'" '
                    else:
                        atom_arg += '-'+k+' '+str(p)+' '
                # print to info file
                print >> self.fh_info, 'atom_arg=', atom_arg
                # runs ED for the atom
                cmd = 'cd '+ self.dir + '; ' + self.atom_exe +' ' +atom_arg+' > nohup.out 2>&1'
                # subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)

            elif (self.l < 3) :
                # executable for ED of the atom
                self.atom_exe = self.env.ROOT + '/atom_d.py'

                if params.has_key('U') and params.has_key('Ud') and params['U'][0]==0.0 and params['U'][0]!=params['Ud'][0]:
                    # Must be cluster-DMFT case. We will assume it is 2-site cluster case (link)
                    self.FromSingleImpurity_to_Link(Eimps, params)
                else:
                    # Old good single impurity calculation
                    atom_arg = ''
                    for k,p in self.spr.items():
                        if k=='mode':
                            if p[0]=='S':
                                atom_arg += '\"HB2=True\" '
                            continue
                        atom_arg += '\"'+k+'='+str(p)+'\" '
                    atom_arg += ' \"Eimp='+str(Eimps)+'\" '
                    
                    print >> self.fh_info, 'atom_arg=', atom_arg
                    print atom_arg
                    
                    # runs ED for the atom
                    cmd = 'cd '+ self.dir + '; ' + self.atom_exe +' ' +atom_arg+' > nohup.out 2>&1'
                    print '.... running ..', cmd
                    # subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)

                # If new cix file is created, the eigenvectors are reordered and hence
                # previous ctqmc kinks can not be directly used. One would need to reorder
                # them. But this is not very simple. Rather remove the status files.
                if (REMOVE_STATUS):
                    for filename in glob.glob(self.dir+'/status.*') :
                        os.remove( filename )
            else:
                print >> self.fh_info, 'l=', self.l, ' not yet implemented!'
        
        self.fh_info.flush()
        '''
        
        if self.l==3:
            params['exe']=['ctqmcf', '#']

        params['mu'] = [mu_QMC, '# QMC chemical potential']
        
        # executable
        #shutil.copy2(self.env.ROOT+'/'+params['exe'][0], self.dir+'/'+params['exe'][0])
        #self.q_exe = './'+params['exe'][0]
        self.q_exe = self.env.ROOT+'/'+params['exe'][0]
        '''
        # Preparing params.dat file
        fp = open(self.dir+self.PARAMS, 'w')
        kparams = params.keys()
        for p in self.sparams.keys():
            if p not in kparams:  # if this parameter is overwritten by the input parameters, skip it
                fp.write(p + '  ' + str(self.sparams[p]) + '\n')
            
        for p in params.keys():
            print >> fp, "%s  %s " % (p, params[p][0]), '\t ', params[p][1]
        fp.close()
        '''
        # Preparing params.dat file
        fp = open(self.dir+self.PARAMS, 'w')
        for p in ['mode', 'Ed', 'deg_idx', 'U', 'J', 'CoulombF', 'beta', 'fit_range', 'fit_points', 'fit_nbaths', 'NOOC', 'restrain', 'distribute', 'pred_gs_deg', 'ful_pcl_sch', 'weight_nooc', 'weight_freze']:
            if p in params:
                print >> fp, "%s  %s " % (p, params[p][0]), '\t ', params[p][1]
        fp.close()

        Sigind = array(self.Sigind_shifted)
        Diagonal={}
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0: Diagonal[Sigind[i,j]-1]= (i==j)

        print 'Diagonal=', Diagonal
        

        # Preparing Delta.inp file
        fp = open(self.dir+self.sparams['Delta'], 'w')
        # Checking that Delta is causal and not regular
        for ib in range(shape(Delta)[0]):
            if Diagonal[ib]:
                for im in range(shape(Delta)[1]):
                    if abs(Delta[ib,im])>maxAc*pi : Delta[ib,im] = -1e-10j
                    if Delta[ib,im].imag>0 : Delta[ib,im] = Delta[ib,im].real-1e-10j
        

        # creating big mesh
        T=1./params['beta'][0]
        m_max = int(om[-1]/(pi*T)+1e-6)
        lom = arange(1,m_max+1,2)*pi*T
        # interpolate on big mesh
        lDelta=[]

        print 'shape(Delta)=', shape(Delta)
        print 'shape(om)=', shape(om)
        print 'shape(lom)=', shape(lom)
        print 'ss=', shape(Delta[0].real)
        print 'om[0]=', om[0], 'lom[0]=', lom[0], 'om[-1]=', om[-1], 'lom[-1]=', lom[-1]
        
        for ib in range(len(Delta)):
            fDr = interpolate.UnivariateSpline(om, Delta[ib].real, s=0)
            fDi = interpolate.UnivariateSpline(om, Delta[ib].imag, s=0)
            lDelta.append( fDr(lom) + fDi(lom)*1j )
            
        for im,ome in enumerate(lom):
            print >> fp, ome,
            for b in range(len(lDelta)):
                print >> fp, "%19.12f  %19.12f " % (lDelta[b][im].real, lDelta[b][im].imag),
            print >> fp
        fp.close()
        
        shutil.copy2(self.dir+self.sparams['Delta'], self.dir+self.sparams['Delta']+'.'+extn)
        '''
        # saving three previous steps of status files, in case we need to restart from previous steps.
        if os.path.exists(self.dir+'status_2.tgz'):
            shutil.move(self.dir+'status_2.tgz', self.dir+'status_3.tgz')
        if os.path.exists(self.dir+'status_1.tgz'):
            shutil.move(self.dir+'status_1.tgz', self.dir+'status_2.tgz')
        if os.path.exists(self.dir+'status.000'):
            cmd = 'cd '+self.dir+'; tar czvf status_1.tgz status.*'
            # subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
        '''

        # Below we execute norg
        cmd = 'cd '+self.dir+'; '+self.mpi_prefix+' '+self.q_exe+' > nohup_imp.out 2>&1 '
        print >> self.fh_info, cmd
        print >> self.fh_info, 'Running ---- norg -----'
        subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
	add_nf(self, params)


    def HighFrequency(self, params, DCs, nl_imp, extn, wbroad=0.0, kbroad=0.0):
        """ Reads impurity output and prepares for further execution
          Output:
               Sig[b][iom]  -- DMFT dynamic self-energy which vanishes at infinity
               sinf[b]      -- DMFT self-energy at infinity
               Edc[b]       -- double counting using DMFT occupancies
        """


        if kbroad>0.0 or wbroad>0.0:
            bexe = self.env.ROOT + '/broad'
            shutil.move(self.dir+self.sparams['Sig'], self.dir+self.sparams['Sig']+'b')
            broad_cmd = 'cd '+self.dir+'; '+bexe+' -w '+str(wbroad)+' -k '+str(kbroad)+' '+self.sparams['Sig']+'b  > '+self.sparams['Sig']
            print broad_cmd
            subprocess.call(broad_cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)

        
        OffDiagonalExist = False
        Sigind = array(self.Sigind_shifted)
        self.Diagonal={}
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0:
                    self.Diagonal[Sigind[i,j]-1]= (i==j)
                    if i!=j: OffDiagonalExist = True


        # Reading Self-energy of the impurity
        Sg_data = loadtxt(self.dir + self.sparams['Sig']).transpose()
        om = Sg_data[0]
        Sigo = Sg_data[1::2,:]+Sg_data[2::2,:]*1j
        
        '''
        # Computing self-energy in alternative way to keep off-diagonal components of self-energy forcing G_{loc}=G_{imp} exactly
        if (self.sparams['OffDiagonal'] and OffDiagonalExist):
            Sign = self.SigmaOffDiagonal(om,Sigo,wbroad,kbroad)
            if Sign is not None: Sigo = Sign
        '''

        # Reading Sig.out which was produced by norg
        fS = open(self.dir + self.sparams['Sig'],'r')
        # Sig.out contains nf, moment, ntot,...
        first_line = fS.next().strip()
        fS.close()
        adat = first_line.split()
        # mom=0
        TrSigmaG=0
        '''
        for par in adat: # setting nf from ctqmc output
            m = re.search('nf=', par) or re.search('TrSigmaG=', par) or re.search('mom=', par)
            if m is not None : exec(par)
        ntot = nf
        '''
        #!##################################################################### 
        #! get the ntot straight from the NORG outputs modify by jmw@ruc.edu.cn
        mom = []
        nf = None  # It's also a good practice to initialize variables like nf
        # (Eimp_old, Olap_old, Edc_old) = loadtxtE(self.dir+'/Eimp.inp')
        total_lines_numbaers = len(params['deg_idx'][0])*2 + 4
        
        with open(self.dir + 'nmat.txt', 'r') as file:
            lines = file.readlines()
                # Skip the first line (header)
        lines = lines[1:]

        i = 0
        while i < total_lines_numbaers:
            line = lines[i].strip()
            if line.startswith('sup') or line.startswith('sdn') or line.startswith('sum'):
                break
            else:
                # Read n_up
                n_up_value = float(line.split()[1])

                i += 1
                if i < total_lines_numbaers:
                    n_dn_line = lines[i].strip()
                    n_dn_value = float(n_dn_line.split()[1])

                    # Compute mom for this orbital
                    mom_value = n_up_value + n_dn_value
                    mom.append(mom_value)
                else:
                    # Handle error if unexpected end of file
                    raise ValueError("Unexpected end of file when reading n_dn value.")
            i += 1  # Move to the next line


        #! The CODE is under testing jmw@ruc.edu.cn
        # Merge the mom array based on the deg_idx array
        deg_idx = params['deg_idx'][0]  # Example deg_idx array
        Max_Deg = max(deg_idx)  # Get the maximum degeneracy value
        filtered_mom = [0] * Max_Deg  # Initialize a new array

        for deg, m in zip(deg_idx, mom):
            filtered_mom[deg - 1] += m  # Accumulate elements with the same degeneracy

        mom = filtered_mom  # Update the mom array
        #!#####################################################################
        # Process the remaining lines to get ntot
        while i < total_lines_numbaers:
            line = lines[i].strip()
            if line.startswith('sum'):
                nf = float(line.split()[1])
                break
            i += 1

        ntot = nf
        #! get the ntot straight from the NORG outputs modify by jmw@ruc.edu.cn
        #

        for b in range(shape(Sigo)[0]):   #print b, 'This is diagonal bath, correcting it'
            for iom in range(shape(Sigo)[1]):
                if self.Diagonal[b] and Sigo[b,iom].imag>0:
                    Sigo[b,iom] = Sigo[b,iom].real
                    print 'Setting something to zero'

        sinf = copy(Sigo[:,-1]).real

        print 'sinf=', sinf
        

        # Reading old Edc and impurity levels
        (Eimp_old,Olap_old,Edc_old) = loadtxtE(self.dir+'/Eimp.inp')
        
        print 'Diagonal=', self.Diagonal
        print 'len(Edc_old)=', len(Edc_old)
        
        if params.has_key('CoulombF'):
            CoulombF = params['CoulombF'][0].strip("'").strip('"')
        else:
            CoulombF = Default_CoulombF
        
        # Here we compute self-energy(infinity) and 
        # double-counting using impurity occupancy
        print >>self.fh_info, 'DCs=', DCs
        Ud = params['U'][0]
        if type(self.J) is list:
            Jd = sum(self.J)/len(self.J)
        else:
            Jd = self.J
        # If we want to have U different than the U used in double counting
        # (useful for cluster-DMFT when U can not be simply added in ctqmc)
        if params.has_key('Ud'): Ud = params['Ud'][0]
        if params.has_key('Jd'): Jd  = params['Jd'][0]
        
        DC_ones = array([int(self.Diagonal[i]) for i in range(len(Edc_old))])
        if DCs[:5]=='exact':  
            # We allow DCs = ['exact', 'exact_l', 'exactd', 'exactdl']
            # exact and exactd use impurity occupation. exact uses yukawa&dielectric, while exactd uses dielectric screening only
            # exact_l and exactdl use lattice occupations. exact_l uses yukawa&dielectric, while exactdl uses dielectric screening only
            lmbda = self.UlamJ[1]
            epsilon = self.UlamJ[2]
            nf_ = mom                         # 'exact': Density determined from impurity occupation
            print 'Here we start nl_imp=', nf_, (list(nl_imp)==True)
            if len(DCs)>6 and DCs[6]=='l' and list(nl_imp): 
                nf_ = nl_imp
            print 'Here we have nl_imp=', nf_
            print >> self.fh_info, 'lmbda=', lmbda, 'epsilon=', epsilon
            (Edc0,Phidc0)=dmft_DC.ComputeExactDC(nf_, lmbda, epsilon, self.icase, self.fh_info, projector='projectorw.dat',trans=self.dir+'Trans.dat')
            print >> self.fh_info, 'The Exact Double-counting is: Vdc=', Edc0, 'PhiDC=', Phidc0

        elif DCs=='FLL' or DCs=='default' or DCs=='FLL_l':
            nf_ = nf
            if DCs=='FLL_l' and nf_: 
                nf_ = sum(nl_imp)
            # This is the fully localized limit double-counting by Anisimov
            Edc0 = DC_ones*( Ud*(nf_-0.5) - Jd*(nf_/2.-0.5) )
            Phidc0 = Ud*0.5*nf_*(nf_-1.) - Jd*0.25*nf_*(nf_-2.)
        elif DCs=='AMF':
            nf_ = nf
            if DCs=='AMF_l' and nf_: 
                nf_ = sum(nl_imp)
            Ueff = Ud*(1-1./(2.*(2.*self.l+1.))) - Jd*self.l/(2*self.l+1.)
            Edc0 = DC_ones*Ueff*nf_
            Phidc0 = Ueff*nf_**2/2.
        elif DCs=='fixn' or DCs=='nominal':         # If the scheme if fixn, we keep double-counting potential fixed
            nf0 = params['nf0'][0]
            Vdc = Ud*(nf0-0.5) - Jd*(nf0/2.-0.5)
            Edc0 = DC_ones * Vdc
            Phidc0 = ntot*Vdc
        else:
            print 'ERROR: Unkown double-counting', DCs
            sys.exit(1)
            
        if params.has_key('dDC'):
            Edc0 += array(params['dDC'][0])
            Phidc0 += sum([params['dDC'][0][i]*nf_[i] for i in range(len(nf_))])
            print >> self.fh_info, 'The Double-counting corrected to: Vdc=', Edc0, 'PhiDC=', Phidc0


        Edc = Edc_old # we want to keep Edc<-1000 for orbitals which are projected out!
        for i in range(len(Edc)):
            if (Edc_old[i]>-1000.): # this is true for all orbitals which are kept!
                Edc[i] = Edc0[i]
            else:                   # these orbitals are projected out -> we want to keep sinf=Edc for orbitals projected out!
                sinf[i] = Edc[i]
                

        # Copying data at each iteration to follow the evolution
        shutil.copy2(self.dir+self.sparams['Sig'], self.dir+self.sparams['Sig']+'.'+extn)
        shutil.copy2(self.dir+self.sparams['Gf'], self.dir+self.sparams['Gf']+'.'+extn)
        shutil.copy2(self.dir+'PARAMS.norg', self.dir+'PARAMS.norg'+'.'+extn)
        shutil.copy2(self.dir+'nohup_imp.out', self.dir+'norg.log'+'.'+extn)
        shutil.copy2(self.dir+'hb_fit.txt', self.dir+'hb_fit'+'.'+extn)
        shutil.copy2(self.dir+'ose_hop', self.dir+'ose_hop'+'.'+extn)
        shutil.copy2(self.dir+'nmat.txt', self.dir+'nmat'+'.'+extn)
        shutil.copy2(self.dir+'h0.txt', self.dir+'h0'+'.'+extn)
        shutil.copy2(self.dir+'ose_hop', self.dir+'ose_hop'+'.'+extn)
        
        print >> self.fh_info, 'nimp, TrSigmaG=', ntot, TrSigmaG
        
        # (Phi_DMFT, lnGimp, Fimp, Vdc_nd, zeorb) = GiveImpFunctional(self.dir,self.PARAMS,Edc,self.fh_info)
        
        self.fh_info.flush()
        
        Ry2eV = 13.60569193
        Epotential = TrSigmaG
        fE = open(self.dir+'/Eorb.dat', 'w')
        
        print >> fE, '#  Phidc=', 0.0  # Phidc0
        print >> fE, '#  Tr(Sigma*G)/2=', 0.0  # Epotential
        print >> fE, '#  Tr(Sigma*G)/2-Phidc=', 0.0  # (Epotential-Phidc0)
        print >> fE, '#  Phi_DMFT=', 0.0  # Phi_DMFT
        print >> fE, '#  Phi_DMFT-Phi_DC=', 0.0  # Phi_DMFT-Phidc0
        print >> fE, '#  Tr(log(-Gimp))=', 0.0  # lnGimp
        print >> fE, '#  Fimpurity+TS=', 0.0  # Fimp
        print >> fE, '#  E_pot+E_kin-Tr(w*dD/dw)=', 0.0  # zeorb
        print >> fE, ':IEORB ', 0.0  # (Epotential-Phidc0)/Ry2eV
        print >> fE, ':XEORB ', 0.0  # (Phi_DMFT-Phidc0)/Ry2eV
        print >> fE, ':ZEORB ', 0.0  # (zeorb-Phidc0)/Ry2eV
        fE.close()

        fE = open(self.dir+'/sig.out', 'w')
        print >> fE, '# s_oo=', sinf.tolist()
        print >> fE, '# Edc=', Edc.tolist()

        for iom in range(len(om)):
            print >> fE, ("%20.15f " % om[iom]),
            for b in range(len(Sigo)):
                print >> fE, ("%20.15f %20.15f  " % (Sigo[b,iom].real-Edc[b], Sigo[b,iom].imag)),
            print >> fE
        fE.close()
        
        return ntot

    def FromSingleImpurity_to_Link(self, Eimps, params):
        # First create local cix for single-impurity .
        Eimps_local=[]
        for i in range(len(Eimps)/2):
            Eimps_local.append( 0.5*(Eimps[2*i]+Eimps[2*i+1]) )
        
        print >> self.fh_info, 'Eimps_local=', Eimps_local

        # Creats local Counterparts for Trans.dat
        Nt = len(self.CF)
        Sigind_local=zeros((Nt/2,Nt/2),dtype=int)
        for i in range(Nt/2):
            for j in range(Nt/2):
                Sigind_local[i,j] = self.Sigind[2*i+1,2*j+1]/2

        floc = open(self.dir+'/Trans.dat_local','w')
        # Creates local CF
        Cp = zeros((Nt/2,Nt),dtype=complex)
        Cm = zeros((Nt/2,Nt),dtype=complex)
        p_is_local=True
        m_is_local=True
        for i in range(0,Nt/2):
            Cp[i,:] = (self.CF[2*i,:]+self.CF[2*i+1,:])/sqrt(2)
            Cm[i,:] = (self.CF[2*i,:]-self.CF[2*i+1,:])/sqrt(2)
            if sum(abs(Cp[i,Nt/2:Nt]))>1e-5: p_is_local=False
            if sum(abs(Cm[i,Nt/2:Nt]))>1e-5: m_is_local=False
        if p_is_local:
            CF_local = Cp[:,:Nt/2]
        elif m_is_local:
            CF_local = Cm[:,:Nt/2]
        else:
            print 'ERROR: Something wrong with transformation T2C given in case.indmfl file'
            print >> fh_info, 'ERROR: Something wrong with transformation T2C given in case.indmfl file'

        print >> floc, Nt/2, Nt/2, '# size of Sigind and CF'
        print >> floc, '#---- Sigind follows'
        for i in range(len(Sigind_local)):
            print >> floc, "%3d "*(len(Sigind_local)) % tuple(Sigind_local[i,:])
        print >> floc, '#---- CF follows'
        for i in range(len(CF_local)):
            for j in range(len(CF_local[i])):
                print >> floc, "%12.8f "*2 % (CF_local[i,j].real, CF_local[i,j].imag),
            print >> floc
        floc.close()

        shutil.move(self.dir+'/Trans.dat', self.dir+'/Trans.dat_cluster')
        shutil.copy(self.dir+'/Trans.dat_local', self.dir+'/Trans.dat')

        atom_arg = ''
        for k,p in self.spr.items():
            atom_arg += '\"'+k+'='+str(p)+'\" '
        atom_arg += ' \"Eimp='+str(Eimps_local)+'\" '
        
        print >> self.fh_info, 'atom_arg=', atom_arg
        print atom_arg
        # runs ED for single site DMFT
        cmd = 'cd '+ self.dir + '; ' + self.atom_exe +' ' +atom_arg+' >> nohup.out 2>&1'
        print >> self.fh_info, '.... running ..', cmd
        print '.... running ..', cmd
        subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
        # runs link.py to compute cix file for 2-site cluster
        link_arg=' -U '+str(params['Ud'][0])+' -E "'+str(Eimps)+'"'
        cmd = 'cd '+ self.dir + '; ' + self.env.ROOT + '/link.py' +' ' +link_arg+' >> nohup.out 2>&1'
        print >> self.fh_info, '.... running ..', cmd
        print '.... running ..', cmd
        subprocess.call(cmd,shell=True,stdout=self.fh_info,stderr=self.fh_info)
        
        shutil.move(self.dir+'/Trans.dat_cluster', self.dir+'/Trans.dat')
        # shutil.move(self.dir+'/link_actqmc.cix', self.dir+'/actqmc.cix')
        
    def SigmaOffDiagonal(self,om,Sigo,wbroad,kbroad,fname='Glatt.imp'):
        glatfile = self.dir + '/'+fname
        if not os.path.isfile(glatfile) or os.path.getsize(glatfile)==0: return None
        
        # Reading Gimp, Glatt, and Delta
        Gimp_dat = loadtxt(self.dir + '/'+self.sparams['Gf'] ).transpose()
        Glat_dat = loadtxt(glatfile).transpose()
        Delta_dat = loadtxt(self.dir+'/Delta.imp').transpose()
        # Two frequency meshes, logarithmic and linear
        om0 = Gimp_dat[0]
        Gimp0 = Gimp_dat[1::2] + Gimp_dat[2::2]*1j
        om1 = Delta_dat[0]
        Glat = Glat_dat[1::2] + Glat_dat[2::2]*1j
        Delta = Delta_dat[1::2] + Delta_dat[2::2]*1j
        # Find correspondence between linear and logarithmic mesh
        iom=0
        Gimp=[] # Gimp on logarithmic mesh
        for i in range(len(om1)):
            while (iom<len(om0) and om0[iom]-1e-6<om1[i] ): iom+=1
            if abs(om1[i]-om0[iom-1])>1e-4:
                print 'Problems! Can not find the correspondence between logarithmic and linear meshin Delta.imp and Gf.out!'
            Gimp.append( Gimp0[:,iom-1] )
        Gimp = array(Gimp).T
            
        # Find which columns and rows are correlated?
        Sigind = array(self.Sigind_shifted)
        ind = []
        for i in range(len(Sigind)):
            if Sigind[i,i]>0:
                ind.append(i)
        dim=len(ind)
        # How many times each correlated orbital appears
        noccur = zeros(len(Delta),dtype=int)
        for i in range(len(Sigind)):
            for j in range(len(Sigind)):
                if Sigind[i,j]>0:
                    noccur[Sigind[i,j]-1]+=1
        # Impurity levels in matrix form
        Eimp = matrix(zeros( (dim,dim), dtype=complex))
        for i,ii in enumerate(ind):
            for j,jj in enumerate(ind):
                if (Sigind[ii,jj]>0): Eimp[i,j] = self.Eimps[Sigind[ii,jj]-1]
        
        #print 'self.Eimps=', self.Eimps
        
        Gc=matrix(zeros((dim,dim),dtype=complex))
        Dlt=matrix(zeros((dim,dim),dtype=complex))
        Id = matrix(identity(dim,dtype=complex))
        Sigma = zeros( (len(Delta), len(om1) ), dtype=complex )

        #print 'len(Delta)=', len(Delta), 'shape(Sigma)=', shape(Sigma), 'len(noccur)=', len(noccur), 'len(ind)=', len(ind), 'ind=', ind
        for iom in range(len(om1)):
            for i,ii in enumerate(ind):
                for j,jj in enumerate(ind):
                    if (Sigind[ii,jj]>0):
                        Dlt[i,j] = Delta[Sigind[ii,jj]-1,iom]
                        Gc[i,j]  = Glat[Sigind[ii,jj]-1,iom]
                    Gc[i,i] = Gimp[Sigind[ii,ii]-1,iom]
                    
            Sigm = Id*(om1[iom]*1j)-Dlt-Eimp-Gc.I
            
            for i,ii in enumerate(ind):
                for j,jj in enumerate(ind):
                    if (Sigind[ii,jj]>0):
                        if (ii!=jj and self.sparams['OffDiagonal']=='real'):
                            Sigma[Sigind[ii,jj]-1, iom] += real(Sigm[i,j])
                        else:
                            Sigma[Sigind[ii,jj]-1, iom] += Sigm[i,j]
                            
            for i in range(len(Sigma)): Sigma[i,iom] *= 1./noccur[i]

        savetxt('S1.dat', vstack( (om, real(Sigo), imag(Sigo)) ).transpose() )
        savetxt('S2.dat', vstack( (om1, real(Sigma), imag(Sigma) ) ).transpose() )

        for ib in range(len(Sigma)):
            Sigma[ib] = brod.Broad(wbroad, kbroad, om1, Sigma[ib])
            
        savetxt('S3.dat', vstack( (om1, real(Sigma), imag(Sigma) ) ).transpose() )
        
        nom = self.sparams['nom']
        # Interpolate om big mesh
        Sign=[]
        for ib in range(len(Sigma)):
            fSr = interpolate.UnivariateSpline(om1, Sigma[ib].real, s=0)
            fSi = interpolate.UnivariateSpline(om1, Sigma[ib].imag, s=0)
            Sg = fSr(om) + fSi(om)*1j
            if not self.Diagonal[ib]:
                Sg[nom:] *= exp(-(om[nom:]-om[nom]))
            Sign.append( Sg )
            
        Sign = array(Sign)

        (Eimp,Olap,Edc) = loadtxtE(self.dir+'/Eimp.inp')
        fE = open(self.dir+'/test.out', 'w')
        #print >> fE, '# s_oo=', sinf.tolist()
        #print >> fE, '# Edc=', Edc.tolist()
        for iom in range(len(om)):
            print >> fE, ("%20.15f " % om[iom]),
            for b in range(len(Sign)):
                print >> fE, ("%20.15f %20.15f  " % (Sign[b,iom].real-Edc[b], Sign[b,iom].imag)),
            print >> fE
        fE.close()
        
        return Sign
    
                
def ferm(x):
    """Fermi function for T=1"""
    if (x>100): return 0
    if (x<-100): return 1
    return 1/(exp(x)+1)


if __name__ == '__main__':
    iparams0={"exe"                : ["norg"              , "# Name of the executable"],
          "U"                  : [0.0                 , "# Coulomb repulsion (F0)"],
          "Ud"                 : [6.0                 , "# Coulomb repulsion (F0)"],
          "J"                  : [1.0                  , "# Coulomb repulsion (F0)"],
          "nf0"                : [1                  , "# Double counting parameter"],
          "beta"               : [35.0                 , "# Inverse temperature"],
          "Nmax"               : [2000                  , "# Maximum perturbation order allowed"],
          "M"                  : [1e8                  , "# Total number of Monte Carlo steps"],
          "Ncout"              : [1000000             , "# How often to print out info"],
          "Naver"              : [1000000000          , "# How often to print out debug info"],
          "nom"                : [50                 , "# Number of Matsubara frequency points sampled"],
          "aom"                : [4                  , "# Number of frequency points used to determin the value of sigma at nom"],
          "sderiv"             : [0.02                , "# Maximum derivative mismatch accepted for tail concatenation"],
          "Ntau"               : [1000                , "# Number of imaginary time points (only for debugging)"],
          "SampleGtau"         : [1000                , "# How often to update G(tau)"],
          "GlobalFlip"         : [500000            , "# How often to try a global flip"],
          "tsample"            : [10                  , "# How often to record measurements"],
          "warmup"             : [100000         , "# Warmup number of QMC steps"],
          "CleanUpdate"        : [100000                , "# How often to make clean update"],
          "minM"               : [1e-10               , "# The smallest allowed value for the atomic trace"],
          "minD"               : [1e-10               , "# The smallest allowed value for the determinant"],
          "PChangeOrder"       : [0.9                 , "# Ratio between trial steps: add-remove-a-kink / move-a-kink"],
          "CoulombF"           : ["'Ising'"         , "# Georges = rough run"],
          "OCA_G"              : [False               , "# Don't compute OCA diagrams for speed"],
          }
    
    dire = '.'
    fi = open(dire+'/Trans.dat', 'r')
    dim = int(fi.next().split()[0])
    fi.next()
    Sigind=[]
    for i in range(dim):
        Sigind.append( map(int,fi.next().split()[:dim]) )
    Sigind = array( Sigind )
    fi.next()
    CF=[]
    for i in range(dim):
        dat = map(float,fi.next().split()[:2*dim])
        CF.append( array(dat[::2])+array(dat[1::2])*1j )
    CF = array( CF )
    
    imp = IMP_CTQMC(23, dire, iparams0, Sigind, CF)

    imp._exe(iparams0, 'nominal', '_bris', 1)
    #imp.non_exe(iparams0, '_bris', 0, 0.03, 0.15)
    
    #print imp.HighFrequency(iparams0, 'fixn', 0.0, '_bris', 0.05, 0.1)
    
