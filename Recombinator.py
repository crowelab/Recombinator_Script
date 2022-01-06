#!/usr/bin/env python 
import numpy as np
import scipy.stats as stats
from subprocess import call
import csv
import re
import sys
import os
import errno
import argparse
import glob
import pickle
#subprocess.call(['python', 'CsvReader3.py', 'list_fastafiles', 'IGHV3-7_IGHJ1_5.fasta','IGHV3-7_IGHJ1_5.fasta','IGHV3-7_IGHJ1_5.fasta'])
try:
    from Bio import SeqIO
    from Bio.Seq  import Seq
    from Bio.Alphabet import IUPAC
    from Bio.SeqRecord import SeqRecord
except ImportError:
        print 'Biopython is not installed or is not in the path of python'
        sys.exit(1)
# for read in Stuff[k]:
#                Reads.append(SeqRecord(Seq(read[1]),id=read[0],description=''))

'''-----------------------------------------------------'''
''' Set the random number generator to validate entries '''
'''-----------------------------------------------------'''
np.random.seed(seed=1010)

class RandomGermline():
    def __init__(self,Commands):
        self.commands=Commands
        self.seq={}
        self.path_to_germ=''
        self.low=0
        self.high=0
        try:
            self.size=self.commands['size']
            self.path_to_germ=self.commands['path_germ']
        except KeyError:
            pass
        '''---------------------'''
        '''  READ IN GERMLINES  '''
        '''---------------------'''
        self._readInGermline(self.path_to_germ,self.seq)
        
    def _numberEntries(self):
        return np.random.randint(self.low,self.high, size=self.size)
    
    def _readInGermline(self,path,dictionary):
        for name in glob.glob(path+'/*.fasta'):
            dictionary[os.path.splitext(os.path.basename(name))[0]]={rec.id : rec.seq for rec in SeqIO.parse(name, "fasta")}
            
    def random(self,germline,size):
        ''' Use regex to search for IGH, IGJ '''
        try:
            reads=[]
            entries=self.seq[germline].keys()
            self.low=0
            self.high=len(entries)
            self.size=size
            random_entries=self._numberEntries()
            #print random_entries
            #print entries
            #print self.seq
            for i in random_entries:
                if self.seq[germline].has_key(entries[i]):
                    reads.append(self.seq[germline][entries[i]])
            return reads
        except KeyError:
            pass
    
            
class Recombine():
    '''Class to read in CSV file containing a header. '''
    def __init__(self, Commands):
        self.commands=Commands
        self.VSegments=[]
        self.DSegments=[]
        self.JSegments=[]
        self.ProductiveReads=[]
        self.ProductiveDnaCdr3=[]
        self.ProductiveProteinCdr3=[]
        self.ProductiveGermlineDivergence = []
        self.Cdr3=[]
        self.GermlineDivergence=[]
        self.Contigs=[]
        self.Nuc=[None] * 4
        self.Nuc[0]='A'
        self.Nuc[1]='T'    
        self.Nuc[2]='C'
        self.Nuc[3]='G'
        try:
            self.vgermline_gene=self.commands['vgerm']
            self.dgermline_gene=self.commands['dgerm']
            self.jgermline_gene=self.commands['jgerm']
            self.recombination()
            self._generateContigs()
        except KeyError:
            pass

    def recombination(self):
        '''-------------------------------'''
        ''' Create the distributions here '''
        '''-------------------------------'''
        v3Distribution=self._generateTruncatedGaussian(self.commands['v3mu'],self.commands['v3sigma'],self.commands['size']).round().astype(int)
        d5Distribution=self._generateTruncatedGaussian(self.commands['d5mu'],self.commands['d5sigma'],self.commands['size']).round().astype(int)
        d3Distribution=self._generateTruncatedGaussian(self.commands['d3mu'],self.commands['d3sigma'],self.commands['size']).round().astype(int)
        j5Distribution=self._generateTruncatedGaussian(self.commands['j5mu'],self.commands['j5sigma'],self.commands['size']).round().astype(int)
        vdDistribution=self._generateTruncatedGaussian(self.commands['vdmu'],self.commands['vdsigma'],self.commands['size']).round().astype(int)
        djDistribution=self._generateTruncatedGaussian(self.commands['djmu'],self.commands['djsigma'],self.commands['size']).round().astype(int)

        for i in np.arange(self.commands['size']):
            Vsequence=list(self.commands['vgerm'][1][i])
            Dsequence=list(self.commands['dgerm'][1][i])
            Jsequence=list(self.commands['jgerm'][1][i])

            '''---------------------------------------'''
            ''' Exonuclease portion of the simulation '''
            '''---------------------------------------'''
            if v3Distribution[i] > 0:
                Vsequence=Vsequence[:-v3Distribution[i]]
            Jsequence=Jsequence[j5Distribution[i]:]
            Dsequence=Dsequence[d5Distribution[i]:]
            if d3Distribution[i] > 0:
                Dsequence=Dsequence[:-d3Distribution[i]]
        

            '''-------------------------------'''
            ''' TdT portion of the simulation '''
            '''-------------------------------'''
            for j,val in enumerate(np.random.randint(4, size=vdDistribution[i])):
                Dsequence.insert(j, self.Nuc[val])
            for k,val in  enumerate(np.random.randint(4, size=djDistribution[i])): 
                Dsequence.append(self.Nuc[val])

            '''--------------------- '''
            '''    Store arrays      '''
            '''----------------------'''
            self.VSegments.append(Vsequence)
            self.JSegments.append(Jsequence)
            self.DSegments.append(Dsequence)

    def shm(self):
        self.Contigs=[]
        self.Cdr3=[]
        ''' ------------------------------------------------------'''
        ''' Allow uniform mutation througout contiguous sequences '''
        '''-------------------------------------------------------'''
        if self.commands['fixed_mut_p'] == False:
            shm_prob=np.random.uniform(0,self.commands['mut_p'],self.commands['size'])

        for i in np.arange(self.commands['size']):
            ''' ------------------------------------ '''
            ''' Create contiguous read to work wtih  '''
            '''------------------------------------- '''
            dna=self.VSegments[i]+self.DSegments[i]+self.JSegments[i]

            '''-----------------------------------------------'''
            ''' Determine the number of nucleotides to mutate '''
            '''-----------------------------------------------'''
            if self.commands['fixed_mut_p']:
                number_nt_mutate=int(len(dna)*self.commands['mut_p'])
            else:
                number_nt_mutate=int(len(dna)*shm_prob[i])

            ''' --------------------------------------------------------'''
            ''' Determine the position in the read that will be mutated '''
            '''---------------------------------------------------------'''
            positions=np.random.randint(len(dna), size=number_nt_mutate)
            #print number_nt_mutate
            #print shm_prob[i]*100.00
            #print positions
            #print ''.join(dna)
            ''' ----------------------------------- '''
            ''' Mutate the positions in the congig  '''
            ''' ----------------------------------- '''
            for j, val in enumerate(np.random.randint(4, size=number_nt_mutate)):
                dna[positions[j]]=self.Nuc[val]

            ''' ------------------------------------------------------------------ '''
            ''' Check to see if cannonical TGT or TGA appears in VSegment+DSegment '''
            ''' ------------------------------------------------------------------ '''
            length=len(self.VSegments[i])
            triplet_length=length/3
            vdj=list(re.findall('..{0,1}.{0,1}',''.join(dna)))
            # vsegment_reversed = list(reversed(re.findall('..{0,1}.{0,1}', ''.join(self.VSegments[i] + self.DSegments[i][:5]))))
            cysteine=-1
            wgxg=-1
            for i in range(triplet_length-4,len(vdj)-4):
                if vdj[i] == 'TGT':
                    cysteine=i
                if vdj[i]=='TGG' and (vdj[i+1][0]=='G' and vdj[i+1][1]=='G') and (vdj[i+3][0]=='G' and vdj[i+3][1]=='G'):
                    wgxg=i

            if (cysteine >0 and wgxg >0):
                #print vdj[cysteine:wgxg+4]
                self.Contigs.append(dna)
                self.Cdr3.append(''.join(vdj[cysteine:wgxg+4]))
                self.GermlineDivergence.append(number_nt_mutate)
    
    def _generateContigs(self):
        if self.commands['shm']:
            self.shm()
        else:
            for i in np.arange(self.commands['size']):

                dna=self.VSegments[i]+self.DSegments[i]+self.JSegments[i]

                ''' ------------------------------------------------------------------ '''
                ''' Check to see if cannonical TGT or TGA appears in VSegment+DSegment '''
                ''' ------------------------------------------------------------------ '''
                length = len(self.VSegments[i]) / 3
                vdj = list(re.findall('..{0,1}.{0,1}', ''.join(dna)))
                cysteine = -1
                wgxg = -1
                for i in range(length - 4, len(vdj) - 4):
                    if vdj[i] == 'TGT':
                        cysteine = i
                    if vdj[i] == 'TGG' and (vdj[i + 1][0] == 'G' and vdj[i + 1][1] == 'G') and (
                            vdj[i + 3][0] == 'G' and vdj[i + 3][1] == 'G'):
                        wgxg = i

                if (cysteine > 0 and wgxg > 0):
                    self.Contigs.append(dna)
                    self.Cdr3.append(''.join(vdj[cysteine:wgxg + 4]))

        if self.commands['productive']:
            self._productive()

    def _stopCodons(self):
        try:
            for i in np.arange(self.commands['size']):
                dna=''.join(self.Contigs[i])
                protein=Seq(dna).translate()
                if '*' not in protein:
                    self.ProductiveReads.append(SeqRecord(Seq(dna),id='%10.10d'%i,description=''))
        except IndexError:
            pass

    def _productive(self):
        try:
            for i,idx in enumerate(self.Contigs):
                dna=''.join(self.Contigs[i])
                protein=Seq(dna).translate()
                if '*' not in protein:
                    #if self._patternScan('WG.G|WG.|WG',str(protein)[-25:]):
                    self.ProductiveReads.append(SeqRecord(Seq(dna),id='%10.10d'%i,description=''))
                    self.ProductiveDnaCdr3.append(self.Cdr3[i])
                    self.ProductiveProteinCdr3.append(str(Seq(self.Cdr3[i]).translate()))
                    self.ProductiveGermlineDivergence.append(self.GermlineDivergence[i])
        except IndexError:
            pass

    def _patternScan(self,AnchorString,Sequence):
        Ranges = []
        Vocabulary = ''
        for record in AnchorString:
            Vocabulary = Vocabulary + record
        for m in re.finditer(r'%s' % Vocabulary, Sequence):
            Ranges.append((m.start(), m.end()))
        return (Ranges)

    def _generateTruncatedGaussian(self,mu,sigma,size):
        lower, upper = 0, mu+3*sigma
        return stats.truncnorm(
            (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma).rvs(size)

    def getContigs(self,selection=''):
        if self.commands.has_key('productive') and self.commands['productive'] or selection=='productive':
            return self.ProductiveReads
        else:
            return self.Contigs

    def getCdr3(self, selection=''):
        if self.commands.has_key('productive cdr3') and self.commands['productive cdr3'] or selection=='productive cdr3':
            return self.ProductiveProteinCdr3
        else:
            return self.Cdr3

def main():
    vseq={}
    sample_size=10000000
    TotalIncrement=10
    base = '/home/sotocs/Scratch/Recombinator/Recombination'
    '''------------------------------------'''
    ''' Generate the V-germline genes here '''
    '''------------------------------------'''
    vgerm = RandomGermline({'size': sample_size,
                            'path_germ': base + '/' + 'IGHV1-69'})
    dgerm = RandomGermline({'size': sample_size,
                            'path_germ': base + '/' + 'IGHD3-3'})
    jgerm = RandomGermline({'size': sample_size,
                            'path_germ': base + '/' + 'IGHJ6'})
    for increment in np.arange(sample_size,TotalIncrement*sample_size+sample_size,sample_size):

        vseq = vgerm.random(germline='IGHV1-69',size=sample_size)

        '''------------------------------'''
        '''Grab the D-germline genes here'''
        '''------------------------------'''
        dseq = dgerm.random(germline='IGHD3-3',size=sample_size)

        '''----------------------------------'''
        '''Generate the J-germline genes here'''
        '''----------------------------------'''
        jseq = jgerm.random(germline='IGHJ6',size=sample_size)

        ''' ----------------------------- '''
        '''  HOW MANY READS DO YOU HAVE   '''
        ''' ----------------------------- '''
        P=Recombine({'vgerm':('IGHV1-69',vseq),
                     'dgerm':('IGHD3-3',dseq),
                     'jgerm':('IGHJ6',jseq),
                     'v3mu':2,
                     'v3sigma':1,
                     'd5mu':5,
                     'd5sigma':1,
                     'd3mu':2,
                     'd3sigma':1,
                     'j5mu':7,
                     'j5sigma':2,
                     'vdmu':2,
                     'vdsigma':1,
                     'djmu':2,
                     'djsigma':1,
           'size':sample_size,
           'mut_p':0.30,
           'fixed_mut_p': False,
           'shm': True,
           'productive':True
                   })
        #print len(P.getCdr3(selection='productive cdr3'))
        #print P.getCdr3(selection='productive cdr3')
        filehandle=open('%d.pickle'%increment,'wb')
        pickle.dump(P.getCdr3(selection='productive cdr3'),filehandle)
        filehandle.close()
        #filehandle=open('%d'%sample_size,'rb')
        #print pickle.load(filehandle)
        del P
        del vseq
        del dseq
        del jseq
    ''' ------------------------------ '''
    ''' Compute some basic statistics  '''
    ''' ------------------------------ '''
    Cdr3dict={}
    for increment in np.arange(sample_size, 10 * sample_size + sample_size, sample_size):
        filehandle = open('%d.pickle' % increment,'rb')
        #Cdr3=Cdr3+pickle.load(filehandle)
        for cdr3 in pickle.load(filehandle):
            if Cdr3dict.has_key(cdr3):
                Cdr3dict[cdr3]=Cdr3dict[cdr3]+1
            else:
                Cdr3dict[cdr3]=1
        filehandle.close()
    print len(Cdr3dict)
    #print len(Cdr3), len(list(set(Cdr3)))
    #print map(len, Cdr3)
    #total_avg = sum( map(len, strings) ) / len(strings)
    #SeqIO.write(P.getContigs(),'Test.fasta','fasta')
if __name__ == "__main__":
     main()

