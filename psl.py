#!/usr/bin/env python

""" PSL - PSL format from the BLAT program """

# Created by Vince Forgetta, Jan. 2010.

import math

class Psl:
    ''' Class to represent PSL output from BLAT.
    Provides basic methods such as score() and calcPercentIdentity().
    Argument: String in PSL format
    '''
    
    def __init__(self, s):
        
        # split and tokenize input
        fields = s.strip().split()
        num_fields = len(fields)
        matches, mismatches, repmatches, ncount, qnuminsert, qbaseinsert, \
            tnuminsert, tbaseinsert, strand, qname, qsize, qstart, qend, \
            tname, tsize, tstart, tend, blockcount, blocksizes, qstarts, \
            tstarts = fields[0:21]
        
        # if pslx format
        self.qblockseqs = None
        self.tblockseqs = None
        if num_fields == 23:
            self.qblockseqs = fields[21].split(',')[0:-1]
            self.tblockseqs = fields[22].split(',')[0:-1]
        
        self.matches = int(matches)
        self.mismatches = int(mismatches)
        self.repmatches = int(repmatches)
        self.ncount = int(ncount)
        self.qnuminsert = int(qnuminsert)
        self.qbaseinsert = int(qbaseinsert)
        self.tnuminsert = int(tnuminsert)
        self.tbaseinsert = int(tbaseinsert)
        self.strand = strand
        self.qname = qname
        self.qsize = int(qsize)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.tname = tname
        self.tsize = int(tsize)
        self.tstart = int(tstart)
        self.tend = int(tend)
        self.blockcount = int(blockcount)
        self.blocksizes = [int(x) for x in blocksizes.split(',')[0:-1]]
        self.qstarts = [int(x) for x in qstarts.split(',')[0:-1]]
        self.tstarts = [int(x) for x in tstarts.strip().split(',')[0:-1]]
        
    ## Private methods
        
    def __lenmul(self):
        ''' Determine length multiplier.
        Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        if self.__isProtein:
            return 3
        else:
            return 1

    def __isProtein(self):
        ''' Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        lastblock = self.blockcount - 1
        return (self.strand[1:1] == '+' and \
                    self.tend == (self.tstarts[lastblock] + (3 * self.blocksizes[lastblock]))) or \
                    ((self.strand[1:1] == '-') and \
                         (self.tstart == (self.tsize - (self.tstarts[lastblock] + 3*self.blocksizes[lastblock]))))
    
    def __calcMilliBad(self, ismrna):
        ''' Return number of non-identical matches. 
        Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        qalisize = self.__lenmul() * self.qspan()
        alisize = min(qalisize, self.tspan())
        millibad = 0
        if alisize <= 0: return 0
        sizediff = alisize - self.tspan()
        if sizediff < 0:
            if ismrna:
                sizediff = 0
            else:
                sizediff = -sizediff
        insertfactor = self.qnuminsert
        if not ismrna: insertfactor += self.tnuminsert
        total = self.__lenmul() *\
            (self.matches + self.repmatches + self.mismatches)
        if total != 0:
            millibad = (1000 * (self.mismatches * self.__lenmul() + insertfactor + \
                                    round(3*math.log(1 + sizediff)))) / total
        return millibad

    # Public methods
    
    def qspan(self):
        ''' Span of alignment on query sequence '''
        return self.qend - self.qstart
    
    def tspan(self):
        ''' Span of alignment on target sequence '''
        return self.tend - self.tstart
    
    def score(self):
        ''' Score as calculated by web-BLAT. 
        Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        return self.matches + (self.repmatches / 2) - self.mismatches - \
            self.qnuminsert - self.tnuminsert
    
    def calcPercentIdentity(self):
        ''' Percent identity as calculated by web-BLAT. 
        Adapted from http://genome.ucsc.edu/FAQ/FAQblat#blat4 '''
        return 100.0 - self.__calcMilliBad(True) * 0.1
    

if __name__ == '__main__':
    
    # Replicate BLAT output from web
    for line in open('line.psl'):
        p = Psl(line)
        print p.qname, p.score(), p.qstart+1, p.qend, p.qsize, \
            "%.1f" % p.calcPercentIdentity(), p.tname, p.strand, \
            p.tstart+1, p.tend, p.tspan()        
    
