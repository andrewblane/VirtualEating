from numpy import *
from pylab import *

import os
import Bio
import re
import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Restriction
from Bio.Restriction import *
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqFeature
from Bio.SeqFeature import *
import itertools
import multiprocessing
import time
from collections import Counter
import operator



#### Plotting output functions

def al_plot_complexity(list):
    '''
    Takes a list of sequences, plots count of sequences vs percent of unique sequences.
    '''
    c = Counter(list)
    freqs = []
    for item in c:
       freqs.append((c[item], c))
    sorted_c = sorted(c.items(), key=operator.itemgetter(1))
    x = xrange(0, len(sorted_c))
    counts = []
    for item, count in sorted_c:
        counts.append(count)
    global rangelist
    rangelist = [a/100.0 for a in xrange(0, 100, 2)]
    y = []
    for p in rangelist:
        a, l = int(len(counts)*p), len(counts)
        y.append(float(sum(counts[l-a:l]))/sum(counts)*100)
    zip(y, rangelist)
    #get_ipython().magic(u'pylab inline')
    figure()
    plot(rangelist, y, "o", color="gray")
    xlabel("Cumulative fraction of unique spacers")
    ylabel("Percent of total reads (molecules) in library")

#### Textual output functions

def al_print_features(inputseq, addpamcutters, directlabel):
    '''
    Takes 3 arguments:
    inputseq: SeqRecord to draw
    addpamcutters: int - 0 or 1. If 1, also draw HpaII/BfaI/ScrFI sites on map.
    directlabel: int, 0 or 1. Recommend 1. Changes position of feature labels (1 = on markings, 0 = below them)
    '''
    if addpamcutters == 1:
        cutline = list(" " * len(inputseq))
        HpaIIsites = HpaII.search(inputseq.seq)
        BfaIsites = BfaI.search(inputseq.seq)
        ScrFIsites = ScrFI.search(inputseq.seq)
        for cut in HpaIIsites:
            cutline[cut-1:cut + len("HpaII")] = "<HpaII"
        for cut in BfaIsites:
            cutline[cut-1:cut + len("BfaI")] = "<BfaI"
        for cut in ScrFIsites:
            cutline[cut-1:cut + len("ScrFI")] = "<ScrFI"
        cutline = "".join(cutline)


    mask = [list((("-" * 9) + "^" )* int(round(len(inputseq.seq)/10.0)))]
    newmaskline = list((("-" * 9) + "^" )* int(round(len(inputseq.seq)/10.0)))

    for feature in inputseq.features:
        # Make a new marker strand if any features overlap. All marker strands can be elements of a list.
        featstart = int(feature.location.start)
        featend = int(feature.location.end)
        if featstart > featend:
            print("Error! Feature end must be after feature start. Use strand to specify direction! Feature " + feature.type +                   " will not be displayed!")
        #if "<" in mask[-1][featstart:featend] or ">" in mask[-1][featstart:featend]:
            #mask.append(newmaskline)

        clean = 0
        for item in mask[-1][featstart:featend]:
            if item == "-":
                clean = 1
            elif item == "^":
                clean = 1
            else:
                clean = 0
                mask.append(newmaskline)
                break

        #print mask[-1][0:50]
        if feature.strand == 1:
            mask[-1] = mask[-1][:featstart-1] + [">"] * int(featend - featstart + 1) + mask[-1][featend:]
            if directlabel == 1:
                mask[-1] = mask[-1][:featstart] + list(str(feature.type)) + mask[-1][featstart + len(str(feature.type)):]
        if feature.strand == -1:
            mask[-1] = mask[-1][:featstart-1] + ["<"] * int(featend+1 - featstart) + mask[-1][featend:]
            if directlabel == 1:
                mask[-1] = mask[-1][:featstart+1] + list(str(feature.type)) + mask[-1][featstart+2 + len(str(feature.type)):]

    #if addpamcutters = 1:
        #cutline = list(" " * len(inputseq)
        #HpaIIsites = HpaII.search(inputseq.seq)



    for index, maskline in enumerate(mask):
        maskline = "".join(maskline)
        mask[index] = maskline
    # add labels
    if directlabel == 0:
        masklab = list(" " * (len(inputseq.seq)))
        for feature in inputseq.features:
            featstart = int(feature.location.start)
            featend = int(feature.location.end)
            featname = str(feature.type)
            masklab = masklab[:featstart] + list(str(feature.type)) + list(" " * (featend-1 - featstart - len(feature.type))) + masklab[featend-1:]
        masklab = "".join(masklab)

    lines = int(round(len(inputseq.seq) / 100)) + 1
    i = 0
    outstring = list(inputseq.name + "\n")
    #print inputseq.name
    outstring = []
    while i < lines:
        indexstart = i*100
        indexend = (i+1) * 100
        if indexend > len(inputseq.seq):
            indexend = len(inputseq.seq)
        if addpamcutters ==1:
            outstring.extend((str(indexstart + 1)) + "  " + cutline[indexstart:indexend] + "   " + str(indexend)+ "\n")
            #print (str(indexstart + 1)) + "  " + cutline[indexstart:indexend] + "   " + str(indexend)
        outstring.extend(str(indexstart+1) + "  " + inputseq.seq[indexstart:indexend] + "   " + str(indexend)+ "\n")
        #print str(indexstart+1) + "  " + inputseq.seq[indexstart:indexend] + "   " + str(indexend)
        for maskline in mask:
            outstring.extend((str(indexstart + 1)) + "  " + maskline[indexstart:indexend] + "   " + str(indexend)+ "\n")
            #print (str(indexstart + 1)) + "  " + maskline[indexstart:indexend] + "   " + str(indexend)
        if directlabel == 0:
            outstring.extend(str(indexstart +1) + "  " + masklab[indexstart:indexend] + "   " + str(indexend)+ "\n")
            #print str(indexstart +1) + "  " + masklab[indexstart:indexend] + "   " + str(indexend)
        outstring.extend("\n")
        #print "\n"
        i = i + 1
    outstring = "".join(outstring)
    return outstring



def al_string2feat(queryseq, ampsdict): #lib5pr is subjectseq; t7 is queryseq
    '''
    This function accepts a query seq and a dictionary of subjectseqs, where the key (amp)
    is contained in a field in queryseq, highlighting the location of queryseq in it.
    Returns a Seq object with new features as specified in ampsdict
    '''
    subjectseq = SeqRecord(ampsdict[queryseq[1][0]])
    #for seqrecord in subjectseq:
    locstart = queryseq[1][1]
    #print queryseq
    locend = queryseq[1][2]
    fwdlocs = []
    revlocs = []
    # Figure out which strand the BLAST hit is on
    if locstart <= locend:
        fwdlocs.append(locstart)
    if locstart > locend:
        revlocs.append(locend)

    for item in fwdlocs:
        start = ExactPosition(int(item))
        end = ExactPosition(int((item) + len(queryseq[0].seq) + 1))
        location = FeatureLocation(start, end)
        feature = SeqFeature(location,type=str("cutsite_fwd"), strand = +1)
        subjectseq.features.append(feature)

    for item in revlocs:
        start = ExactPosition(int(item))
        end = ExactPosition(start + len(queryseq[0].seq))
        location = FeatureLocation(start, end)
        feature = SeqFeature(location,type=str("cutsite_rev"), strand = -1)
        subjectseq.features.append(feature)
    #print subjectseq.features
    return subjectseq
