from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqFeature
from Bio.SeqFeature import *
import re

# Define the filter that delineates the end of a target 20mer. This uses the first 9nt of the sgRNA hairpin.
# Makes a Feature for each of the marker 9nt sgRNA signatures in the FASTQ sequencing files.

def spacersonly(seqs):
    sgRNAconst = SeqRecord(Seq("GTTTAAGAG"))
    while True:
        seqrecord = seqs.next()
        #for seqrecord in seqs:
        fwdlocs = []
        revlocs = []
        fwdlocs = [tloc.start() for tloc in re.finditer(str(sgRNAconst.seq), str(seqrecord.seq))]
        for item in fwdlocs:
            start = ExactPosition(int(item) + 1)
            end = ExactPosition(int((item) + len(sgRNAconst)))
            location = FeatureLocation(start, end)
            feature = SeqFeature(location,type="sgRNAconst", strand = +1)
            seqrecord.features.append(feature)
        revlocs = [tloc.start() for tloc in re.finditer(str(sgRNAconst.reverse_complement().seq), str(seqrecord.seq))]
        for item in revlocs:
            start = ExactPosition(int(item) + 1)
            end = ExactPosition(start + len(sgRNAconst) - 1)
            location = FeatureLocation(start, end)
            feature = SeqFeature(location,type="sgRNAconst", strand = -1)
            seqrecord.features.append(feature)
        for feat in seqrecord.features:
            if feat.strand == 1:
                tgtstart = int(feat.location.start) - 36 # -21
                tgtend = int(feat.location.start) - 1
                sgtgt = seqrecord[tgtstart:tgtend]
                #yield sgtgt
                #alltgts.append(sgtgt)
                #print "pos \n \n"
            if feat.strand == -1:
                tgtend = int(feat.location.end) + 36 # +21
                tgtstart = int(feat.location.end)
                sgtgt = seqrecord[tgtstart:tgtend].reverse_complement()
                sgtgt.name=seqrecord.name
                #yield sgtgt
                #alltgts.append(sgtgt)
            bad = 0
            try:
                l = [tloc.end() for tloc in re.finditer("ACTCACTATAG", str(sgtgt.seq))]
                sgtgt = sgtgt[int(l[0]):]
            except:
                None
            for score in sgtgt.letter_annotations["phred_quality"]:
                if score < 30:
                    bad = 1
            if bad == 0 and len(sgtgt) > 10:
                yield sgtgt
                break
            #except:
                #yield None
            #sys.stdout.flush() # These two lines produce the line counter as the loop runs
            #sys.stdout.write('\r' + str(index) + " "),
