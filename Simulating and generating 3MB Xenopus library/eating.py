from __future__  import print_function
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
import random
import itertools
import multiprocessing
import time
from collections import Counter
import operator
#import MySQLdb
from numpy import *
from pylab import *


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

def al_string2feat(queryseq, ampsdict): #lib5pr is subjectseq; t7 is queryseq
    '''
    This function accepts a query seq and a dictionary of subjectseqs, where the key (amp)
    is contained in a field in queryseq, highlighting the location of queryseq in it.
    Returns a string.
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

def al_digesttarget(list_of_targets, process_to_file=0):
    '''
    Substrate should be a list of SeqRecords (i.e. scaffolds, or whatever).
    Process_to_file [optional] is an integer; if 1 will avoid building list in memory and will write to file instead.
        Defaults to 0, building list in memory.
    Output list has information in it:
        ID is a count of the potential spacer along the scaffold, starting from zero.
        Name is the base position on the scaffold.
        Description is the enzyme that generated the spacer-end (ScrFI/HpaII/BfaI)
        dbxrefs is the name the scaffold was given in the input SeqRecord
    '''
    f = open("genome_allspacers.fa","w")
    global substrate
    scaffoldcuts = []
    for substrate in list_of_targets:
        cutslist = []
        pos = HpaII.search(substrate.seq)
        # Positions in this list correspond to the right boundaries of fragments;
        # last one is thus the sequence end
        pos.append(len(substrate.seq))
        pos = iter(pos)
        cuts = HpaII.catalyze(substrate.seq)
        for item in cuts:
            cutslist.append([item, "HpaII", int(pos.next())])

        cuts = BfaI.catalyze(substrate.seq)
        pos = BfaI.search(substrate.seq)
        pos.append(len(substrate.seq))
        pos = iter(pos)
        for item in cuts:
            cutslist.append([item, "BfaI", int(pos.next())])

        cuts = ScrFI.catalyze(substrate.seq)
        pos = ScrFI.search(substrate.seq)
        pos.append(len(substrate.seq))
        pos = iter(pos)
        for item in cuts:
            cutslist.append([item, "ScrFI", int(pos.next())])

        #The above is all to get the results of a catalyze operation (i.e. tuples) into
        # a list format. Next part makes them into SeqRecords.

        cutslistrecords = []
        for i, item in enumerate(cutslist):
            cutslistrecords.append(SeqRecord(item[0], id = str(i), description = str(item[1]), name=str(item[2]),                                              dbxrefs=[str(substrate.id)]))
        cutslist = cutslistrecords[:]

        # This part takes the 3' 20nt of each fragment and makes a new sequence with it.
        # For the 5' end, the Mung-Bean treatment is simulated by removing two more nt (for HpaII and BfaI), or one nt for ScrFI;
        # these would be the 5' overhang. Then we take the reverse-complement of the sequence.
        # The Restriction module just returns sequences as if the top strand only was being cut. In other words,
        # no bases are deleted from consecutive fragments.

        from Bio.Seq import MutableSeq
        twentymers = []
        #record2 = []
        for record2 in cutslist:
                try: # This is because a second run of this code on already mutable seqs seems to fail. Not sure how to flush out and revert back to non-mutables...
                    record2.seq = record2.seq.tomutable()
                except:
                    pass
                if record2.description == "ScrFI":
                    #offset here (e.g. 1:21 is for simulating MBN digeston)
                    # because the entry.names are rooted on the right of each fragment, the length
                    # of the entry.name has to be subtracted to get the desired left position for the "reverse"
                    # tgts
                    entry = record2[1:21].reverse_complement                    (description=True, id=True, name=True)
                    entry.name = str(int(record2.name)+1 - len(record2.seq))
                    entry.id = str("%s_R" % record2.id) ##
                    twentymers.append(entry)
                else: # Should work for HpaII/BfaI
                    entry = record2[2:22].reverse_complement                    (description=True, id=True, name=True)
                    entry.name = str(int(record2.name)+2 - len(record2.seq))
                    entry.id = str("%s_R" % record2.id) ##
                    twentymers.append(entry)
                record2.seq = record2.seq.toseq()
                record2.id = str("%s_F" % record2.id)
                entry = record2[-20:]
                entry.name = str(int(record2.name)-20)
                twentymers.append(entry)

        for item in twentymers:
            item.dbxrefs = [substrate.id]

        # The ends of the fragments aren't bonafide CRISPR targets; these can be removed:
        noends = []
        twentymerstr = [item for item in twentymers if item.description == "HpaII"]
        trimmed = twentymerstr[1:-1] # removes first and last 20mer
        noends.append(trimmed)
        twentymerstr = [item for item in twentymers if item.description == "BfaI"]
        trimmed = twentymerstr[1:-1]
        noends.append(trimmed)
        twentymerstr = [item for item in twentymers if item.description == "ScrFI"]
        trimmed = twentymerstr[1:-1]
        noends.append(trimmed)

        # index = str(substrate.id + str(random.random()))
        if process_to_file != 1:
            scaffoldcuts.append((item for sublist in noends for item in sublist))
        if process_to_file == 1:
            f = open("genome_allspacers.fa","a") #opens file with name of "test.txt"
            for item in [item for sublist in noends for item in sublist]:
                f.write(">lcl|" + str(substrate.id) + "|" + str(item.description) + "|" + str(item.name) + "|" + str(item.id) + "\n")
                f.write(str(item.seq) + "\n")
            f.close()
    if process_to_file != 1:
        return itertools.chain.from_iterable(scaffoldcuts) #dammit, the from_iterable part is important here. that took a while.

def al_scoreguides(guides, genome, genomedict, hits=50, return_blast =1):
    '''

    To produce a score for a given guide in a given BLAST database. Returns an int, the score.
    Only evaluates a max of 500 hit sequences (i.e. really bad guides will not necessarily have a super accurate score.)
    This version only searches the plus strand of the BLAST db, meant to be used with PAM BLAST DBs generated such that
    all potential guide hits are on the plus strand.
    Arguments:  guide = a SeqRecord object containing a proposed guide sequence (or a list of SeqRecords)
                genome = a BLAST database set up on this machine.
                genomedict = a dict made from a FASTA file identical to the one used to make the BLAST DB. Dict keys should be
                    BLAST db hit_def.
    Returns: Tuple containing (finalscore, a list of the individual hits, their locations and their subscores)
    '''
    if isinstance(guides, list) == False:
        guides = [guides]

    name = multiprocessing.current_process().name
    M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0 ,0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
    #import random
    #random.jumpahead(1)
    #filesuffix = str(int(random.random()*100000000))
    filesuffix = str(guides[0].id)
    filename = str("currseq" + filesuffix + ".tmp")
    Bio.SeqIO.write(guides, filename, "fasta")
    # Added dust="no" to stop inflating scores from polyN repeats...
    blastn_cline = NcbiblastnCommandline(query=filename, db=genome,     task = "blastn-short",outfmt=5, out=filename + ".blast", max_target_seqs=100, num_threads = 7, evalue = 10, dust="no")
    #timeit.timeit(blastn_cline, number =1)
    blastn_cline()
    result_handle = open(filename + ".blast")
    blasts = NCBIXML.parse(result_handle)
    # This generates an generator of BLAST objects, each corresponding to a guide query. Loop through those.
    for guideindex, item in enumerate(blasts):
        firstmatch = 0
        notfound = 1
        scorelist = []
        pamlist = []
        # Within each queried guide, loop through the aligned scaffolds
        for scaffold in item.alignments:
            # Within each aligned scaffold, loop through invidivual HSPs on that scaffold.
            for pair in scaffold.hsps:
                hit_threeprime_offset = len(guides[0]) - pair.query_end
                start = pair.sbjct_start
                end = pair.sbjct_end
                if end > start:
                    pamstart = end + hit_threeprime_offset
                    pam = genomedict[scaffold.hit_def][pamstart:pamstart+3]
                elif start > end:
                    pamstart = end - hit_threeprime_offset
                    pam = genomedict[scaffold.hit_def][pamstart-4:pamstart-1].reverse_complement()

                # Construct a bar-based match string, padded to query length,
                # where bar = match position and space is non-match
                mmloc = []
                mmstr = list(pair.match)
                if pair.query_start > 1:
                    mmstr = list(" " * (pair.query_start - 1)) + mmstr
                if pair.query_end < 20:
                    mmstr = mmstr + list(" " * (20 - pair.query_end))
                mmstr = "".join(mmstr)

                # Test for PAM adjacency
                if len(pam) == 3:
                    if (pam[1] == "G" or pam[1] == "A" or pam[1] == "g" or pam[1] == "a") and (pam[2] == "G" or pam[2] == "g"):
                        if (pair.positives >16 and pair.positives < 20):
                            pos = 20
                            for linestring in mmstr:
                                if linestring != "|":
                                    mmloc.append(pos)
                                pos = pos - 1

                            # Actually implement Zhang lab algorithm
                            mmscore = [21 -x for x in mmloc]
                            t1 = 1
                            for mismatchlocation in mmscore:
                                t1 = t1 * (1.0 - float(M[mismatchlocation - 1]))
                            if len(mmscore) > 1:
                                d = (float(max(mmscore)) - float(min(mmscore))) / float((len(mmscore) - 1))
                            else:
                                d = 19
                            t2 = 1 / ( (((19.0 - d)/19.0) * 4) + 1)
                            t3 = float(1)/ float(pow(len(mmloc), 2))
                            scorelist.append({"match_score": float(t1 * t2 * t3 * 100), "scaffold": scaffold.hit_def,                                               "hit_location":pair.sbjct_start, "hit_sequence": pair.sbjct, "pam":str(pam.seq),                                              "match_bars":mmstr})
                            #pamlist.append(pam)
                        # Zhang lab algorithm doesn't handle perfect matches (I think?): give it a 50 if it's perfect
                        if pair.positives >= 20: #changed from == 20; miight be worth keeping in mind for bugs
                            if firstmatch != 0:
                                scorelist.append({"match_score": float(50), "scaffold": scaffold.hit_def,                                                   "hit_location":pair.sbjct_start, "hit_sequence": pair.sbjct, "pam":str(pam.seq),                                                  "match_bars":mmstr})
                                #pamlist.append(pam)
                            firstmatch = 1
                            notfound = 0
                    else:
                        scorelist.append({"match_score": float(0), "scaffold": scaffold.hit_def,                                           "hit_location":pair.sbjct_start, "hit_sequence": pair.sbjct, "pam":str(pam.seq),                                          "match_bars":mmstr})
                        #pamlist.append(pam)
        # Only include the scorelist items that have a score greater than zero
        if return_blast == 1:
            guides[guideindex].annotations["blastdata"] = [s for s in scorelist if s["match_score"] > 0]
        finalscore = int(10000.000 / (100.000 + float(sum(item["match_score"] for item in scorelist))))
        guides[guideindex].annotations["score"] = finalscore
        guides[guideindex].annotations["blast_filesize"] = os.path.getsize(filename + ".blast")
        guides[guideindex].annotations["order_in_xml"] = guideindex
    result_handle.close()
    os.remove(filename) # want to keep BLAST results? comment out these two lines.
    #os.remove(filename + ".blast")
    #return hits
    return guides





def dbconnect(dbname):
    '''
    dbname: the filename of the sqlite db to store guide scores in. Can be an existing database, in which case 
        new entries are appended. If the filename isn't found, a new db with this name is created.
    '''
    # Generate a DictCursor (i.e.
    import sqlite3
    def dict_factory(cursor, row):
        d = {}
        for idx, col in enumerate(cursor.description):
            d[col[0]] = row[idx]
        return d

    class DictConnnection(sqlite3.Connection):
        def __init__(self, *args, **kwargs):
            sqlite3.Connection.__init__(self, *args, **kwargs)

        def cursor(self):
            return DictCursor(self)

    class DictCursor(sqlite3.Cursor):
        def __init__(self, *args, **kwargs):
            sqlite3.Cursor.__init__(self, *args, **kwargs)
            self.row_factory = lambda cur, row: dict_factory(self, row)

    if os.path.isfile(str(dbname + '.db')):
        print("Database " + dbname +".db exists. Appending to this database.")
    # 
    con = sqlite3.connect(str(dbname + '.db'), factory=DictConnnection)
    con.row_factory = dict_factory
    cur = con.cursor()

    try:
        print("Creating tables.")
        cur.execute('''CREATE TABLE scores
                 (sequence text, genome text, score real, version real, filesize real, order_in_xml real, id integer PRIMARY KEY)''')
        cur.execute('''CREATE TABLE locations
                 (id integer, sequence text, genome text, loc_start integer, scaffold test, enzyme text, rel_name text)''')
    except:
        print("Tables already set up.")
    return cur

def


def guide2DB(guide, genome, log_blast = 0, relativeserver="localhost", version=2):
    '''
    A function to enter a SeqRecord (with score) into a mySQL database.
    Need to pass it a DictCursor con.cursor(MySQLdb.cursors.DictCursor) as cur
    Contains an optional "version" tag for if you are trying multiple scoring strategies and
    would like to distinguish between them in the database.
    '''
    # Try to add a new entry into the table for unique scores. Will fail if there's already an entry for same
    # sequence in same genome with same enzyme as the cut site
    con = MySQLdb.connect(relativeserver, 'root', 'pass', 'guidescores');
    cur = con.cursor(MySQLdb.cursors.DictCursor)
    try:
        cur.execute("INSERT INTO scores(sequence, genome, score, version, filesize, order_in_xml) VALUES('{}', '{}', '{}', '{}', '{}', '{}')".format(guide.seq, genome, guide.annotations['score'], version, guide.annotations['blast_filesize'], guide.annotations['order_in_xml']))
        con.commit()
    # except MySQLdb.Error as e:
    #     print("{}".format(e.args))
    except:
        None
    # Pull the specified score record out to grab its ID
    cur.execute("SELECT * FROM scores WHERE sequence = '{}' AND genome = '{}' AND version = '{}'".format(guide.seq, genome, version))
    con.commit()
    data = cur.fetchall()
    # Using its ID, add into locations a reference for this score
    try:
        cur.execute("INSERT INTO locations(id, sequence, genome, loc_start, scaffold, enzyme, rel_name) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}', '{5}', '{6}')".format(data[0]["id"], guide.seq, genome, guide.name, guide.dbxrefs[0], guide.description, guide.id))
        con.commit()
    except MySQLdb.Error, e:
        None
        #print "%s" %e
        #print("location not inserted {} \n".format(guide.id))
    # Insert the annotation blast data if present
    if log_blast ==1:
        try:
            qstring = "INSERT INTO blastdata (id, hit_sequence, match_bars, scaffold, hit_location, match_score, pam, genome) VALUES"
            qstring=list(qstring)
            try:
                for item in guide.annotations["blastdata"]:
                    qstring.extend(str("(" + str(int(data[0]["id"])) + ",\"" + item["hit_sequence"] + "\",\"" + item["match_bars"] + "\",\"" +  item["scaffold"] + "\"," +                      str(item["hit_location"]) + ",\"" +  str(item["match_score"]) + "\",\"" +  item["pam"] + "\",\"" + genome + "\"),"))
                qstring = "".join(qstring)[:-1]
            except e:
                print("%s") %e
            #print qstring
            if len(guide.annotations["blastdata"]) > 0:
                try:
                    cur.execute(qstring)
                    #print "success!"
                except MySQLdb.Error, e:
                    print("%s") %e
            else:
                #print("No BLAST data to enter")
                None
            con.commit()
            #print guide.annotations["blastdata"]
        except:
            #print("No BLAST data found")
            None
    return data

class ProgressBar:
    def __init__(self, iterations):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 50
        self.__update_amount(0)

    def animate(self, iter):
        print('\r', self, end='')
        sys.stdout.flush()
        self.update_iteration(iter + 1)

    def update_iteration(self, elapsed_iter):
        self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
        self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] +             (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)

def dbpresencecheck(guide, genome, con, cur, version="0"):
    try:
        cur.execute("SELECT * FROM scores WHERE sequence = '{}' AND genome = '{}' AND version = '{}'".format(guide.seq, genome, version))
        if len(cur.fetchall()) > 0:
            return 1
        else:
            return 0
    except:
        return 0

def dbpresencecheck_location(guide, genome, con, cur):
    try:
        cur.execute("SELECT * FROM locations WHERE rel_name='{}' AND genome = '{}' AND scaffold='{}'"                    .format(guide.dbxrefs[0], genome, guide.id))
        if len(cur.fetchall()) > 0:
            return 1
        else:
            return 0
    except:
        return 0

class Amplicon:
    '''
    A complete Amplicon for use in EATING; the class contains methods for designing
    amplifying primers and describing the guides within it.

    properties:
    .start, .end: count of runs along chromosome
    .fiveprimeabut, .threeprimeabut: count of guides along chromosome that abut good guides
    .guides: list of good guides in run
    .length: permissible region distance
    .guidecount: number of guides in list (=len(self.guides))
    .genomename: genomename
    .chromosome
    .left_outside: the id of the guide to the left of the outside
    .left_inside: ...
    .permissible_start: absolute numbering vs chromosome
    .required_start_absolute: absolute numbering vs chromosome
    .required_start_relative
    .permissible_region
    '''
    def __init__(self, run):
        self.start = run[0]
        self.end = run[1]
        try:
            self.fiveprimeabut = scores_and_details[self.start-1]
            self.threeprimeabut= scores_and_details[self.end]
        except:
            self.fiveprimeabut = "0"
            self.threeprimeabut = "None"
        self.guides = scores_and_details[self.start:self.end]
        self.length = int(self.threeprimeabut.name)-int(self.fiveprimeabut.name)
        self.guidecount = self.end - self.start
        self.genomename = genomename
        self.chromosome = chromosome


    # Define the sequences necessary
    def pullsequences(self):
        self.left_outside = self.fiveprimeabut.id[-1]
        self.left_inside = self.guides[0].id[-1]
        if self.left_outside == "F" and self.left_inside == "R":
            self.permissible_start = int(self.fiveprimeabut.name) + 10
            self.required_start_absolute = int(self.guides[0].name) +14
        elif self.left_outside == "R" and self.left_inside == "R":
            self.permissible_start = int(self.fiveprimeabut.name) + 1
            self.required_start_absolute = int(self.guides[0].name) +14
        elif self.left_outside == "R" and self.left_inside == "F":
            self.permissible_start = int(self.fiveprimeabut.name) + 1
            self.required_start_absolute = int(self.guides[0].name) +18
        elif self.left_outside == "F" and self.left_inside == "F":
            self.permissible_start = int(self.fiveprimeabut.name) + 10
            self.required_start_absolute = int(self.guides[0].name) +18
        else:
            print("error on left")
        # (fiveprimeabuttingguide, threeprimeabuttingguide), end-start, scores_and_details[start:end]))
        #self.right_inside = item[2][-1][1].id[-1]
        self.right_inside = self.guides[-1].id[-1]
        self.right_outside = self.threeprimeabut.id[-1]
        if self.right_outside == "F" and self.right_inside == "R":
            self.permissible_end = int(self.threeprimeabut.name) + 19
            self.required_end_absolute = int(self.guides[-1].name) + 2
        elif self.right_outside == "R" and self.right_inside == "F":
            self.permissible_end = int(self.threeprimeabut.name) + 10
            self.required_end_absolute = int(self.guides[-1].name) + 8
        elif self.right_outside == "R" and self.right_inside == "R":
            self.permissible_end = int(self.threeprimeabut.name) + 10
            self.required_end_absolute = int(self.guides[-1].name) + 2
        elif self.right_outside == "F" and self.right_inside == "F":
            self.permissible_end = int(self.threeprimeabut.name) + 19
            self.required_end_absolute = int(self.guides[-1].name) + 8
        else:
            print("error on right")

        self.permissible_region = annotchrom[self.permissible_start:self.permissible_end]
        # Bounds that need to be included in PCR product :
        self.required_start_relative = self.required_start_absolute-self.permissible_start
        self.required_end_relative = self.required_end_absolute - self.permissible_start
        #self.amp.dbxrefs=((self.required_start_relative, self.required_end_relative))
    # Set up some other stuff:
        self.permissible_region.name =str(self.fiveprimeabut.name)
        self.permissible_region.id =str(self.fiveprimeabut.name)
        self.permissible_region.description=str(self.guidecount)
        self.permissible_region.seq.alphabet = IUPACAmbiguousDNA()

'''
The following three functions are used to design primers
'''
def al_collect_good_primers(template, primerdict):
    i = 0
    badlist = []
    try:
        while i < primerdict["PRIMER_PAIR_NUM_RETURNED"]:
            bad = 0
            leftprimer = primerdict[str("PRIMER_LEFT_" + str(i) + "_SEQUENCE")]
            leftprimer_start =  str(primerdict[str("PRIMER_LEFT_"+ str(i))].split(",")[0])
            leftprimer_length =  str(primerdict[str("PRIMER_LEFT_"+ str(i))].split(",")[1])
            leftprimer_gc =  str(primerdict[str("PRIMER_LEFT_"+ str(i) + "_GC_PERCENT")])
            leftprimer_tm =  str(primerdict[str("PRIMER_LEFT_"+ str(i) + "_TM")])

            rightprimer = primerdict[str("PRIMER_RIGHT_" + str(i) + "_SEQUENCE")]
            rightprimer_start =  str(primerdict[str("PRIMER_RIGHT_"+ str(i))].split(",")[0])
            rightprimer_length =  str(primerdict[str("PRIMER_RIGHT_"+ str(i))].split(",")[1])
            rightprimer_gc =  str(primerdict[str("PRIMER_RIGHT_"+ str(i) + "_GC_PERCENT")])
            rightprimer_tm =  str(primerdict[str("PRIMER_RIGHT_"+ str(i) + "_TM")])

            product_len = int(rightprimer_start) + int(rightprimer_length) - int(leftprimer_start)
            left_bad = al_screen_primer(leftprimer)
            right_bad = al_screen_primer(rightprimer)
            #print bad
            if left_bad == 0 and right_bad == 0:
                with open("primerlist.txt", "a") as primerlist:
                    primerlist.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                         (template.name,leftprimer,leftprimer_start,leftprimer_length,leftprimer_tm,leftprimer_gc,\
                         rightprimer,rightprimer_start,rightprimer_length,rightprimer_tm,rightprimer_gc,\
                         len(template.seq),str(product_len),template.description))
                    primerlist.close()
                print("Success!")
                break
            if left_bad ==1:
                print("iteration" + str(i) + "left primer" + leftprimer + "is bad")
            if right_bad == 1:
                print("iteration" + str(i) + "right primer" + rightprimer + "is bad")
            i = i +1
            if left_bad ==1 and right_bad ==1 and i ==primerdict["PRIMER_PAIR_NUM_RETURNED"]:
                with open("primerlist.txt", "a") as primerlist:
                    primerlist.write("All the primers were bad for this amplicon!\n")
                    primerlist.close()
    except:
        with open("primerlist.txt", "a") as primerlist:
            primerlist.write("Primer3 failed to find any primers for this amplicon! " + primerdict["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] + "\n")
            primerlist.close()
        print("Primer3 failed to find any primers for this amplicon! " + primerdict["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] + "\n")
        print(sys.exc_info())


def al_screen_primer(primer):
    '''
    Input is a primer as a string.
    '''
    currfile = open("currprimer.fa", "w")
    currfile.write(">" + str(primer) + "\n")
    currfile.write(str(primer))
    currfile.close()
    blastn_cline = NcbiblastnCommandline(query="currprimer.fa", db="xl71", \
    task = "blastn-short",outfmt=5, out="primerblast.tmp", max_target_seqs=100, num_threads = 8)
    blastn_cline
    result = blastn_cline()
    badprimer = 0
    # Parse data
    result_handle = open("primerblast.tmp")
    blast_record = NCBIXML.read(result_handle) # if there were multiple queries, use NCBIXML.parse(result_handle)
    # How many matches are there with more than 14 or matching bases?
    match14 = 0
    for x in blast_record.alignments:
        for y in x.hsps:
            if y.positives > 14:
                match14 = match14 + 1
    match15 = 0
    for x in blast_record.alignments:
        for y in x.hsps:
            if y.positives > 15:
                match15 = match15 + 1
    #print(primer.description)
    #print(match14)
    #print(match15)
    # Set a cutoff of
    if match14 > 40:
        badprimer = 1
    elif match15 > 10:
            badprimer = 1
    return badprimer

def al_primersearch(current_amp):
    '''
    Returns a dict of primer parameters.
    '''
    length_of_required_region = current_amp.dbxrefs[1]-current_amp.dbxrefs[0]
    start_of_required_region = current_amp.dbxrefs[0]
    end_of_required_region = current_amp.dbxrefs[1]
    primeableregionleft_start = str(0)
    primeableregionleft_length = str(start_of_required_region)
    primeableregionright_start = str(end_of_required_region)
    primeableregionright_length = str(len(current_amp)-end_of_required_region)
    boulder = open("current_amp.boulder", "w")
    boulder.write("SEQUENCE_ID=" + current_amp.id + "\n")
    boulder.write("SEQUENCE_TEMPLATE=" + str(current_amp.seq) + "\n")
    #boulder.write("SEQUENCE_INCLUDED_REGION=" + "0," + str(len(current_amp.seq)) + "\n")
    #boulder.write("SEQUENCE_TARGET=" + str(current_amp.dbxrefs[0]) + "," + str(current_amp.dbxrefs[1] - current_amp.dbxrefs[0]) + "\n")
    boulder.write("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=" + primeableregionleft_start + "," + primeableregionleft_length+","\
                  +primeableregionright_start+"," + primeableregionright_length + "\n")
    boulder.write("PRIMER_PRODUCT_SIZE_RANGE=" +str(length_of_required_region) + "-"  + str(len(current_amp)) + "\n")
    boulder.write("PRIMER_PRODUCT_OPT_SIZE=" + str(length_of_required_region) + "\n")
    #boulder.write("P3_FILE_FLAG=1\n")
    boulder.write("=\n")
    boulder.close()
    primer_output = subprocess.check_output(["primer3_core", "current_amp.boulder",\
                              "-p3_settings_file=primer3_global_parameters.txt"])
    primerdict = {}
    for item in primer_output.split("\n")[0:-3]:
        val = item.split("=")[1]
        try:
            val = float(val)
        except: pass
        primerdict[item.split("=")[0]]=val
    return primerdict
