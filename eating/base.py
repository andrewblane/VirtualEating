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
import itertools
import multiprocessing
import time
from collections import Counter
import operator
import eating.base

import subprocess
import cStringIO

def digest_target(list_of_targets, process_to_file=0):
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
    if isinstance(list_of_targets, list) == False:
        list_of_targets = [list_of_targets]
    if process_to_file == 1:
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
        return list(itertools.chain.from_iterable(scaffoldcuts)) #the from_iterable part is important here. that took a while.

def score_guides(guides, genome, genomedict, hits=25, max_hsps=100, return_blast =0, sqlitedb="", pam_blast_lookup = 0):
    '''
    To produce a score for a given guide in a given BLAST database. Returns an int, the score.
    Only evaluates a max of 500 hit sequences (i.e. really bad guides will not necessarily have a super accurate score.)
    This version only searches the plus strand of the BLAST db, meant to be used with PAM BLAST DBs generated such that
    all potential guide hits are on the plus strand.
    Arguments:  guide = a SeqRecord object containing a proposed guide sequence (or a list of SeqRecords)
                genome = a BLAST database set up on this machine.
                genomedict = a dict made from a FASTA file identical to the one used to make the BLAST DB. Dict keys should be
                    BLAST db hit_def.
    Returns: List of guides with score in guide.annotations["score"], and if
                return_blast == 1, then blast hits are in guide.annotations["blastdata"]
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

    # Trim the list of guides to those that aren't in the sqlite DB
    if sqlitedb != "":
        guides = [item for item in guides if eating.sqlite.presence_check(item, sqlitedb) == 0]

    # Only do all this if there are unscored guides in the list:
    if len(guides) > 0:
        Bio.SeqIO.write(guides, filename, "fasta")
        # Added dust="no" to stop inflating scores from polyN repeats...
        blastn_cline = NcbiblastnCommandline(query=filename, db=genome, task = "blastn-short",outfmt=5, out=filename + ".blast", max_target_seqs=hits, num_threads = 7, evalue = 10, dust="no", max_hsps = max_hsps)
        #timeit.timeit(blastn_cline, number =1)
        blastn_cline()
        result_handle = open(filename + ".blast")
        blasts = NCBIXML.parse(result_handle)
        # This generates an generator of BLAST objects, each corresponding to a guide query. Loop through those.
        for guideindex, item in enumerate(blasts) :
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
                    if pam_blast_lookup == 0:
                        lookup_context = genomedict[scaffold.hit_def]
                        if end > start:
                            pamstart = end + hit_threeprime_offset
                            pam = lookup_context[pamstart:pamstart+3]
                            #print(pam.seq)
                        elif start > end:
                            pamstart = end - hit_threeprime_offset
                            pam = lookup_context[pamstart-4:pamstart-1].reverse_complement()
                            #print(pam.seq)
                    elif pam_blast_lookup == 1:
                        if end > start:
                            pamstart = end + hit_threeprime_offset
                            context_lookup_command = str("blastdbcmd -db " + genome + " -dbtype nucl -entry " + \
                                                         scaffold.accession + " -range " + str(pamstart+1) + "-" + str(pamstart+3))
                            context_lookup_command_output = subprocess.Popen(context_lookup_command, stdout=subprocess.PIPE, shell = True)
                            fstring = context_lookup_command_output.communicate()
                            fstring = cStringIO.StringIO(fstring[0])
                            pam = SeqIO.read(fstring, "fasta")
                            print(pam.seq)

                        elif start > end: #Get PAMs on the opposite strand
                            pamstart = end - hit_threeprime_offset
                            pamstart = end + hit_threeprime_offset
                            context_lookup_command = str("blastdbcmd -db " + genome + " -dbtype nucl -entry " + \
                                                         scaffold.accession + " -range " + str(pamstart-3) + "-" + str(pamstart-1))
                            context_lookup_command_output = subprocess.Popen(context_lookup_command, stdout=subprocess.PIPE, shell = True)
                            fstring = context_lookup_command_output.communicate()
                            fstring = cStringIO.StringIO(fstring[0])
                            pam = SeqIO.read(fstring, "fasta").reverse_complement()
                            print(pam.seq)


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
                                scorelist.append({"match_score": float(t1 * t2 * t3 * 100), "scaffold": scaffold.hit_def, "hit_location":pair.sbjct_start, "hit_sequence": pair.sbjct, "pam":str(pam.seq),                                              "match_bars":mmstr})
                                #pamlist.append(pam)
                            # Zhang lab algorithm doesn't handle perfect matches: give it a 50 if it's perfect
                            if pair.positives >= 20: #changed from == 20; miight be worth keeping in mind for bugs
                                if firstmatch != 0:
                                    scorelist.append({"match_score": float(50), "scaffold": scaffold.hit_def, "hit_location":pair.sbjct_start, "hit_sequence": pair.sbjct, "pam":str(pam.seq),                                                  "match_bars":mmstr})
                                    #pamlist.append(pam)
                                firstmatch = 1
                                notfound = 0
                        else:
                            scorelist.append({"match_score": float(0), "scaffold": scaffold.hit_def, "hit_location":pair.sbjct_start, "hit_sequence": pair.sbjct, "pam":str(pam.seq),                                          "match_bars":mmstr})
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
        os.remove(filename + ".blast")
    if sqlitedb != "":
        eating.sqlite.db_format(sqlitedb)
        for guide_for_sql in guides:
            try:
                eating.sqlite.put_score(guide_for_sql, sqlitedb, genome)
            except:
                print(str("Failed to enter a guide: " + guide_for_sql.id))
    return guides

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

def count_non_overlapping_guides(guidelist, binding_interference_spacing=20):
    '''
    Sequence position information must be encoded in sequence name attribute,
    e.g. name="100" indicates that left edge of guide (regardless of strand)
    starts 100nt along the target scaffold.

    Optional argument binding_interference_spacing specifies the number of
    nucleotides apart that a guides must be to contribute to unique
    labeling events (for example, of two guides only 10nt apart, it's impossible
    for both to bind at once)

    From the spCas9 crystal structure, it looks like guides will need to be
    at least 24-25 nucleotides apart to bind simultaneously, and possibly
    more.
    '''
    prev_position = 0
    guidecount = 0
    for item in guidelist:
        distance_from_last = int(item.name) - prev_position
        if distance_from_last > binding_interference_spacing:
            guidecount = guidecount + 1
        prev_position = int(item.name)
    return guidecount
