# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 12:11:37 2014

@author: andylane

This script is designed to loop through generated BLAST files and given a 
particular off-target-match cutoff, plot out the specificity profiles by scaffold


"""

"""
First, load in the BLAST files one by one and make a table that summarizes hits:
FOR EACH BLAST FILE
    FOR EACH SEQUENCE(TARGET) IN EACH BLAST FILE
        GET THE TOTAL NUMBER OF HITS
        GET THE HIT LOCATIONS (TUPLE)
        ADD THIS INFO TO ONE_BIG_FILE AS FOLLOWS
            (BLAST FILE, SEQUENCE, [(HIT SCAFFOLD, HIT START, HIT END), (HIT SCAFFOLD, HIT START, HIT END), (HIT SCAFFOLD, HIT START, HIT END)])
        
SEPARATELY:
        IF TOTAL LEN(TUPLE(3)) IS GREATER THAN SOME THRESHOLD:
            TRIM THE LIST INTO FILTERED_LIST

FROM FILTERED_LIST:
    FOR EACH SEQUENCE
            HASH (HIT SCAFFOLD, HIT START, HIT END) AND REMOVE DUPES
    MAKE A FLAT FILE WITH ALL THE (HIT SCAFFOLD, HIT START, HIT END) TUPLES
    COUNT THE SCAFFOLD OCCURRENCES OF EACH TARGET(DICTIONARY)
    PLOT OCCURRENCES PER SCAFFOLD
    
"""
from Bio.Blast.Applications import NcbiblastnCommandline
import Bio 
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Restriction 
from Bio.Restriction import *
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seprepseq

fasta_file = "../../Genomic Data/LAEVIS_7.1.repeatMasked.fa"  # Input fasta file
longestseqo = extractlongest.extractlongest(fasta_file)

# Write out the longest sequence
output_handle = open("longestseq.fasta", "w")
Bio.SeqIO.write(longestseqo, output_handle, "fasta")
output_handle.close()
""" 
This is what I ran 6/21 on server
"""
# Read it in again
longestseq = []
for record in Bio.SeqIO.parse(open("longestseq.fasta", "rU"), "fasta"):
    longestseq.append(record)
longestseq = longestseq[0]

nonrepseqs = seprepseq.seprepseq(longestseq)

#Define a list of sequences to exclude from the top144 (i.e ones you already have oligos for)
#excludelist = ["12709050", "8261124", "5570943", "14351343", "8489316", "9593441", "19747367", "14421422"]

# Sort nonrepseqs by length, take longest 144
sortednonrepseqs = []
for record in nonrepseqs:
    #if record.id not in excludelist:
    sortednonrepseqs.append(record)
        
#sortednonrepseqs = (record for record in longestseq if record.id  in excludelist)
sortednonrepseqs.sort(cmp=lambda x,y: cmp(len(y),len(x)))
amps144 = sortednonrepseqs[0:300] # now 300 instead of 144...

# Write out a FASTA file of amps144
records = amps144
SeqIO.write(records, "amps144.fasta", "fasta")

"""
Next part: cut up 300 amps
"""
amps_300 = []
master_featlist = []

for record in Bio.SeqIO.parse("amps144.fasta", "fasta"):
    record.seq.alphabet = IUPACAmbiguousDNA()
    amps144_300.append(record)

for individual_sequence in amps_300: # Use the 300 amps that has ids matching the BLAST file names.
    matchlist= []
    offtgt = []
    title = ""
    blast_records = []
    
    result_handle = open("blastresults_"+ str(individual_sequence.id) + ".blast")
    blast_records = NCBIXML.parse(result_handle) # use NCBIXML.parse(result_handle) for multiple queries here
    
    twentymers = seprepseq.maketwentymers(individual_sequence)     # Re-generating cut 20mers only because this is the only source of info for locations of blast query results    
    twentymers_meta_tuples = [(item.id+" "+item.description, int(item.name)+int(item.dbxrefs[0]), item.seq) for item in twentymers]
    
    
    blast_records_list = [] # Harness the generator to make a list...
    for blast_record in blast_records:
        blast_records_list.append(blast_record)     
    result_handle.close()
    
    annotated_blast_records = []
    for blast_record in blast_records_list:
        for y in twentymers_meta_tuples:  
            if blast_record.query in y[0]:
                annotated_blast_records.append([blast_record, y])

        
    # Next: get necessary info from blast_records_list. (total hits per seq, hit locations)    
    for blast_record in annotated_blast_records:
        featlist = []        
        for x in blast_record[0].alignments:
            for y in x.hsps:
                    featlist.append([individual_sequence.id, x.title, y.sbjct_start, y.sbjct_end, y.positives])
        blast_record.append(featlist)

    master_featlist.append(annotated_blast_records[:])
    print(individual_sequence.id)
#master_featlist_file = open("master_featlist.csv", "w")
#for item in master_featlist:
#    master_featlist_file.write("{0} \n".format(item))
#master_featlist_file.close()

import cPickle as pickle
with open('master_featlist_file.pkl', 'wb') as output:
    pickle.dump(master_featlist, output, -1)
output.close()

with open('master_featlist_file.pkl', 'rb') as input:
    master_featlist = pickle.load(input)

"""
This probably isn't useful now.
    # Figure out the target yield as a function of tolerance of mismatches  
    specificity_vs_tolerance = []
    for i in range(300)[1:]:
        specificity_vs_tolerance.append(len(blast_records_list) - len([item[0:3] for item in matchlist if item[5] == i])) # Returns the number of targets with more than i hits in the genome
        
    # Can you figure out the lengths of the scaffolds, badness plot and dynamically choose a cutoff?
    offtgt20 = [item[0:3] for item in matchlist if item[5] == 20]
    offtgt25 = [item[0:3] for item in matchlist if item[5] == 25]
    offtgt50 = [item[0:3] for item in matchlist if item[5] == 50]
    offtgt100 = [item[0:3] for item in matchlist if item[5] == 100]

    # Make a database matching target to its own off-target hits list:
    ## The NAME attribute of the twentymers list that was the input to FILE 
    ## would have been useful here, but it was lost when output to FASTA and sent to BLAST
        
    import string
    matchlist_with_tgt_abs_locations = matchlist[:]    
    for item in matchlist_with_tgt_abs_locations:
        
        
        item_relative_loc_within_amp = string.find(str(item[2]), str(individual_sequence.seq))
        itemloc = tuple(int(item[0]) + int(string.find(str(item[2]), str(individual_sequence.seq))))
        item = item + tuple(int(item[0]) + int(string.find(str(item[2]), str(individual_sequence.seq))))



"""    
"""
Guess I need two things: for each cutoff, figure out what the targeting plot looks like (can we make a metric?) 
AND figure out the number of targetable regions on longscaffold after masking.

Going to need to: 
- generate a cut=id (and sequence) scaffold hits list to refer to later: DONE - this is master_featlist

- for each offtgt amp of a cutoff:
    mask on offtgts -DONE (amps144_masked_cutoff)
    generate new amps, top 144 -DONE

- for new amps:
    re-cut
        make cut ids that can be matched to already-made cut id scaffold hits list
    do specificity analysis()

"""

#cutoffsrange = range(5, 201, 10) # Testing a range from 5 to 100 in steps of 5
cutoffsrange = []
expgenerator = (2**exp for exp in range(1, 11))
for n in expgenerator:
    print(n)
    cutoffsrange.append(n)

replstr = "N" * 20
targets_per_threshold = []
amps144_300_current_cutoff = []
for cutoff in cutoffsrange: # generate a masked amps144 for a given offtarget threshold
    amps_144_300_current_cutoff_masked = []
    offtgt_summary = []    
    for superitem in master_featlist: #because it's split into the original amps144 amps at a top level
        for twmer in superitem: # first, pick out the targets that are below threshold
            n_target_hits = len([subitem for subitem in [item[4] for item in twmer[-1]] if subitem > 17]) # this gets the count of hits for each target (over 18 nt matches)
            if n_target_hits > cutoff:
                # format of offtarget_summary: item, originating amp, position on originating amp
                offtgt_summary.append(((twmer[1:-1], int(twmer[-1:][0][0][0]), int(twmer[1][1]) - int(twmer[-1:][0][0][0])))) # if the count of hits is greater than the cutoff, add the item to a new list for masking. 
    # Find the locations of poorly-specific seqs and mask them
    amps144_300_current_cutoff = amps144_300[:]
    for amp in amps144_300_current_cutoff:
        amp = SeqRecord(seq=amp.seq.tomutable(), id = amp.id)
        for tgt in offtgt_summary: 
                if int(amp.id) == tgt[-2]:
                    amp.seq[int(tgt[-1])-1:int(tgt[-1])+19] = replstr
        amps_144_300_current_cutoff_masked.append(amp)
        
    with open('amps144_300_masked_cutoff'+str(cutoff)+'.pkl', 'wb') as output:
        pickle.dump(amps_144_300_current_cutoff_masked, output, -1)
    output.close()
    
    # next, use the items in cutoff_summary to generate individual sequences out amps144 using seprepseqN()
    amps144_300_masked_cutoff_splitonNs = []
    for item in amps_144_300_current_cutoff_masked:
        split = seprepseq.seprepseqN(item) #altered seprepseq to mask on Ns instead of lwrcase
        for x in split:
            x.id = int(item.id) + int(x.id) # haven't totally checked that ids are base-pair accurate, but they should be close
        amps144_300_masked_cutoff_splitonNs.append(split)    
    amps144_300_masked_cutoff_splitonNs = [item for sublist in amps144_300_masked_cutoff_splitonNs for item in sublist] # Flatten the nested list

    amps144_300_masked_cutoff_splitonNs.sort(cmp=lambda x,y: cmp(len(y),len(x))) # Generate a sort... somehow
    #specific_amps144 = offtgt_mask_split_results[0:144] # Don't just pick longest
    # Here, we find the target density for each of 300 amps and sort by that:
    for item in amps144_300_masked_cutoff_splitonNs:
        cuts = [ScrFI.search(item.seq)]
        cuts.append(list(HpaII.search(item.seq)))
        cuts.append(list(BfaI.search(item.seq)))
        cuts = [t for sublist in cuts for t in sublist]
        cuts.sort()
        
        # This seems to be good code for removing dupes from a list, or nearby numbers from a list
        # Here, it's being used to prevent dupes from increasing apparent target density on amps
        result = []
        num= 0
        for i in cuts:
            if abs(i-num)<20:
                    num=i
            else:
                result.append(i)
                num=i
        cuts = result[:]
        amp_tgt_count = len(cuts) * 2
        item.description = amp_tgt_count
    
    amps144_300_masked_cutoff_splitonNs = sorted(amps144_300_masked_cutoff_splitonNs, key=lambda x: x.description)
    specific_amps144 = amps144_300_masked_cutoff_splitonNs[-144:]
    
    # Write out the best amps
    for item in specific_amps144: #gb needs the ids to be strings...
        item.id=str(item.id)
        item.description=str(item.description)

    # specific_amps144 is now ready to be written as gb...
    
    # This part just determines the number of targets in a given set of SeqRecords    
    total_strategy_targets = 0
    for item in specific_amps144:
        cuts = [ScrFI.search(item.seq)]
        cuts.append(list(HpaII.search(item.seq)))
        cuts.append(list(BfaI.search(item.seq)))
        cuts = [t for sublist in cuts for t in sublist]
        cuts.sort()
        
        # This seems to be good code for removing dupes from a list, or nearby numbers from a list
        result = []
        num= 0
        for i in cuts:
            if abs(i-num)<3:
                    num=i
            else:
                result.append(i)
                num=i
        cuts = result[:]
        amp_tgt_count = len(cuts) * 2
        total_strategy_targets = total_strategy_targets + amp_tgt_count
    print(str("There will be {0} total targets with this set of amplicons".format(total_strategy_targets)))
    targets_per_threshold.append((cutoff, total_strategy_targets)) #this doesn't take into account specificity, just target yield

    with open('secondarily_cut_amps_using_cutoff_'+str(cutoff)+'.pkl', 'wb') as output:
        pickle.dump(specific_amps144, output, -1)
    output.close()

targets_per_threshold
    

"""
END CODE DESIGNED TO DETERMINE THE TOTAL LABELING AT EACH CUTOFF, WHEN CONSTRAINED TO 144 AMPS
"""


"""
- for new amps:
    re-cut
        make cut ids that can be matched to already-made cut id scaffold hits list
    do specificity analysis() - almost done!
    next(weds): get the hits dictionary for all amps. there's something wrong with the hits dict now, though, 
    so it seems to have 50% of the hits be a single locus... probably writing out the hit name wrong. this will have upstream effects...
    or maybe it was just that there were a lot of BfaI sites that were all the same, and they come out to hit the same spot?
    so, this is going to under-count when a locus repeats: is it only going to get the first place the  hit? it shouldn't...
"""
from Bio.SeqFeature import SeqFeature, FeatureLocation

for cutoff in cutoffsrange: ##pick up from here; make loop of below for all iters
    with open('secondarily_cut_amps_using_cutoff_'+str(cutoff)+'.pkl', 'rb') as input:
        specific_amps144 = pickle.load(input)
    longhitlist = []
    for item in specific_amps144:
        paired_cuts_list =[]
        # Remove duplicate BLAST hits, since no one BLAST locus can be labeled twice
        # First, generate the string to be checked for duplicates (Scaffold, start_match, end_match)
        cuts = []
        if len(ScrFI.search(item.seq)) > 0:
            cuts.append(list((ScrFI.search(item.seq), "ScrFI")))
        if len(HpaII.search(item.seq)) > 0:
            cuts.append(list((HpaII.search(item.seq), "HpaII")))
        if len(BfaI.search(item.seq)) > 0:
            cuts.append(list((BfaI.search(item.seq), "BfaI")))
        #cuts = [t for sublist in cuts for t in sublist]
        for c in cuts:
            for u in c[0]:
                paired_cuts_list.append((u, c[-1]))        
        paired_cuts_list.sort()

        # This seems to be good code for removing dupes from a list, or nearby numbers from a list
        # Next: have cut sites in amps. Using the correct adjustments, add them as seqfeatures
        absolute_locations_adjusted = []    
        i = 0
        for m in paired_cuts_list:
            if i < len(paired_cuts_list)-1:
                if m[1] == "ScrFI":
                    #offset here (e.g. 1:21 is for simulating MBN digeston)
                    # because the entry.names are rooted on the right of each fragment, the length
                    # of the entry.name has to be subtracted to get the desired left position for the "reverse"
                    # tgts
                    absolute_locations_adjusted.append(((int(item.id) + 1 + m[0], m[1]), "rev"))
                else: # Should work for HpaII/BfaI
                    absolute_locations_adjusted.append(((int(item.id) + 2 + m[0], m[1]), "rev"))
                absolute_locations_adjusted.append(((int(item.id) + int(paired_cuts_list[i+1][0]) - 20, m[1]), "fwd")) # amp start + next cut point(i+1) - 2
        i = i + 1
        
        #    hitlist = []
        #    for item in absolute_locations_adjusted:
        #        for blast_result in master_featlist:
        #            for subblast in blast_result:            
        #                if subblast[1][1] == item[0][0]:
        #                    scaffold_hits = [y for y in subblast[-1] if y[-1] > 17]
        #                    hitlist.append((item, \
        #                    len([subitem for subitem in \
        #                    [z[4] for z in subblast[-1]] if subitem > 17]), \
        #                    [(x[1][:-19], x[2], x[3], x[4]) for x in scaffold_hits]))
    
        for v in absolute_locations_adjusted:
            for blast_result in master_featlist:
                for subblast in blast_result:            
                    if subblast[1][1] == v[0][0]:
                        scaffold_hits = [y for y in subblast[-1] if y[-1] > 17]
                        hitlist.append((v, \
                        len([subitem for subitem in \
                        [z[4] for z in subblast[-1]] if subitem > 17]), \
                        [(x[1][:-19], x[2], x[3], x[4]) for x in scaffold_hits]))
    
    
        # Ok, so: have a hitlist with info for each 20mer on its enzyme site,  
        # direction, number of 18+mer hits and the scaffold locations of all hits
        
        # Next: extract scaffolds only and deduplicate    
        hitlist_flat = [k for k in [h[-1] for h in [j for j in hitlist]]]
        from itertools import chain
        hitlist_flat = list(chain(*hitlist_flat))
        
        longhitlist.append(list(zip(hitlist_flat, [0] * len(hitlist_flat))))
        
    longhitlist_flat = list(chain(*longhitlist))
    master_hitlist_dedup = dict(longhitlist_flat)
    x = 0
    for x in longhitlist_flat: # This increments the hitlist_dict for each item in the flattened hitlist
        master_hitlist_dedup[x[0]] += 1
    
    # Have: list of hits per locus. Don't really care about that because each locus can't be labeled twice, but it had to be deduplicated in this way.
    # Next, need to count the occurrence of each scaffold in the list
    deduped_on_scaffold = master_hitlist_dedup.copy()
    master_hitlist_dedup_list = deduped_on_scaffold.items()
    master_hitlist_dedup_list = [(x[0][0], 0) for x in master_hitlist_dedup_list] #Need instances of scaffold...
    
    deduped_on_scaffold = dict(master_hitlist_dedup_list)
    for x in master_hitlist_dedup_list: # This increments the hitlist_dict for each item in the flattened hitlist
        deduped_on_scaffold[x[0]] = deduped_on_scaffold[x[0]] + 1
    
    #Write out the match hits; remember to call matchlistfile.close() or random shit happens!
    deduped_on_scaffold_list = deduped_on_scaffold.items() #otherwise just the key, not the value, gets output
    deduped_on_scaffold_file = open("deduped_on_scaffold_"+str(cutoff)+".csv", "w")
    for d in deduped_on_scaffold_list:
        deduped_on_scaffold_file.write("{0} \n".format(d))
    deduped_on_scaffold_file.close()
    
    # Next up: collate the info (remove annoying parentheses etc)
    import csv
    scaffhits = []
    with open("deduped_on_scaffold_"+str(cutoff)+".csv",'rb') as csvin:
        csvin = csv.reader(csvin, delimiter=',')
        for row in csvin:
            scaff = row
            scaffhits.append(scaff)
    
    for n in scaffhits:
        n[0] = n[0][3:-1]
        n[1] = int(n[1][1:-2])
        n = str(n[0] + "," + str(n[1]))
    deduped_on_scaffold_file = open("deduped_on_scaffold_"+str(cutoff)+".csv", "w")
    for d in scaffhits:
        deduped_on_scaffold_file.write("{0} \n".format(d))
    deduped_on_scaffold_file.close()


"""
END CODE RELATED TO *ANALYSIS* OF OFF TARGET HITS (RATHER THAN REFINEMENT OF AMPLICONS)
"""

