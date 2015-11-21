# Pre-process to align mutations data from SwissProt (humsavar.txt) data with the ProteomeScout dataset. Specifically,
# this generates a pickle file corresponding to the list of unique human accessions which allows for complete description
# of the ProteomeScout dataset which is also cross-referable to the SwissProt mutations dataset. This list is then
# used by the StudyBias code

# Set up the workspace
from proteomeScoutAPI import ProteomeScoutAPI
from pylab import *
import pickle

proteomeScoutFile = '../../data/proteomescout_everything_20151118.tsv'

PTM_API = ProteomeScoutAPI(proteomeScoutFile)
SwissProtMuts = pickle.load(open('swissprot_mutations.p','rb'))

# uniqueKeys returns an arbitrary list of unique keys such that each protein is represented exactly once in the list. Because
# many proteins have multiple IDs, the uniqueKeys list provides a simple way to iterate over all the proteins loaded from 
# the file without worrying about redundancy. However, the actual keys are basically random. We want to correlate our UniProt
# accession IDs from the SwissProtMuts dataset with the ProteomeScout data, but to do that we must ensure each accession in
# the SwissProtMuts dataset corresponds to an entry in the uniqueKeys list. To do this we use the following code segment

# get naieve list of unique keys and extract only those which are human
uniKeys = PTM_API.uniqueKeys
humanIds = []
for key in uniKeys:
    species = PTM_API.get_species(key)
    if species == 'homo sapiens':
        humanIds.append(key)

# Now, for each mutations key
totalMutKeys=len(SwissProtMuts)
count=0

# get a list of all the human IDs
all_human_IDs = []
for key in PTM_API.database:
    species = PTM_API.get_species(key)
    if species == 'homo sapiens':
        all_human_IDs.append(key)


# This is really slow to carry out, but we only need to do it once to construct a definitive list of accessions which
# work with both SwissProt mutations and the ProteomeScout dataset. Basically, this loop will replace all the ProteomeScoutAPI
# uniq IDs with (unique)  SwissProt IDs such that the humanIDs list becomes a list which is a mixture 
for mutKey in SwissProtMuts:
    count=count+1
    if count % 500 == 0:
        print count

    # if the SwissProt mutations derived accession is in the database
    if mutKey in PTM_API.database:
        # then we want to replace whatever key(s) are in the unique_keys lists that points to this record with the 
        # mutations key to ensure consistency between the two databases. 
            
        # get the object ID that all the accessions are point to
        commonID = id(PTM_API.database[mutKey])

        # scan through the whole bloody database looking for other IDs which point to
        # the same object
        for localKey in all_human_IDs:

            # if we found another ID pointing at the object we care about..
            if id(PTM_API.database[localKey]) == commonID:

                # if THAT ID is in the unique keys then THIS is the ID we want
                # to replace from the humanIds list
                if localKey in humanIds:
                    humanIds =[mutKey if x==localKey else x for x in humanIds]
                    continue
    else:
        print "mutKey %s not found in ProteomeScout dataset..." % mutKey

pickle.dump(humanIds, open('parsed_unique_IDs.p', 'wb'))
