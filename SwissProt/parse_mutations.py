## Python code for reading in a humsvars data file after all the header and footer 
## lines have been removed. 
##
## The function read_swissprot_vars will return a dictionary where the keys are Swissprot
## IDs and each value is a superlist comprising of 1 or more sublists. The sublists each
## have four elements which correspond to
##
## [0] - Original amino acid (1 letter)
## [1] - Mutated amino acid (1 letter)
## [2] - Position in sequence (indexing from 1 onwards - i.e. bioinformatics indexing)
## [3] - The designated status of this mutation. Will be one of 'Disease', 'Polymorphism',
##       or 'Unclassified'
import pickle


THREE_TO_ONE = {'ALA':'A', 
                'CYS':'C',
                'ASP':'D',
                'GLU':'E',
                'PHE':'F',
                'GLY':'G',
                'HIS':'H', 
                'ILE':'I',
                'LYS':'K',
                'LEU':'L',
                'MET':'M',
                'ASN':'N',
                'PRO':'P',
                'GLN':'Q',
                'ARG':'R',
                'SER':'S',
                'THR':'T',
                'VAL':'V',
                'TRP':'W',
                'TYR':'Y'}


def read_swissprot_vars(filename):
    """
    Function which reads in a previously downloaded SwissProt mutations file and extract out
    the mutations and their assigned status

    """

    KWs=set([])

    with open(filename, 'r') as fh:
        content = fh.readlines()
        
    ACC_lookup = {}

    for line in content:
        cleanline=line.strip()
        splitline=cleanline.split()

        
        ACC = splitline[1]
        raw_mut = splitline[3]
        raw_mut = raw_mut.split('p.')[1]


        try:
            org = THREE_TO_ONE[str(raw_mut[0:3]).upper()]
            new = THREE_TO_ONE[str(raw_mut[len(raw_mut)-3:len(raw_mut)]).upper()]
        except KeyError:
            # skip this one..
            continue

        # mutation position
        pos   = int(raw_mut[3:len(raw_mut)-3])
        state =  splitline[4]
        KWs.add(state)

        mutation_record = (pos, org, new, state)

        if ACC not in ACC_lookup:
            ACC_lookup[ACC] = []
        
        ACC_lookup[ACC].append(mutation_record)

    return (ACC_lookup)


mutations_dict = read_swissprot_vars('data_no_header.tsv')
pickle.dump(mutations_dict, open('swissprot_mutations.p', 'wb'))
                    

