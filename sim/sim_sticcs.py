import msprime
import tskit
import tsinfer
import demes

import numpy as np

from sticcs import sticcs
from twisst2 import TopologySummary
from twisst2 import parse_newick

import time, sys
import random
import itertools
import math
from multiprocessing import Pool
import subprocess
import tempfile
import gzip
import os

def n_choose_k(n, k):
    f = math.factorial
    return f(n) // f(k) // f(n-k)

def show_tskit_tree(tree, node_labels=None): print(tree.draw(format="unicode", node_labels=node_labels))


def filter_to_infinite_sites(ts):
    tables = ts.dump_tables()
    
    # Group mutations by site
    site_mutations = {}
    for mut in tables.mutations:
        site_id = mut.site
        if site_id not in site_mutations:
            site_mutations[site_id] = []
        site_mutations[site_id].append(mut)
    
    # Clear mutations table and rebuild with one mutation per site
    tables.mutations.clear()
    
    for site_id, mutations in site_mutations.items():
        if len(mutations) == 1:
            # Single mutation - keep it
            tables.mutations.append(mutations[0])
        else:
            # Multiple mutations - keep one (e.g., earliest, or random choice)
            chosen_mut = max(mutations, key=lambda x: x.time)  # Keep earliest
            tables.mutations.append(chosen_mut)
    
    return tables.tree_sequence()

def two_pop_model(N=1e5, t=1e5):   #population sizes
    
    #use the deme builder to set up the demographic history
    b = demes.Builder(time_units="generations")
    b.add_deme("anc", epochs=[{"start_size":N, "end_time":t}])
    b.add_deme("A",   ancestors=["anc"],  epochs=[{"start_size":N}])
    b.add_deme("B",   ancestors=["anc"],  epochs=[{"start_size":N}])
    
    graph = b.resolve()
    
    return(graph)

def two_pop_OG_model(N=1e5, tAB=0, tABC=10e5):   #population sizes
    
    #use the deme builder to set up the demographic history
    b = demes.Builder(time_units="generations")
    b.add_deme("ABC", epochs=[{"start_size":N, "end_time":tABC}])
    b.add_deme("AB", ancestors=["ABC"], epochs=[{"start_size":N, "end_time":tAB}])
    b.add_deme("C", ancestors=["ABC"], epochs=[{"start_size":N}])
    b.add_deme("A",   ancestors=["AB"],  epochs=[{"start_size":N}])
    b.add_deme("B",   ancestors=["AB"],  epochs=[{"start_size":N}])
    
    graph = b.resolve()
    
    return(graph)


def run_tsinfer(haps, positions, seq_len):
    
    sample_data = tsinfer.SampleData(sequence_length=seq_len)
    
    for i in range(haps.shape[0]):
        #genotyps for this site
        gts = np.array(haps[i,], dtype=np.int8)
        n_alleles = max(gts)+1
        #alleles for this site
        alleles = ["A", "T", "C", "G"][:n_alleles]
        if -1 in gts:
            #if there is missing data, add that
            gts[gts==-1] = n_alleles
            alleles += ["-1"]
        
        sample_data.add_site(positions[i], gts, alleles) # note that for tsinfer versions from 0.3.0 onwards, you can add 'ancestral_allele=0', but that's irrelevant if it's always the first allele)
    
    sample_data.finalise()
    
    #run tsinfer
    ts_tsinfer = tsinfer.infer(sample_data)
    
    #simplify is needed for kc_distance to work
    return ts_tsinfer.simplify(reduce_to_site_topology=True)


def run_sticcs(haps, positions, seq_len, forced_ploidy=1, second_chances=False, silent=False):
    # In all cases we simulate haploids. These can be represented as diploids or tetraploids by manipulating the data
    # This is preferred over simulating different ploidies in msprime, because that changes the effective population size
    # Our way ensures that the same exact simulation can be analysed at different ploidy levels without changing any population parameters
    # Make sure that n is always a multiple of the ploidy you want to impose on it.
    
    n = haps.shape[1]
    
    #check that you have a multiple of forced ploidy haplotypes to work with!
    assert n % forced_ploidy == 0
    
    ploidies = np.array([forced_ploidy]*int(n/forced_ploidy))
    
    der_counts, der_positions = sticcs.get_dac_from_haploid_matrix(haps, positions, ploidies)
    
    #currently doesn't accept missing data, so remove any sites with missing data
    no_missing = np.where(der_counts.min(axis=1) >= 0)
    
    der_counts = der_counts[no_missing]
    der_positions = der_positions[no_missing]
    
    #sticcs has three steps - find patterns, make clusters, infer trees.
    #These could theoretically be combined into a single finctions within sticcs, but currently I prefer the transparency
    patterns, matches, n_matches = sticcs.get_patterns_and_matches(der_counts)
    
    clusters = sticcs.get_clusters(patterns, matches, der_positions, ploidies=ploidies,
                                    second_chances=second_chances,
                                    seq_start=1, seq_len=seq_len, silent=silent)
    
    ts_sticcs = sticcs.infer_ts(patterns, ploidies, clusters, silent=silent)
    
    #simplify is needed for kc_distance to work
    return ts_sticcs #.simplify(reduce_to_site_topology=True)


# a class that just holds a unch of sticcs trees, but has a few attributes and methods that resemble a treesequence
# This allos us to import newick trees and treet the list as a tree sequence for a few things (like my quartet distance calculation)
class TreeList:
    def __init__(self, trees):
        self.num_trees = len(trees)
        self.tree_list = trees
    
    def trees(self):
        for tree in self.tree_list:
            yield tree


def run_argweaver(haps, positions, N, rec, mu, sequence_start, sequence_end, burnin=500, sample_n_args=10, spacing=10):
    n_sites, n_samples = haps.shape
    iters = burnin + sample_n_args*spacing
    
    #tmpdir = tempfile.mkdtemp(dir=".") #for debugging
    
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        
        sites_file=f"{tmpdir}/temp.sites"
        region_name = "chr1"
        
        # Write .sites file in ARGweaver format
        with open(sites_file, "w") as f:
            # Header line with sample names
            sample_names = [f"n{i}" for i in range(n_samples)]
            header = "NAMES\t" + "\t".join(sample_names)
            f.write(header + "\n")
            
            # Region line
            region_line = f"REGION\t{region_name}\t{sequence_start}\t{sequence_end}"
            f.write(region_line + "\n")
            
            # Site data - one line per segregating site
            for site_idx, pos in enumerate(positions):
                
                # Get alleles for this site
                alleles = []
                for sample_idx in range(n_samples):
                    allele_idx = haps[site_idx, sample_idx]
                    alleles.append("ACGT"[allele_idx])
                
                # Write position and alleles
                allele_string = "".join(alleles)
                f.write(f"{pos}\t{allele_string}\n")
        
        #here the command for running argweaver
        cmd = ["arg-sample", f"-s {sites_file}", f"-N {N/2}", f"-r {rec}", f"-m {mu}", f"--maxtime {N*10}", f"--iters {iters}", f"--sample-step {spacing}", f"-o {tmpdir}/out"]
        
        result = subprocess.run(" ".join(cmd), capture_output=True, text=True, check=True, shell=True)
        
        #1000 reps with 100 samples
        #sample the last 50 (last 500 reps), so that the first 500 reps is burn-in
        #read stats file to find best arg
        sample_indices = list(range(burnin, burnin+sample_n_args*spacing-1, spacing))
        
        #return sample_n_args args in a list
        tree_lists = []
        for idx in sample_indices:
            trees = []
            with gzip.open(f"{tmpdir}/out.{idx}.smc.gz", "rt") as argfile:
                names = argfile.readline().split()[1:]
                leaf_idx_dict = dict([(str(i), int(names[i].lstrip("n"))) for i in range(n_samples)])
                chom, chrom_start, chrom_len = argfile.readline().split()[1:]
                for line in argfile:
                    if line.startswith("TREE"):
                        elements = line.split()
                        interval=(int(elements[1]), int(elements[2]),)
                        tree_newick = elements[3]
                        tree, node_labels = parse_newick.newick_to_sticcs_Tree(tree_newick, leaf_idx_dict=leaf_idx_dict, interval=interval)
                        trees.append(tree)
            tree_lists.append(TreeList(trees))
    
    return tree_lists


def run_singer(haps, positions, N, rec, mu, sequence_start, sequence_end, singer_path, burnin=500, sample_n_args=10, spacing=20):
    
    #tmpdir = tempfile.mkdtemp(dir=".") #for debugging
    
    assert burnin % spacing ==0, "Burnin must be a multiple of the spacing step"
    
    first_to_keep = burnin//spacing
    total_args = first_to_keep + sample_n_args
    
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
    
        #write vcf
        with open(f"{tmpdir}/temp.vcf", "wt") as vcf:
            haplotypes_to_diploid_vcf(haps, positions, contig_length=sequence_end, output_file=vcf)
        
        singer_cmd = [singer_path + "singer_master", f"-m {mutation_rate}", f"-Ne {N/2}", f"-n {total_args}", f"-thin {spacing}", "-polar 0.99", f"-v {tmpdir}/temp", f"-output {tmpdir}/temp", "-start 0", f"-end {int(sequence_end)}"]
        #" ".join(singer_cmd) #for debugging
        
        result = subprocess.run(" ".join(singer_cmd), capture_output=True, text=True, check=True, shell=True)
        
        convert_cmd = [singer_path + "convert_to_tskit", f"-input {tmpdir}/temp", f"-output {tmpdir}/temp", f"-start {first_to_keep}", f"-end {total_args}", "-step 1"]
        #" ".join(convert_cmd) #for debugging
        
        result = subprocess.run(" ".join(convert_cmd), capture_output=True, text=True, check=True, shell=True)
        
        ts_list = [tskit.load(f"{tmpdir}/temp_{i}.trees") for i in range(first_to_keep, total_args)]
        
        return ts_list




def run_relate(haps, positions, N, mu, rec, relate_bin):
    
    #only zero and 1 are accepted
    if haps.max() > 1:
        n_alleles = np.apply_along_axis(lambda x: np.count_nonzero(np.bincount(x, minlength=4)), 1, haps)
        
        biallelic = n_alleles == 2
        
        haps = haps[biallelic,:]
        positions = positions[biallelic]
        
        haps[haps > 1] = 1
    
    n_variants, n_samples = haps.shape
    
    sample_ids = [f"n{i}" for i in range(n_samples)]
    
    #tmpdir = tempfile.mkdtemp(dir=".") #for debugging because this temdir won't be automatically deleted
    #if True:
    
    #Create temporary directory
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        
        relate_prefix = tmpdir.split("/")[-1]
        
        # File paths
        haps_file = f"{tmpdir}/temp.haps"
        sample_file = f"{tmpdir}/temp.sample"
        map_file = f"{tmpdir}/temp.map"
        
        with open(haps_file, 'w') as f:
            for i, pos in enumerate(positions):
                
                # Write basic info: chr, snp_id, position, allele0, allele1
                line = f"chr1 {pos} {pos} A T "
                
                # add haps
                line = line + " ".join(haps[i, ].astype(str))
                
                f.write(line + "\n")
        
        # Write .sample file
        with open(sample_file, 'w') as f:
            # Header
            f.write("ID_1 ID_2 missing\n")
            f.write("0 0 0\n")
            
            # Sample information
            for sample_id in sample_ids:
                f.write(f"{sample_id} NA 0\n")
        
        # write .map file (assumes uniform map)
        cM_pos = 100 * rec * positions
        rate_cM = rec*1e8
        
        with open(map_file, 'w') as f:
            for i in range(n_variants):
                f.write(f"{positions[i]} {rate_cM} {cM_pos[i]}\n")
        
        # For now, let's use --dist to specify that positions are in base pairs
        relate_cmd = [
            f"{relate_bin}/Relate",
            "--mode", "All",
            "-m", str(mu),
            "-N", str(N),
            "--haps", haps_file,
            "--sample", sample_file,
            "--map", map_file,
            "-o", f"{relate_prefix}_relate"
        ]
        
        # Run Relate
        #print(f"Running Relate with command: {' '.join(relate_cmd)}")
        result = subprocess.run(relate_cmd, capture_output=True, text=True, check=False)
        
        #convert to treesequence
        convert_cmd = [f"{relate_bin}/RelateFileFormats", "--mode", "ConvertToTreeSequence",
                       "-i", f"{relate_prefix}_relate", "-o", f"{relate_prefix}_relate"]
        
        #print(f"Converting with command: {' '.join(convert_cmd)}")
        result = subprocess.run(convert_cmd, capture_output=True, text=True, check=True)
        
        ts_relate = tskit.load(f"{relate_prefix}_relate.trees")
        
        os.rmdir(f"{relate_prefix}_relate/chunk_0/paint")
        os.rmdir(f"{relate_prefix}_relate/chunk_0/")
        os.rmdir(f"{relate_prefix}_relate/")
        os.remove(f"{relate_prefix}_relate.anc")
        os.remove(f"{relate_prefix}_relate.mut")
        os.remove(f"{relate_prefix}_relate.trees")
        
        return ts_relate

def haplotypes_to_diploid_vcf(haplotypes, positions, contig_length, output_file):
    # Validate inputs
    if haplotypes.shape[1] % 2 != 0:
        raise ValueError("Number of haplotypes must be even")
    
    if len(positions) != haplotypes.shape[0]:
        raise ValueError("Length of positions must match number of sites")
    
    # Mapping from numeric encoding to nucleotides
    allele_map = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    
    n_sites, n_haplotypes = haplotypes.shape
    n_individuals = n_haplotypes // 2
    
    # Generate individual IDs
    individual_ids = [f'n{i}' for i in range(n_individuals)]
    
    # Write VCF header directly to file
    output_file.write('##fileformat=VCFv4.2\n')
    output_file.write('##reference=unknown\n')
    output_file.write('##contig=<ID=chr1,length={}>\n'.format(contig_length))
    output_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    
    # Write column header
    header_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + individual_ids
    output_file.write('\t'.join(header_cols) + '\n')
    
    # Process each site
    for site_idx in range(n_sites):
        site_haplotypes = haplotypes[site_idx, :]
        
        # Get unique alleles at this site
        unique_alleles = np.unique(site_haplotypes)
        alleles = [allele_map[allele] for allele in sorted(unique_alleles)]
        
        # ALT alleles are everything except A
        alt_alleles = [allele for allele in alleles if allele != 'A']
        alt_string = ','.join(alt_alleles) if alt_alleles else '.'
        
        # Create allele to index mapping for genotype encoding
        allele_to_idx = {'A': 0}
        idx_counter = 1
        for allele in alt_alleles:
            allele_to_idx[allele] = idx_counter
            idx_counter += 1
        
        # Generate genotypes for each individual
        genotypes = []
        for ind_idx in range(n_individuals):
            hap1_idx = ind_idx * 2
            hap2_idx = ind_idx * 2 + 1
            
            hap1_allele = allele_map[site_haplotypes[hap1_idx]]
            hap2_allele = allele_map[site_haplotypes[hap2_idx]]
            
            hap1_gt = allele_to_idx[hap1_allele]
            hap2_gt = allele_to_idx[hap2_allele]
            
            # Unphased genotype
            genotype = f'{hap1_gt}|{hap2_gt}'
            genotypes.append(genotype)
        
        # Build VCF line
        vcf_line = [
            'chr1',  # CHROM
            str(int(positions[site_idx])),  # POS
            '.',  # ID
            "A",  # REF (always A)
            alt_string,  # ALT
            '.',  # QUAL
            'PASS',  # FILTER
            f'AA=A',  # INFO - Ancestral allele (A is always ancestral)
            'GT'  # FORMAT
        ] + genotypes
        
        # Write VCF line directly to file
        output_file.write('\t'.join(vcf_line) + '\n')


#def get_kc_dist(tree1, tree2, tree2_leaf_reorder=None):
    
    #leaves1 = range(tree1.num_samples())
    #leaves2 = tree2_leaf_reorder if tree2_leaf_reorder is not None else range(tree2.num_samples())
    
    #dists_to_root_1 = []
    #dists_to_root_2 = []
    
    ##tree 1 distances to root
    #for pair in itertools.combinations(leaves1, 2):
        #node = tree1.mrca(*pair)
        #dist = 0
        #while node != tree1.root:
            #dist += 1
            #node = tree1.parent(node)
        
        #dists_to_root_1.append(dist)
    
    ##tree2 distances to root
    #for pair in itertools.combinations(leaves2, 2):
        #node = tree2.mrca(*pair)
        #dist = 0
        #while node != tree2.root:
            #dist += 1
            #node = tree2.parent(node)
        
        #dists_to_root_2.append(dist)
    
    ##euclidean distance (apparently this is correct)
    #return np.linalg.norm(np.array(dists_to_root_1) - np.array(dists_to_root_2))


## An iterative function to find the minimum kc distance by changing the ordering of each individual (i.e. individual leaves are interchangeable)
#def get_min_kc_dist(tree1, tree2, inds, max_itr = 10):
    
    #new_inds = inds[:]
    
    #for itr in range(max_itr):
        
        #for i in range(len(inds)):
            
            #ind_orderings = list(itertools.permutations(inds[i]))
            
            #dists = []
            
            #for ind_ordering in ind_orderings:
                #current_inds = new_inds[:]
                #current_inds[i] = ind_ordering
                #dists.append(get_kc_dist(tree1, tree2, tree2_leaf_reorder=[i for ind in current_inds for i in ind]))
            
            #new_inds[i] = ind_orderings[np.argmin(dists)]
        
        #if itr > 0 and new_inds == previous_new_inds: break
        #else:
            #previous_new_inds = new_inds[:]
    
    ##finally, return the optimised kc distance
    #return get_kc_dist(tree1, tree2, tree2_leaf_reorder=[i for ind in new_inds for i in ind])


#def get_kc_dist_ts(ts1, ts2, minimize=False, inds=None):
    ##a ts-based KC distance function exists in tskit.
    ## But this function gets one for each marginal tree (along with the weighting)
    ## This is useful for dealing with diploid interchangeable treesequences
    ## where we might want to pick the lowest kc for each tree from a set of possible ones 
    
    ##if minimizing, check that individuals are provided
    #if (minimize and inds is None) or (inds is not None and not minimize):
        #raise ValueError("Set minimize to true and provide individual indices if you want to get minimal kc distances")
    
    #overlaps = []
    #dists = []
    
    #counter1 = 0
    #counter2 = 0
    
    #n_trees1 = ts1.num_trees
    #n_trees2 = ts2.num_trees
    
    #trees1 = ts1.trees(sample_lists=True) #sample_lists=True is apparently needed to be able to apply kc_distance on marginal trees
    #trees2 = ts2.trees(sample_lists=True)
    
    #tree1 = next(trees1)
    #tree2 = next(trees2)
    
    #while True:
        
        ##check for overlap
        #if tree2.interval[0] < tree1.interval[1] and tree2.interval[1] > tree1.interval[0]:
            #overlaps.append(min([tree1.interval[1], tree2.interval[1]]) - max([tree1.interval[0], tree2.interval[0]]))
            #if minimize:
                #dists.append(get_min_kc_dist(tree1, tree2, inds))
            #else:
                #dists.append(tree1.kc_distance(tree2))
        
        ##get next tree from whichever ends sooner
        #if tree1.interval[1] <= tree2.interval[1]:
            #counter1 += 1
            #if counter1 == n_trees1: break
            #tree1 = next(trees1)
        #else:
            #counter2 += 1
            #if counter2 == n_trees2: break
            #tree2 = next(trees2)
    
    #return np.average(dists, weights=overlaps)


#def make_individual_indices(ploidies):
    #partitions = np.cumsum(ploidies)[:-1]
    #return np.split(np.arange(sum(ploidies)), partitions)


#a function that allows calculation of minimum quartet distance between tree sequences (for higher ploidies).
# I'm not using this now because it's too slow.
# The other function I wrote is also faster because it calculates quartet IDs once for each tree, not for each overlap
#def get_quartet_dist_ts_OLD(ts1, ts2, unrooted=False, minimize=False, inds = None, approximation_subset_size=None):
    
    #assert not (minimize and approximation_subset_size is not None), "Cannot find minimal dist using approximate method"
    
    #overlaps = []
    #dists = []
    
    #counter1 = 0
    #counter2 = 0
    
    #n_trees1 = ts1.num_trees
    #n_trees2 = ts2.num_trees
    
    #trees1 = ts1.trees()
    #trees2 = ts2.trees()
    
    #tree1 = next(trees1)
    #tree2 = next(trees2)
    
    #while True:
        
        ##check for overlap
        #if tree2.interval[0] < tree1.interval[1] and tree2.interval[1] > tree1.interval[0]:
            #overlaps.append(min([tree1.interval[1], tree2.interval[1]]) - max([tree1.interval[0], tree2.interval[0]]))
            #if not minimize:
                #dists.append(TopologySummary.get_quartet_dist(tree1, tree2, unrooted=unrooted,
                                                            #approximation_subset_size=approximation_subset_size))
            #else:
                #dists.append(TopologySummary.get_min_quartet_dist(tree1, tree2, inds, unrooted=unrooted))
        
        ##get next tree from whichever ends sooner
        #if tree1.interval[1] <= tree2.interval[1]:
            #counter1 += 1
            #if counter1 == n_trees1: break
            #tree1 = next(trees1)
        #else:
            #counter2 += 1
            #if counter2 == n_trees2: break
            #tree2 = next(trees2)
    
    #return np.average(dists, weights=overlaps)


def get_hamming_dist(arr1, arr2):
    dif = arr1 - arr2
    return np.mean(dif != 0)

def get_quartet_dist_ts(ts1, ts2, unrooted=False, return_n_tree_shapes=False, approximation_subset_size=None):
    
    overlaps = []
    dists = []
    
    counter1 = 0
    counter2 = 0
    
    n_trees1 = ts1.num_trees
    n_trees2 = ts2.num_trees
    
    if return_n_tree_shapes:
        n_changes1 = 0
        n_changes2 = 0
    
    trees1 = ts1.trees()
    trees2 = ts2.trees()
    
    tree1 = next(trees1)
    tree2 = next(trees2)
    
    topoSummary1 = TopologySummary.TopologySummary(tree1)
    topoSummary2 = TopologySummary.TopologySummary(tree2)
    
    leaves = list(topoSummary1.leavesRetained)
    
    if approximation_subset_size is None:
        quartetIDs1 = np.array(topoSummary1.get_all_quartet_IDs(unrooted=unrooted))
        quartetIDs2 = np.array(topoSummary2.get_all_quartet_IDs(unrooted=unrooted))
    else:
        
        if unrooted:
            quartets =  [random.sample(leaves, 4) for i in range(approximation_subset_size)] #the same set of quartets will be used throughout the tree sequence
            quartetIDs1 = np.array([topoSummary1.get_quartet_ID(quartet) for quartet in quartets], dtype=int)
            quartetIDs2 = np.array([topoSummary2.get_quartet_ID(quartet) for quartet in quartets], dtype=int)
        else:
            trios = [random.sample(leaves, 3) for i in range(approximation_subset_size)] # use trios if doing rooted distance
            quartetIDs1 = np.array([topoSummary1.get_quartet_ID(trio + [tree1.root]) for trio in trios], dtype=int)
            quartetIDs2 = np.array([topoSummary2.get_quartet_ID(trio + [tree2.root]) for trio in trios], dtype=int)
    
    while True:
        
        #check for overlap
        if tree2.interval[0] < tree1.interval[1] and tree2.interval[1] > tree1.interval[0]:
            overlaps.append(min([tree1.interval[1], tree2.interval[1]]) - max([tree1.interval[0], tree2.interval[0]]))
            dists.append(get_hamming_dist(quartetIDs1, quartetIDs2))
        
        #get next tree from whichever ends sooner
        if tree1.interval[1] <= tree2.interval[1]:
            counter1 += 1
            if counter1 == n_trees1: break
        
            quartetIDs1_previous = quartetIDs1
            
            tree1 = next(trees1)
            topoSummary1 = TopologySummary.TopologySummary(tree1)
            
            if approximation_subset_size is None:
                quartetIDs1 = np.array(topoSummary1.get_all_quartet_IDs(unrooted=unrooted))
            else:
                if unrooted:
                    quartetIDs1 = np.array([topoSummary1.get_quartet_ID(quartet) for quartet in quartets], dtype=int)
                else:
                    quartetIDs1 = np.array([topoSummary1.get_quartet_ID(trio + [tree1.root]) for trio in trios], dtype=int)
            
            if return_n_tree_shapes:
                if not np.all(quartetIDs1 == quartetIDs1_previous): n_changes1 += 1
        
        else:
            counter2 += 1
            if counter2 == n_trees2: break
            
            quartetIDs2_previous = quartetIDs2
            
            tree2 = next(trees2)
            topoSummary2 = TopologySummary.TopologySummary(tree2)
            
            if approximation_subset_size is None:
                quartetIDs2 = np.array(topoSummary2.get_all_quartet_IDs(unrooted=unrooted))
            else:
                if unrooted:
                    quartetIDs2 = np.array([topoSummary2.get_quartet_ID(quartet) for quartet in quartets], dtype=int)
                else:
                    quartetIDs2 = np.array([topoSummary2.get_quartet_ID(trio + [tree2.root]) for trio in trios], dtype=int)
            
            if return_n_tree_shapes:
                if not np.all(quartetIDs2 == quartetIDs2_previous): n_changes2 += 1
    
    dist = np.average(dists, weights=overlaps)
    
    if return_n_tree_shapes:
        return (dist, n_changes1+1, n_changes2+1,)
    
    return dist


# error rates from Albers and McVeabn 2020 https://doi.org/10.1371/journal.pbio.3000586
# The genotype is 0 (hom ref), 1(het) or 2 (hom alt)
# convert_to is self explanatory
# some rates are divided by 2. This is when there was one rate given in the paper but two possible phases, so we halve them

error_model = {0:{"convert_to":[(0,1), (1,0), (1,1)], "rates":[0.00041/2, 0.00041/2, 0.00017]},
               1:{"convert_to":[(0,0), (1,1)], "rates":[0.00550, 0.00168]},
               2:{"convert_to":[(0,0), (0,1), (1,0)], "rates":[0.00034, 0.00231/2, 0.00231/2]}}


def add_genotyping_errors(haps, error_model, seed=None):
    
    #pick random haplotype to represent the reference
    ref = haps[:,[0]]
    
    #represent genotypes as 0 for ref or 1 for (any) alt
    alt_match = (haps!=ref)*1
    
    #sum the first and second, third and fourth, etc to make diploid genos
    diploid_geno = np.apply_along_axis(lambda x: x[::2] + x[1::2] , 1, alt_match)
    
    #because there are different kinds of errors with different rates,
    #we need to sample the indices to alter without overlaps
    #so we randomly shuffle, then pick the n first values, where n depends on each rate
    #function to give multiple non-overlapping indices when there are multiple things to convert to
    def sample_indices(n_indices, rates, seed=None):
        
        rng = np.random.RandomState(seed)
        
        all_indices = np.arange(n_indices)
        rng.shuffle(all_indices)  # Shuffle once
        
        start_idx = 0
        
        for rate in rates:
            if start_idx < n_indices:
                n_errors = rng.binomial(n_indices, rate)
                error_indices = all_indices[start_idx:start_idx + n_errors]
                start_idx += n_errors
            else:
                error_indices = np.array([])
            
            yield error_indices
    
    #for each genotype, get the indices of each possible convert_to, and convert them
    
    gts = error_model.keys()
    
    #a 2D set of indices over the whole datafram indicate where each genotype state is
    matches_by_genotype = dict([(g, np.where(diploid_geno == g),) for g in gts])
    
    for g in gts:
        matches = matches_by_genotype[g]
        n = matches[0].shape[0]
        indices_to_alter = sample_indices(n, error_model[g]["rates"], seed)
        for i, indices in enumerate(indices_to_alter):
            convert_to = error_model[g]["convert_to"][i]
            for j in indices:
                row = matches[0][j]
                columns = [matches[1][j]*2, matches[1][j]*2 + 1]
                haps[row,columns] = convert_to


def infer_polarity(haps, OG):
    
    assert OG.shape[1] == 1 and OG.shape[0] == haps.shape[0], "OG should be a single column array giving the (consensus) outgroup allele"
    
    matches_og = haps == OG          # Where haps[i,j] == OG[i]
    is_zero = haps == 0             # Where haps[i,j] == 0
    
    # Apply the swaps
    haps[matches_og] = 0
    
    #make a version of OG that is same dimensions of MAT and use that to replace
    haps[is_zero] = np.broadcast_to(OG, haps.shape)[is_zero]


def make_haplotype_array(mts, columns=None, OG_column_for_repolarisation=None):
    
    positions = np.array([variant.position for variant in mts.variants()], dtype=float)
    all_haps = mts.genotype_matrix()
    
    if columns is not None:
        haps = all_haps[:,columns]
    else:
        haps = all_haps
    
    #variants only
    var = ~np.all(haps == haps[:, [0]], axis=1)
    positions = positions[var]
    haps = haps[var]
    
    if OG_column_for_repolarisation is not None:
        OG = all_haps[:,[OG_column_for_repolarisation]][var]
        infer_polarity(haps, OG)
    
    return (positions, haps)


def tree_mean_children(t):
    return np.mean([len(t.children(n)) for n in t.nodes() if len(t.children(n)) != 0])

def mean_children(ts):
    return np.mean([tree_mean_children(t) for t in ts.trees()])



####################################################################################################
################################### once-off sim and comparison ####################################
####################################################################################################


n=100 # MAKE SURE THAT n IS A MULTIPLE OF ALL DESIRED PLOIDLY LEVELS.
N=100000
t=0.001
seq_len=1e5
recombination_rate=1e-8
gene_conversion_rate=None
gene_conversion_tract_length=None
seed=133

#graph = two_pop_model(N=N, t=t)
graph = two_pop_OG_model(N=1e5, tAB=0.001, tABC=2e5)

ts = msprime.sim_ancestry(samples={"A": n/2, "B": n/2, "C":1},
                          demography = msprime.Demography.from_demes(graph),
                          sequence_length = seq_len, discrete_genome=True,
                          recombination_rate = recombination_rate,
                          gene_conversion_rate=gene_conversion_rate,
                          gene_conversion_tract_length=gene_conversion_tract_length,
                          ploidy=1, #we always simulate haploids. We can make any ploidy in the genotype data later
                          random_seed=seed)

mutation_rate = 2e-8

mts = msprime.sim_mutations(ts, rate = mutation_rate, discrete_genome=True, random_seed=seed)

#mts = filter_to_infinite_sites(mts)

positions, haps = make_haplotype_array(mts, columns=range(n))







#positions, haps = make_haplotype_array(mts, columns=range(n), OG_column_for_repolarisation=n)

#1 - np.all(haps==haps_repol, axis=1).sum() / haps.shape[0]

#np.where(haps.max(axis=1) > 1)

#add errors if desired
add_genotyping_errors(haps, error_model, seed)


ts_ingroup = ts.simplify(samples=range(n))
####################################################### tsinfer ############################################################
ts_tsinfer = run_tsinfer(haps, positions, seq_len)
print(ts_tsinfer.num_trees)
print(mean_children(ts_tsinfer))

print(get_quartet_dist_ts(ts_ingroup, ts_tsinfer, return_n_tree_shapes=True, approximation_subset_size=None))

####################################################### #sticcs ############################################################
t = time.time()
ts_sticcs = run_sticcs(haps, positions, seq_len, forced_ploidy=1, second_chances=False, silent=True)
ts_sticcs_I = run_sticcs(hapsI, positionsI, seq_len, forced_ploidy=1, second_chances=False, silent=True)
time.time() - t
print(ts_sticcs.num_trees)
print(mean_children(ts_sticcs))

print(get_quartet_dist_ts(ts_ingroup, ts_sticcs, return_n_tree_shapes=True, approximation_subset_size=100))
print(get_quartet_dist_ts(ts_ingroup, ts_sticcs_I, return_n_tree_shapes=True, approximation_subset_size=100))


#ts_sticcs_2 = run_sticcs(haps, positions, seq_len, forced_ploidy=2, second_chances=False, silent=True)

#ts_sticcs_4 = run_sticcs(haps, positions, seq_len, forced_ploidy=4, second_chances=False, silent=True)

####################################################### #relate ############################################################

relate_bin="/home/smarti11/software/relate_v1.2.2_x86_64_static/bin/"

ts_relate = run_relate(haps, positions, N=N, mu=mutation_rate, rec=recombination_rate, relate_bin=relate_bin)

print(get_quartet_dist_ts(ts_ingroup, ts_relate, return_n_tree_shapes=True, approximation_subset_size=None))

###################################################### #argweaver #########################################################

t = time.time()
tl_list_argweaver = run_argweaver(haps, positions, N=N, rec=recombination_rate, mu=mutation_rate, sequence_start=0, sequence_end=seq_len)
time.time() - t

print(np.array([get_quartet_dist_ts(ts_ingroup, tl, return_n_tree_shapes=True, approximation_subset_size=None) for tl in tl_list_argweaver]).mean(axis=0))



####################################################### #singer ############################################################

singer_path = "/home/smarti11/software/singer-0.1.8-beta-linux-x86_64/"


t = time.time()
ts_list_singer = run_singer(haps, positions, N=N, rec=recombination_rate, mu=mutation_rate, sequence_start=0, sequence_end=seq_len, singer_path=singer_path, burnin=500, sample_n_args=10, spacing=20)
time.time() - t


print(np.array([get_quartet_dist_ts(ts_ingroup, ts, return_n_tree_shapes=True, approximation_subset_size=None) for ts in ts_list_singer]).mean(axis=0))




#############################################################################################################
###########################  multiple comparison quartet dists  #############################################
#############################################################################################################


def infer_and_get_quartet_dists(ts, n, N, rec, mu,
                                infinite_sites=False,
                                infer_pol=False, gt_error=False, seed=None,
                                approximation_subset_size=None):
    result = []
    
    #add mutations
    mts = msprime.sim_mutations(ts, rate = mutation_rate, random_seed=seed)
    
    #make infinite sites if required
    if infinite_sites:
        mts = filter_to_infinite_sites(mts)
    
    if infer_pol:
        positions, haps = make_haplotype_array(mts, columns=range(n), OG_column_for_repolarisation=n)
    else:
        positions, haps = make_haplotype_array(mts, columns=range(n))
    
    if len(positions) < 10: #if there are too few SNPs
        return [np.nan]*25
    
    #remove any SNP at position zero since argweaver can't handle that
    if positions[0] == 0:
        positions = positions[1:]
        haps = haps[1:,:]
    
    #add errors if desired
    if gt_error:
        add_genotyping_errors(haps, error_model, seed)
    
    seq_len = mts.sequence_length
    
    ts_ingroup = ts.simplify(samples=range(n))
    
    #decide whether approximation is required
    n_quartets = n_choose_k(haps.shape[1], 3) #three because these are rooted quartets
    subset_size = approximation_subset_size if n_quartets > approximation_subset_size else None
    
    #sticcs
    t = time.time()
    ts_sticcs = run_sticcs(haps, positions, seq_len, forced_ploidy=1, second_chances=False, silent=True)
    elapsed = time.time() - t
    dist_data = get_quartet_dist_ts(ts_ingroup, ts_sticcs, return_n_tree_shapes=True, approximation_subset_size=None)
    result.append(dist_data[1]) # number of changes in topology in the simulated ts (record this only once)
    result.append(dist_data[0]) #quartet distance
    result.append(dist_data[2]) #number of changes in topology in infered ts
    result.append(mean_children(ts_sticcs)) # mean number of children per node
    result.append(elapsed)
    
    #sticcs_sc
    t = time.time()
    ts_sticcs_sc = run_sticcs(haps, positions, seq_len, forced_ploidy=1, second_chances=True, silent=True)
    elapsed = time.time() - t
    dist_data = get_quartet_dist_ts(ts_ingroup, ts_sticcs_sc, return_n_tree_shapes=True, approximation_subset_size=None)
    result.append(dist_data[0]) #quartet distance
    result.append(dist_data[2]) #number of changes in topology in infered ts
    result.append(mean_children(ts_sticcs_sc)) # mean number of children per node
    result.append(elapsed)
    
    #tsinfer
    t = time.time()
    ts_tsinfer = run_tsinfer(haps, positions, seq_len)
    elapsed = time.time() - t
    dist_data = get_quartet_dist_ts(ts_ingroup, ts_tsinfer, return_n_tree_shapes=True, approximation_subset_size=None)
    result.append(dist_data[0]) #quartet distance
    result.append(dist_data[2]) #number of changes in topology in infered ts
    result.append(mean_children(ts_tsinfer)) # mean number of children per node
    result.append(elapsed)
    
    #ts_relate
    relate_bin="/home/smarti11/software/relate_v1.2.2_x86_64_static/bin/"
    t = time.time()
    ts_relate = run_relate(haps, positions, N=N, mu=mu, rec=rec, relate_bin=relate_bin)
    elapsed = time.time() - t
    dist_data = get_quartet_dist_ts(ts_ingroup, ts_relate, return_n_tree_shapes=True, approximation_subset_size=None)
    result.append(dist_data[0]) #quartet distance
    result.append(dist_data[2]) #number of changes in topology in infered ts
    result.append(mean_children(ts_relate)) # mean number of children per node
    result.append(elapsed)
    
    #singer
    singer_path = "/home/smarti11/software/singer-0.1.8-beta-linux-x86_64/"
    t = time.time()
    ts_list_singer = run_singer(haps, positions, N=N, rec=rec, mu=mu, sequence_start=0, sequence_end=seq_len, singer_path=singer_path, burnin=500, sample_n_args=10, spacing=20)
    elapsed = time.time() - t
    dist_data = np.array([get_quartet_dist_ts(ts_ingroup, ts, return_n_tree_shapes=True, approximation_subset_size=None) for ts in ts_list_singer]).mean(axis=0)
    result.append(dist_data[0]) #quartet distance
    result.append(dist_data[2]) #number of changes in topology in infered ts
    result.append(np.mean([mean_children(ts) for ts in ts_list_singer]))
    result.append(elapsed)
    
    #argweaver
    t = time.time()
    tl_list_argweaver = run_argweaver(haps, positions, N=N, rec=rec, mu=mu, sequence_start=0, sequence_end=seq_len, burnin=500, sample_n_args=10, spacing=10)
    elapsed = time.time() - t
    dist_data = np.array([get_quartet_dist_ts(ts_ingroup, tl, return_n_tree_shapes=True, approximation_subset_size=None) for tl in tl_list_argweaver]).mean(axis=0)
    result.append(dist_data[0]) #quartet distance
    result.append(dist_data[2]) #number of changes in topology in infered ts
    result.append(np.mean([mean_children(tl) for tl in tl_list_argweaver]))
    result.append(elapsed)
    
    return result


#############################################################################################################


threads = 20

nreps=20

inf_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5], "mutation_rate": [1e-8],
               "infinite_sites":[True],
          #"n":[4, 8],
          "n":[4, 8, 16, 32, 64],
          #"n":[64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None],
            "gt_error":[False], "infer_pol": [False],
            "t_split": [0.001],}

fin_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5], "mutation_rate": [1e-8],
               "infinite_sites":[True],
          #"n":[4, 8],
          "n":[4, 8, 16, 32, 64],
          #"n":[64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None],
            "gt_error":[False], "infer_pol": [False],
            "t_split": [0.001],}

error_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5], "mutation_rate": [1e-8],
               "infinite_sites":[False],
          #"n":[4, 8],
          "n":[4, 8, 16, 32, 64],
          #"n":[64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None],
            "gt_error":[True], "infer_pol": [False],
            "t_split": [0.001],}

pol_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5], "mutation_rate": [1e-8],
               "infinite_sites":[False],
          #"n":[4, 8],
          "n":[4, 8, 16, 32, 64],
          #"n":[64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None],
            "gt_error":[False], "infer_pol": [True],
            "t_split": [0.001],}

mu_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5], "mutation_rate": [1e-9],
               "infinite_sites":[False],
          #"n":[4, 8],
          "n":[4, 8, 16, 32, 64],
          #"n":[64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None],
            "gt_error":[False], "infer_pol": [False],
            "t_split": [0.001],}

gconv_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5], "mutation_rate": [1e-8],
               "infinite_sites":[False],
          #"n":[4, 8],
          "n":[4, 8, 16, 32, 64],
          #"n":[64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [1e-8],
            "gene_conversion_tract_length": [300],
            "gt_error":[False], "infer_pol": [False],
            "t_split": [0.001],}

struc_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5], "mutation_rate": [1e-8],
               "infinite_sites":[False],
          #"n":[4, 8],
          "n":[4, 8, 16, 32, 64],
          #"n":[64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None],
            "gt_error":[False], "infer_pol": [False],
            "t_split": [2e5],} #corresponds to fst ~ 0, 0.5



combos = (list(itertools.product(*inf_params.values())) +
          list(itertools.product(*fin_params.values())) +
          list(itertools.product(*error_params.values())) +
          list(itertools.product(*pol_params.values())) +
          list(itertools.product(*mu_params.values())) +
          list(itertools.product(*gconv_params.values())) +
          list(itertools.product(*struc_params.values())))

n_combos = len(combos)

results = []

for combo_number,combo in enumerate(combos):
    
    print(f"Parameter combination {combo_number + 1} out of {n_combos}...", file=sys.stderr)
    
    seed, nreps, N, seq_len, mutation_rate, infinite_sites, n, recombination_rate, gconvrate, gconvlen, gt_error, infer_pol, t_split = combo
    
    graph = two_pop_OG_model(N=N, tAB=t_split, tABC=t_split+2e5)
    
    print(f"Simulating...", file=sys.stderr)
    
    ts_reps = list(msprime.sim_ancestry(samples={"A": n/2, "B": n/2, "C":1},
                   demography = msprime.Demography.from_demes(graph),
                   sequence_length = seq_len,
                   recombination_rate = recombination_rate, ploidy=1, #we always simulate haploids. We can make any ploidy in the genotype data later
                   gene_conversion_rate=gconvrate,
                   gene_conversion_tract_length=gconvlen,
                   random_seed=seed, num_replicates=nreps))
    
    # Make a temporary function to do the inferring and getting distances for each
    # This is because pool.map requires a function with only 1 argument
    print(f"Inferring tree sequences and computing distances...", file=sys.stderr)
    
    def _infer_and_get_dists_(i):
        return infer_and_get_quartet_dists(ts=ts_reps[i], n=n, N=N, rec=recombination_rate, mu=mutation_rate, infinite_sites=infinite_sites, infer_pol=infer_pol, gt_error=gt_error, seed=seed, approximation_subset_size=100)
    
    #for debugiing -> run in series
    #for i in range(nreps):
        #print(i)
        #_infer_and_get_dists_(i)
    
    with Pool(threads) as p:
        result_reps = np.array(p.map(_infer_and_get_dists_, range(nreps)))
    
    results.append(result_reps.mean(axis=0))



params_path = "/mnt/loki/martin/simon/sticcs_twisst2/sim/params.csv"
header = list(fin_params.keys())
np.savetxt(params_path, np.array(combos), fmt='%s', delimiter = ",", header=",".join(header), comments="")

inference_path = "/mnt/loki/martin/simon/sticcs_twisst2/sim/arg_inference_data.csv"
header = ["sim_trees"]
for method in ["sticcs", "sticcs_sc", "tsinfer", "relate", "singer", "argweaver"]:
    header += [method + "_dist", method + "_trees", method + "_mean_children", method + "_time"]

np.savetxt(inference_path, np.array(results), fmt='%s', delimiter = ",", header=",".join(header), comments="")



#############################################################################################################
###########################  multiple comparison speed by seq length #######################################
#############################################################################################################

def infer_and_get_time(ts, n, N, rec, mu, seed=None):
    
    #add mutations
    mts = msprime.sim_mutations(ts, rate = mutation_rate, random_seed=seed)
    
    positions, haps = make_haplotype_array(mts, columns=range(n))
    
    #remove any SNP at position zero since argweaver can't handle that
    if positions[0] == 0:
        positions = positions[1:]
        haps = haps[1:,:]
    
    result = [mts.num_trees, len(positions)]
    
    if len(positions) < 10: #if there are too few SNPs
        result += [np.nan]*6
        return result
    
    seq_len = mts.sequence_length
    
    #sticcs
    t = time.time()
    ts_sticcs = run_sticcs(haps, positions, seq_len, forced_ploidy=1, second_chances=False, silent=True)
    elapsed = time.time() - t
    result.append(elapsed)
    
    #sticcs_sc
    t = time.time()
    ts_sticcs_sc = run_sticcs(haps, positions, seq_len, forced_ploidy=1, second_chances=True, silent=True)
    elapsed = time.time() - t
    result.append(elapsed)
    
    #tsinfer
    t = time.time()
    ts_tsinfer = run_tsinfer(haps, positions, seq_len)
    elapsed = time.time() - t
    result.append(elapsed)
    
    #ts_relate
    relate_bin="/home/smarti11/software/relate_v1.2.2_x86_64_static/bin/"
    t = time.time()
    ts_relate = run_relate(haps, positions, N=N, mu=mu, rec=rec, relate_bin=relate_bin)
    elapsed = time.time() - t
    result.append(elapsed)
    
    #singer
    singer_path = "/home/smarti11/software/singer-0.1.8-beta-linux-x86_64/"
    t = time.time()
    ts_list_singer = run_singer(haps, positions, N=N, rec=rec, mu=mu, sequence_start=0, sequence_end=seq_len, singer_path=singer_path, burnin=500, sample_n_args=10, spacing=20)
    elapsed = time.time() - t
    result.append(elapsed)
    
    #argweaver
    if seq_len <= 1e5:
        t = time.time()
        tl_list_argweaver = run_argweaver(haps, positions, N=N, rec=rec, mu=mu, sequence_start=0, sequence_end=seq_len, burnin=500, sample_n_args=10, spacing=10)
        elapsed = time.time() - t
        result.append(elapsed)
    else:
        result.append(np.nan)
    
    return result


threads = 20

nreps=20

equiv_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5, 1e6, 1e7], "mutation_rate": [1e-8],
               "infinite_sites":[True],
            "n":[4, 8, 16, 32, 64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None]}

lowMu_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5, 1e6, 1e7], "mutation_rate": [1e-9],
               "infinite_sites":[True],
            "n":[4, 8, 16, 32, 64],
            "recombination_rate": [1e-8],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None]}

lowRho_params = {"seed":[123], "nreps":[nreps], "N": [100000], "seq_len": [1e5, 1e6, 1e7], "mutation_rate": [1e-8],
               "infinite_sites":[True],
            "n":[4, 8, 16, 32, 64],
            "recombination_rate": [1e-9],
            "gene_conversion_rate": [None],
            "gene_conversion_tract_length": [None]}


combos = (list(itertools.product(*equiv_params.values())) +
          list(itertools.product(*lowMu_params.values())) +
          list(itertools.product(*lowRho_params.values())))

n_combos = len(combos)

results = []

for combo_number,combo in enumerate(combos):
    
    print(f"Parameter combination {combo_number + 1} out of {n_combos}...", file=sys.stderr)
    
    seed, nreps, N, seq_len, mutation_rate, infinite_sites, n, recombination_rate, gconvrate, gconvlen = combo
    
    print(f"Simulating...", file=sys.stderr)
    
    ts_reps = list(msprime.sim_ancestry(samples=n,
                   population_size=N,
                   sequence_length = seq_len,
                   recombination_rate = recombination_rate, ploidy=1, #we always simulate haploids. We can make any ploidy in the genotype data later
                   gene_conversion_rate=gconvrate,
                   gene_conversion_tract_length=gconvlen,
                   random_seed=seed, num_replicates=nreps))
    
    # Make a temporary function to do the inferring and getting distances for each
    # This is because pool.map requires a function with only 1 argument
    print(f"Inferring tree sequences and computing distances...", file=sys.stderr)
    
    def _infer_and_get_time_(i):
        return infer_and_get_time(ts=ts_reps[i], n=n, N=N, rec=recombination_rate, mu=mutation_rate, seed=seed)
    
    #for debugiing -> run in series
    #for i in range(nreps):
        #print(i)
        #_infer_and_get_time_(i)
    
    with Pool(threads) as p:
        result_reps = np.array(p.map(_infer_and_get_time_, range(nreps)))
    
    results.append(result_reps.mean(axis=0))




params_path = "/mnt/loki/martin/simon/sticcs_twisst2/sim/timetest_params.csv"
header = list(equiv_params.keys())
np.savetxt(params_path, np.array(combos), fmt='%s', delimiter = ",", header=",".join(header), comments="")

time_path = "/mnt/loki/martin/simon/sticcs_twisst2/sim/times_data.csv"
header = ["trees", "variants", "sticcs", "sticcs_sc", "tsinfer", "relate", "singer", "argweaver"]
np.savetxt(time_path, np.array(results), fmt='%s', delimiter = ",", header=",".join(header), comments="")

