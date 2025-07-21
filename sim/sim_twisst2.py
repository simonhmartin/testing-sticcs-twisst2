import time, random, gzip, sys, os

import msprime
import tskit
import tsinfer
import demes

import numpy as np

from sticcs import sticcs
from twisst2 import twisst2
from twisst2 import parse_newick
from twisst2 import bionj

from multiprocessing import Pool
import subprocess
import tempfile
from collections import defaultdict


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


def run_argweaver(haps, positions, N, rec, mu, sequence_start, sequence_end, burnin=500, sample_n_args=10, spacing=10, unphased=False):
    n_sites, n_haps = haps.shape
    iters = burnin + sample_n_args*spacing
    
    #tmpdir = tempfile.mkdtemp(dir=".") #for debugging
    
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        
        sites_file=f"{tmpdir}/temp.sites"
        region_name = "chr1"
        
        # Write .sites file in ARGweaver format
        with open(sites_file, "w") as f:
            # Header line with sample names
            if unphased:
                assert n_haps % 2 == 0, "There must be an even number of haplotypes"
                hap_names = [f"ind{i}_{j}" for i in range(n_haps//2) for j in [1,2]]
            else:
                hap_names = [f"ind{i}" for i in range(n_haps)]
            
            hap_number_dict = dict([(hap_names[i], i) for i in range(n_haps)])
            
            header = "NAMES\t" + "\t".join(hap_names)
            f.write(header + "\n")
            
            # Region line
            region_line = f"REGION\t{region_name}\t{sequence_start}\t{sequence_end}"
            f.write(region_line + "\n")
            
            # Site data - one line per segregating site
            for site_idx, pos in enumerate(positions):
                
                # Get alleles for this site
                alleles = []
                for hap_idx in range(n_haps):
                    allele_idx = haps[site_idx, hap_idx]
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
                #annoyingly, the trees have numeric leaves, but these are NOT in the order of the haps
                #The order if given by th efirst line, so we can figure out which hap each leaf points to
                names = argfile.readline().split()[1:]
                leaf_idx_dict = dict([(str(i), hap_number_dict[names[i]]) for i in range(n_haps)])
                
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


def bionj_windows(haps, positions, window_size, seq_len):
    nsites = haps.shape[0]
    assert nsites == len(positions)
    for window_start in range(0, nsites, window_size):
        window_end = window_start + window_size
        if window_end >= nsites: window_end = nsites-1
        current_hap_rows = haps[window_start:window_end,:].transpose()
        distance_matrix = bionj.distMatrix(current_hap_rows)
        tree = bionj.bionj(distance_matrix)
        #intervals are defined by gaps between SNPs
        startpos = 1 if window_start==0 else int((positions[window_start-1]+1 + positions[window_start]+1)/2)
        endpos = seq_len if window_end == nsites-1 else int((positions[window_end]+1 + positions[window_end+1]+1)/2)
        tree.interval = [startpos, endpos]
        yield tree


def admix_model(t_C=1e4, t_ABC=5e4, t_ABCD=2e5,                     #split times
                    N_A=1e5, N_B=1e5, N_Cstart = 1e4, N_Cend = 1e4, N_D=1e5, N_ABC=1e5, N_ABCD=1e5):   #population sizes
    
    #use the deme builder to set up the demographic history
    b = demes.Builder(time_units="generations")
    b.add_deme("ABCD",                      epochs=[{"start_size":N_ABCD, "end_time":t_ABCD}])
    b.add_deme("ABC", ancestors=["ABCD"], epochs=[{"start_size":N_ABC,  "end_time":t_ABC}])
    b.add_deme("D",     ancestors=["ABCD"],  epochs=[{"start_size":N_D}])
    b.add_deme("A",     ancestors=["ABC"],  epochs=[{"start_size":N_A}])
    b.add_deme("B",     ancestors=["ABC"],  epochs=[{"start_size":N_B}])
    b.add_deme("C",     ancestors=["A", "B"], proportions=[0.999, 0.001], start_time=t_C, epochs=[{"start_size":N_Cstart, "end_size":N_Cend}])
    
    graph = b.resolve()
    
    return(graph)


def fourPop_model(t_AB=1e4, t_ABC=5e4, t_ABCD=2e5, N=1e5):
    
    #use the deme builder to set up the demographic history
    b = demes.Builder(time_units="generations")
    b.add_deme("ABCD",                      epochs=[{"start_size":N, "end_time":t_ABCD}])
    b.add_deme("ABC", ancestors=["ABCD"], epochs=[{"start_size":N,  "end_time":t_ABC}])
    b.add_deme("D",     ancestors=["ABCD"],  epochs=[{"start_size":N}])
    b.add_deme("AB",     ancestors=["ABC"],  epochs=[{"start_size":N,  "end_time":t_AB}])
    b.add_deme("C",     ancestors=["ABC"],  epochs=[{"start_size":N}])
    b.add_deme("A",     ancestors=["AB"], epochs=[{"start_size":N}])
    b.add_deme("B",     ancestors=["AB"], epochs=[{"start_size":N}])
    
    graph = b.resolve()
    
    return(graph)


def admix3(t_C=1e4, t_ABC=1e5, t_ABCD=5e5,                     #split times
                    N_A=1e5, N_C=1e5, N_Bstart = 1e4, N_Bend = 1e4, N_ABC=1e5):   #population sizes
    
    #use the deme builder to set up the demographic history
    b = demes.Builder(time_units="generations")
    b.add_deme("ABC", epochs=[{"start_size":N_ABC,  "end_time":t_ABC}])
    b.add_deme("A",     ancestors=["ABC"],  epochs=[{"start_size":N_A}])
    b.add_deme("C",     ancestors=["ABC"],  epochs=[{"start_size":N_C}])
    b.add_deme("B",     ancestors=["A", "C"], proportions=[0.7, 0.3], start_time=t_C, epochs=[{"start_size":N_Bstart, "end_size":N_Bend}])
    
    graph = b.resolve()
    
    return(graph)


def add_missing_data(mat, missing_rate, mean_run):
    
    target = mat.shape[0] * (missing_rate + missing_rate**2)
    
    for ind in range(mat.shape[1]):
        missing = 0 
        while missing < target:
            start = np.random.randint(mat.shape[0])
            run = int(np.ceil(np.random.exponential(mean_run)))
            try:
                mat[start:(start+run), ind] = -1
                missing += run
            except: pass


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


def ds_euc_dist(valuesA, valuesB, max_values, min_values, A_is_target=False, independent=True):
    #calculate maximum distance for each variable
    if not A_is_target: md = (max_values - min_values)**2 
    #if A is target, use maximum distance from A for md
    else: md = np.absolute(np.column_stack([valuesA-max_values,valuesA-min_values])).max(axis=1)**2 
    
    if independent: return np.sqrt( (((valuesA - valuesB)**2)/md).sum() / valuesA.shape[0])
    else: return np.sqrt( (((valuesA - valuesB)**2)/md).sum() / (valuesA.shape[0]-1))


def weighted_nanmean(a, w):
    good = np.where(np.invert(np.isnan(a)))[0]
    return np.average(a[good], weights=w[good])


def weights_error_rate(topocounts_true, topocounts_query, window_size = None):
    
    if not window_size:
        new_intervals = topocounts_true.intervals
        topocounts_true_norm = topocounts_true
        topocounts_query_norm = topocounts_query.recast(new_intervals=new_intervals)
    else:
        topocounts_true_norm = topocounts_true.recast(window_size)
        topocounts_query_norm = topocounts_query.recast(window_size)
    
    n = topocounts_true_norm.intervals.shape[0]
    
    k = topocounts_true_norm.counts.shape[1]
    
    lengths = np.diff(topocounts_true_norm.intervals, axis=1)[:,0] + 1
    
    weights_true = topocounts_true_norm.counts / topocounts_true_norm.totals.reshape((n,1))
    
    weights_query = topocounts_query_norm.counts / topocounts_query_norm.totals.reshape((n,1))
    
    dists = np.array([ds_euc_dist(weights_true[i,],
                                  weights_query[i,],
                                  np.zeros(k), np.array([1]*k),
                                  A_is_target=True, independent=False) for i in range(n)])
    
    return weighted_nanmean(dists, lengths)



def get_pattern_counts(tree):
    n = tree.num_samples()
    pattern_counts = defaultdict(int)
    pattern_counts_dip = defaultdict(int)
    for node in tree.nodes():
        pattern = np.zeros(n, dtype=int)
        for i in tree.leaves(node):
            pattern[i] = 1 
        
        if 1 < sum(pattern) < n:
            pattern_dip = pattern[::2] + pattern[1::2]
            
            pattern_counts[tuple(pattern)] += 1
            pattern_counts_dip[tuple(pattern_dip)] += 1
    
    return (pattern_counts, pattern_counts_dip,)


def analyse_unique_nodes(ts):
    totals_hap = []
    totals_dip = []
    
    for t in ts.trees():
        counts_hap, counts_dip = get_pattern_counts(t)
        totals_hap.append(np.bincount(list(counts_hap.values()), minlength=3))
        totals_dip.append(np.bincount(list(counts_dip.values()), minlength=3))
    
    totals_hap = np.array(totals_hap)[:,1:]
    totals_dip = np.array(totals_dip)[:,1:]
    
    print(f"Trees: {totals_dip.shape[0]}")
    print(f"Nodes: {totals_dip.sum()}")
    print(f"Unique nodes: {totals_dip[:,0].sum()}, ({totals_dip[:,0].sum() / totals_dip.sum()})")
    print(f"Completely unique trees: {np.sum(totals_dip[:,1] == 0)}, ({np.mean(totals_dip[:,1] == 0)})")


####################################################################################################
################################ 4 taxa, 100kb recomnining sequence ################################
####################################################################################################

# we will always simulate haploidas, and later coerce the genotypes into individuals of ploidy 1, 2, or 4.
# This ensures that sim params and Ne stay the same in all tests.
# MAKE SURE THAT n IS A MULTIPLE OF ALL DESIRED PLOIDLY LEVELS.

n=20
sim_ploidy = 1 # NOTE - keep this as 1. We later set ploidy manually by joining individuals.
seq_len=2e5
recombination_rate=1e-8

seed=123
graph = fourPop_model()
run_name = "fourPop_n4_123_l2e5_r1e8"
n_samples={"A": n, "B": n, "C":n, "D":n}


seed=123
graph = admix_model()
run_name = f"admix_n{n}_seed{seed}_l{seq_len}_r{recombination_rate}"
n_samples={"A": n, "B": n, "C":n, "D":n}


seed=6
graph = admix3()
run_name = f"admix3_n{n}_seed{seed}_l{seq_len}_r{recombination_rate}"
n_samples={"A": n, "B": n, "C":n}


total_samples = sum(n_samples.values())
total_subtrees = np.prod(list(n_samples.values()))

ts = msprime.sim_ancestry(samples=n_samples,
                          demography = msprime.Demography.from_demes(graph), sequence_length = seq_len,
                          recombination_rate = recombination_rate, ploidy=sim_ploidy, random_seed=seed)


###################### weights from ts #############################

unrooted = False

n_groups = len(n_samples)

groups_haploid = twisst2.make_numeric_groups([n]*n_groups)


#tskit
#t = time.time()
#topocounts_tskit = twisst2.get_topocounts_tskit(ts, groups = groups_haploid)
#time.time() - t

max_subtrees = n**n_groups

t = time.time()
topocounts = twisst2.get_topocounts(ts.trees(), leaf_groups = groups_haploid, max_subtrees=max_subtrees, unrooted=unrooted)
time.time() - t


#correct intervals to 1 based inclusive
topocounts.intervals[:,0] += 1


outdir = "/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/"
outdir = "./"

with open(outdir + run_name + ".true.topocounts.tsv", "wt") as tcFile:
    topocounts.write(tcFile)

with open(outdir + run_name + ".true.intervals.tsv", "wt") as intsFile:
    topocounts.write_intervals(intsFile)


############################### make mutation marix ####################

mutation_rate = 1e-8
mut_name = f"mu{mutation_rate}"

mts = msprime.sim_mutations(ts, rate = mutation_rate, discrete_genome=True, random_seed=seed)



positions, haps = make_haplotype_array(mts, columns=range(total_samples))
#positions, haps = make_haplotype_array(mts, columns=range(total_samples), OG_column_for_repolarisation=n)

 #how many were incorrectly flipped?
#1 - np.all(haps==haps_repol, axis=1).sum() / haps.shape[0]


#add errors if desired
#add_genotyping_errors(haps, error_model, seed)
#mut_name += "_err"


### optionally add missing data

#if missing_rate > 0:
    #add_missing_data(mat, missing_rate=missing_rate, mean_run=2)
    #mut_name += f"_mis{missing_rate}"


############################## relate ##############################

relate_bin="/home/smarti11/software/relate_v1.2.2_x86_64_static/bin/"
t = time.time()
ts_relate = run_relate(haps, positions, N=1e5, mu=mutation_rate, rec=recombination_rate, relate_bin=relate_bin)
topocounts_relate = twisst2.get_topocounts(ts_relate.trees(), leaf_groups=groups_haploid, max_subtrees=total_subtrees, unrooted=False)
time.time() - t

# relate tree sequence intervals need to be corrected to be 1-based inclusive
topocounts_relate.intervals[:,0] += 1

#round to make all integers
topocounts_relate.intervals = np.floor(topocounts_relate.intervals)

#and extend last interval to the sequence length
topocounts_relate.intervals[-1,1] = seq_len


with open(f"{run_name}_{mut_name}.relate.topocounts.tsv", "wt") as tcFile:
    topocounts_relate.write(tcFile)

with open(f"{run_name}_{mut_name}.relate.intervals.tsv", "wt") as intsFile:
    topocounts_relate.write_intervals(intsFile)


############################## tsinfer ##############################

t = time.time()
ts_tsinfer = run_tsinfer(haps, positions, seq_len)

topocounts_tsinfer = twisst2.get_topocounts(ts_tsinfer.trees(), leaf_groups=groups_haploid, max_subtrees=max_subtrees, unrooted=unrooted)
time.time() - t

#correct intervals to 1 based inclusive
topocounts_tsinfer.intervals[:,0] += 1



with open(f"{run_name}_{mut_name}.tsinfer.topocounts.tsv", "wt") as tcFile:
    topocounts_tsinfer.write(tcFile)

with open(f"{run_name}_{mut_name}.tsinfer.intervals.tsv", "wt") as intsFile:
    topocounts_tsinfer.write_intervals(intsFile)


analyse_unique_nodes(ts_tsinfer)

################################# print some trees to look at #####################################################

start, end = 0, 2000  # Adjust as needed
sub_ts = ts_tsinfer.keep_intervals([[start, end]])
with open("/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/" + run_name + mut_name + mod_name + ".tsinfer.subset.svg", "w") as f:
    f.write(sub_ts.draw_svg(node_labels={}, mutation_labels={}, omit_sites=True, symbol_size=0, x_axis=False, size=(2000,500)))


#single tree
for i in [1,2,ts_tsinfer.num_trees - 1]:
    with open(f"/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/{run_name}{mut_name}{mod_name}.tsinfer.tree{i}.svg", "w") as f:
        f.write(ts_tsinfer.at_index(i).draw_svg(node_labels={}, mutation_labels={}, omit_sites=True, symbol_size=0))


############################ sticcs #################################

# In all cases we simulate diploids. These can be represented as haploids or tetraploids by manipulating the data
# This is preferred over simulating different ploidies in msprime, because that changes the effective population size
# The way we do it ensures that the same exact simulation can be analysed at different ploidy levels without changing any population parameters
# Make sure that n is always a multiple of the ploidy you want to impose on it.

topocounts_twisst2_list = []

#topocounts_twisst2_ave5k_list = []
max_subtrees = 512

for forced_ploidy in [1,2,4]:
    print(forced_ploidy)
    #check that you have a multiple of forced ploidy haplotypes to work with!
    assert n % forced_ploidy == 0
    
    t = time.time()
    
    ploidies = np.array([forced_ploidy]*int(n/forced_ploidy)*n_groups)
    
    der_counts, der_positions = sticcs.get_dac_from_haploid_matrix(haps, positions+1, ploidies)
    
    groups = twisst2.make_numeric_groups([int(n/forced_ploidy)]*n_groups)
    
    sticcs_name = f"sticcs_ploidy{forced_ploidy}_subtre{max_subtrees}"
    
    topocounts_twisst2_list.append(twisst2.get_topocounts_stacking_sticcs(der_counts, der_positions, ploidies, groups,
                                                                chrom_start=1, chrom_len = seq_len, max_subtrees=max_subtrees, second_chances=True, unrooted=unrooted, silent=False))
    
    print(time.time() - t)
    
    with open(f"{run_name}_{mut_name}.{sticcs_name}.topocounts.tsv", "wt") as outfile:
        topocounts_twisst2_list[-1].write(outfile)
    
    with open(f"{run_name}_{mut_name}.{sticcs_name}.intervals.tsv", "wt") as outfile:
        topocounts_twisst2_list[-1].write_intervals(outfile)

window_sizes = [None] + [int(x) for x in [1e3, 2e3, 5e3, 10e3, 20e3, 50e3]]

[weights_error_rate(topocounts, topocounts_twisst2_list[0], window_size=x) for x in window_sizes]
[weights_error_rate(topocounts, topocounts_twisst2_list[1], window_size=x) for x in window_sizes]
[weights_error_rate(topocounts, topocounts_twisst2_list[2], window_size=x) for x in window_sizes]


#or just single ts

ts_sticcs_2 = run_sticcs(haps, positions, seq_len, forced_ploidy=2, second_chances=True, silent=False)

analyse_unique_nodes(ts_sticcs_2)

###################################### stacking sticcs example with single iterations

forced_ploidy = 2
assert n % forced_ploidy == 0

ploidies = np.array([forced_ploidy]*int(n/forced_ploidy)*n_groups)

der_counts, der_positions = sticcs.get_dac_from_haploid_matrix(haps, positions+1, ploidies)

groups = twisst2.make_numeric_groups([int(n/forced_ploidy)]*n_groups)

sticcs_name = f".sticcs_ploidy{forced_ploidy}_subtre8"


i = 2

combo = [g[i-1] for g in groups]

der_counts_sub = der_counts[:, combo]

ploidies_sub = ploidies[(combo,)]

no_missing = np.all(der_counts_sub >= 0, axis=1)

site_sum = der_counts_sub.sum(axis=1)
variable = (1 < site_sum) & (site_sum < sum(ploidies_sub))

usable_sites = np.where(no_missing & variable)

der_counts_sub = der_counts_sub[usable_sites]

positions_sub = positions[usable_sites]

patterns, matches, n_matches  = sticcs.get_patterns_and_matches(der_counts_sub)

clusters = sticcs.get_clusters(patterns, matches, positions_sub, ploidies=ploidies_sub, second_chances=True,
                                seq_start=1, seq_len=seq_len, silent=True)

ts_sticcs = sticcs.infer_ts(patterns, ploidies_sub, clusters, multi_pass = True, silent=True)



for j in [1,2,ts_sticcs.num_trees - 1]:
    with open(f"/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/{run_name}{mut_name}{mod_name}{sticcs_name}.tree.{j}.svg", "w") as f:
        f.write(ts_sticcs.at_index(j).draw_svg(node_labels={}, mutation_labels={}, omit_sites=True, symbol_size=0))



##################### singer ##################
singer_path = "/home/smarti11/software/singer-0.1.8-beta-linux-x86_64/"
t = time.time()
ts_list_singer = run_singer(haps, positions, N=100000, rec=recombination_rate, mu=mutation_rate, sequence_start=0, sequence_end=seq_len, singer_path=singer_path, burnin=500, sample_n_args=10, spacing=20)
#Need to average over all inferred args
topocounts_list_singer = [twisst2.get_topocounts(ts.trees(), leaf_groups=groups_haploid, max_subtrees=total_subtrees, unrooted=False) for ts in ts_list_singer]

for tc in topocounts_list_singer:
    tc.intervals[:,0] += 1

topocounts_list_singer_unrooted = [tc.unroot() for tc in topocounts_list_singer]

[np.mean([weights_error_rate(topocounts, topocounts_singer, window_size=x) for topocounts_singer in topocounts_list_singer]) for x in window_sizes]

#stack the topocounts for export (this is effectively the same as averaging them)

topocounts_singer_stacked = twisst2.stack_topocounts(topocounts_list_singer)


with open(f"{run_name}_{mut_name}.singer.topocounts.tsv", "wt") as tcFile:
    topocounts_singer_stacked.write(tcFile)

with open(f"{run_name}_{mut_name}.singer.intervals.tsv", "wt") as intsFile:
    topocounts_singer_stacked.write_intervals(intsFile)


##################### argweaver ##################

#argweaver
t = time.time()
tl_list_argweaver = run_argweaver(haps, positions+1, N=100000, rec=recombination_rate, mu=mutation_rate, sequence_start=1, sequence_end=seq_len, burnin=500, sample_n_args=10, spacing=10, unphased=False)
time.time() - t

topocounts_list_argweaver = [twisst2.get_topocounts(tl.trees(), leaf_groups=groups_haploid, max_subtrees=total_subtrees, unrooted=False) for tl in tl_list_argweaver]

#unrooted
topocounts_list_argweaver_unrooted = [tc.unroot() for tc in topocounts_list_singer]

[np.mean([weights_error_rate(topocounts, topocounts_argweaver, window_size=x) for topocounts_argweaver in topocounts_list_argweaver]) for x in window_sizes]


topocounts_argweaver_stacked = twisst2.stack_topocounts(topocounts_list_argweaver)

with open(f"{run_name}_{mut_name}.argweaver.topocounts.tsv", "wt") as tcFile:
    topocounts_argweaver_stacked.write(tcFile)

with open(f"{run_name}_{mut_name}.argweaver.intervals.tsv", "wt") as intsFile:
    topocounts_argweaver_stacked.write_intervals(intsFile)


#######################################################################################
######################  bionj windows  ################################################
#######################################################################################

haps_with_anc = np.column_stack([haps, np.zeros(shape=(haps.shape[0], 1))])

bionj_trees = bionj_windows(haps_with_anc, positions, 50, seq_len=seq_len)

t = time.time()
topocounts_bionj = twisst2.get_topocounts(bionj_trees, leaf_groups=groups_haploid, max_subtrees=max_subtrees)
time.time() - t

with open("/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/" + run_name + mut_name + mod_name + ".bionj50.topocounts.tsv", "wt") as tcFile:
    topocounts_bionj.write(tcFile)

with open("/home/smarti11/Dropbox/Research/sticcs_and_twisst2/sim/" + run_name + mut_name + mod_name + ".bionj50.intervals.tsv", "wt") as intsFile:
    topocounts_bionj.write_intervals(intsFile)



#######################################################################################
########################  error rates  ################################################
#######################################################################################

window_sizes = [None] + [int(x) for x in [1e3, 2e3, 5e3, 10e3, 20e3, 50e3]]

[weights_error_rate(topocounts, topocounts_tsinfer, window_size=x) for x in window_sizes]
[weights_error_rate(topocounts, topocounts_relate, window_size=x) for x in window_sizes]
[weights_error_rate(topocounts, topocounts_twisst2_list[0], window_size=x) for x in window_sizes]
[weights_error_rate(topocounts, topocounts_twisst2_list[1], window_size=x) for x in window_sizes]
[weights_error_rate(topocounts, topocounts_twisst2_list[2], window_size=x) for x in window_sizes]
[weights_error_rate(topocounts, topocounts_bionj, window_size=x) for x in window_sizes]


[weights_error_rate(topocounts.unroot(), topocounts_tsinfer.unroot(), window_size=x) for x in window_sizes]
[weights_error_rate(topocounts.unroot(), topocounts_twisst2_list[0].unroot(), window_size=x) for x in window_sizes]
[weights_error_rate(topocounts.unroot(), topocounts_twisst2_list[1].unroot(), window_size=x) for x in window_sizes]
[weights_error_rate(topocounts.unroot(), topocounts_twisst2_list[2].unroot(), window_size=x) for x in window_sizes]
[weights_error_rate(topocounts.unroot(), topocounts_bionj.unroot(), window_size=x) for x in window_sizes]


########################################################################################
##################  loop over parameters to compare tools  ########################
########################################################################################

from multiprocessing import Pool

def infer_and_get_dists(ts, n, N, rec, mu, n_samples, window_sizes,
                                infer_pol=False, gt_error=False,
                                missing_rate=0,
                                max_subtrees_sticcstack=None, seed=None):
    
    result = []
    
    n_groups = len(n_samples)
    
    groups_haploid = twisst2.make_numeric_groups(list(n_samples.values()))
    
    total_subtrees = np.prod(list(n_samples.values()))
    
    topocounts = twisst2.get_topocounts(ts.trees(), leaf_groups = groups_haploid, max_subtrees=total_subtrees, unrooted=False)
    
    topocounts.intervals[:,0] += 1
    
    topocounts_unrooted = topocounts.unroot()
    
    #add mutations
    mts = msprime.sim_mutations(ts, rate = mu, random_seed=seed)
    
    if infer_pol:
        positions, haps = make_haplotype_array(mts, columns=range(n), OG_column_for_repolarisation=n)
    else:
        positions, haps = make_haplotype_array(mts, columns=range(n))
    
    #if len(positions) < 10: #if there are too few SNPs
        #return [np.nan]*189
    
    #remove any SNP at position zero since argweaver can't handle that
    if positions[0] == 0:
        positions = positions[1:]
        haps = haps[1:,:]
    
    #add errors if desired
    if gt_error:
        add_genotyping_errors(haps, error_model, seed)
    
    #if missing data is desired
    if missing_rate > 0:
        add_missing_data(haps, missing_rate=missing_rate, mean_run=2)
    
    seq_len = mts.sequence_length
    
    #sticcs
    ### sticcsStack       
    for forced_ploidy in [1,2,4]:
        #check that you have a multiple of forced ploidy haplotypes to work with!
        assert n % forced_ploidy == 0
        
        ploidies = np.array([forced_ploidy]*int(n/forced_ploidy)*n_groups)
        
        t = time.time()
        
        der_counts, der_positions = sticcs.get_dac_from_haploid_matrix(haps, positions+1, ploidies)
        
        groups = twisst2.make_numeric_groups([n//forced_ploidy for n in n_samples.values()])
        
        topocounts_twisst2 = twisst2.get_topocounts_stacking_sticcs(der_counts, der_positions, ploidies, groups,
                                                                    chrom_start=1, chrom_len = seq_len, max_subtrees=max_subtrees_sticcstack, unrooted=False,
                                                                    second_chances=True, random_seed=seed, silent=True)
        elapsed = time.time() - t
        
        result += [weights_error_rate(topocounts, topocounts_twisst2, window_size=x) for x in window_sizes]
        result += [weights_error_rate(topocounts_unrooted, topocounts_twisst2.unroot(), window_size=x) for x in window_sizes]
        
        result.append(elapsed)
    
    ##################### tsinfer #################
    t = time.time()
    ts_tsinfer = run_tsinfer(haps, positions, seq_len)
    
    topocounts_tsinfer = twisst2.get_topocounts(ts_tsinfer.trees(), leaf_groups=groups_haploid, max_subtrees=total_subtrees, unrooted=False)
    
    topocounts_tsinfer.intervals[:,0] += 1
    
    result += [weights_error_rate(topocounts, topocounts_tsinfer, window_size=x) for x in window_sizes]
    result += [weights_error_rate(topocounts_unrooted, topocounts_tsinfer.unroot(), window_size=x) for x in window_sizes]
    elapsed = time.time() - t
    result.append(elapsed)
    
    ##################### relate #################
    relate_bin="/home/smarti11/software/relate_v1.2.2_x86_64_static/bin/"
    t = time.time()
    ts_relate = run_relate(haps, positions+1, N=N, mu=mu, rec=rec, relate_bin=relate_bin)
    topocounts_relate = twisst2.get_topocounts(ts_relate.trees(), leaf_groups=groups_haploid, max_subtrees=total_subtrees, unrooted=False)
    elapsed = time.time() - t
    
    # relate tree sequence intervals need to be corrected to be 1-based inclusive
    topocounts_relate.intervals[:,0] += 1
    
    #round to make all integers
    topocounts_relate.intervals = np.floor(topocounts_relate.intervals)
    
    #and extend last interval to the sequence length
    topocounts_relate.intervals[-1,1] = seq_len
    
    result += [weights_error_rate(topocounts, topocounts_relate, window_size=x) for x in window_sizes]
    result += [weights_error_rate(topocounts_unrooted, topocounts_relate.unroot(), window_size=x) for x in window_sizes]
    result.append(elapsed)
    
    ##################### singer ##################
    singer_path = "/home/smarti11/software/singer-0.1.8-beta-linux-x86_64/"
    t = time.time()
    ts_list_singer = run_singer(haps, positions, N=N, rec=rec, mu=mu, sequence_start=0, sequence_end=seq_len, singer_path=singer_path, burnin=500, sample_n_args=10, spacing=20)
    #Need to average over all inferred args
    topocounts_list_singer = [twisst2.get_topocounts(ts.trees(), leaf_groups=groups_haploid, max_subtrees=total_subtrees, unrooted=False) for ts in ts_list_singer]
    elapsed = time.time() - t
    
    #adjust zero-based to 1-based
    for tc in topocounts_list_singer:
        tc.intervals[:,0] += 1
    
    topocounts_list_singer_unrooted = [tc.unroot() for tc in topocounts_list_singer]
    
    result += [np.mean([weights_error_rate(topocounts, topocounts_singer, window_size=x) for topocounts_singer in topocounts_list_singer]) for x in window_sizes]
    
    result += [np.mean([weights_error_rate(topocounts_unrooted, topocounts_singer_unrooted, window_size=x) for topocounts_singer_unrooted in topocounts_list_singer_unrooted]) for x in window_sizes]
    
    result.append(elapsed)
    
    ##################### argweaver ##################
    
    #argweaver phased
    t = time.time()
    tl_list_argweaver = run_argweaver(haps, positions+1, N=N, rec=rec, mu=mu, sequence_start=1, sequence_end=seq_len, burnin=500, sample_n_args=10, spacing=10, unphased = False)
    
    topocounts_list_argweaver = [twisst2.get_topocounts(ts.trees(), leaf_groups=groups_haploid, max_subtrees=total_subtrees, unrooted=False) for ts in tl_list_argweaver]
    
    elapsed = time.time() - t
    #Need to average over all inferred args
    topocounts_list_argweaver_unrooted = [tc.unroot() for tc in topocounts_list_argweaver]
    
    result += [np.mean([weights_error_rate(topocounts, topocounts_argweaver, window_size=x) for topocounts_argweaver in topocounts_list_argweaver]) for x in window_sizes]
    
    result += [np.mean([weights_error_rate(topocounts_unrooted, topocounts_argweaver_unrooted, window_size=x) for topocounts_argweaver_unrooted in topocounts_list_argweaver_unrooted]) for x in window_sizes]
    
    result.append(elapsed)
    
    #argweaver unphased
    t = time.time()
    tl_list_argweaver = run_argweaver(haps, positions+1, N=N, rec=rec, mu=mu, sequence_start=1, sequence_end=seq_len, burnin=500, sample_n_args=10, spacing=10, unphased=True)
    
    topocounts_list_argweaver = [twisst2.get_topocounts(ts.trees(), leaf_groups=groups_haploid, max_subtrees=total_subtrees, unrooted=False) for ts in tl_list_argweaver]
    
    elapsed = time.time() - t
    #Need to average over all inferred args
    topocounts_list_argweaver_unrooted = [tc.unroot() for tc in topocounts_list_argweaver]
    
    result += [np.mean([weights_error_rate(topocounts, topocounts_argweaver, window_size=x) for topocounts_argweaver in topocounts_list_argweaver]) for x in window_sizes]
    
    result += [np.mean([weights_error_rate(topocounts_unrooted, topocounts_argweaver_unrooted, window_size=x) for topocounts_argweaver_unrooted in topocounts_list_argweaver_unrooted]) for x in window_sizes]
    
    result.append(elapsed)
    
    ######## bionj windows
    
    haps_with_anc = np.column_stack([haps, np.zeros(shape=(haps.shape[0], 1))])
    
    t = time.time()
    bionj_trees = bionj_windows(haps_with_anc, positions+1, 50, seq_len=seq_len)
    
    topocounts_bionj = twisst2.get_topocounts(bionj_trees, leaf_groups=groups_haploid, max_subtrees=total_subtrees)
    elapsed = time.time() - t
    
    topocounts_bionj = topocounts_bionj.simplify(fill_gaps=True, seq_start=1, seq_end = seq_len)
    
    result += [weights_error_rate(topocounts, topocounts_bionj, window_size=x) for x in window_sizes]
    result += [weights_error_rate(topocounts_unrooted, topocounts_bionj.unroot(), window_size=x) for x in window_sizes]
    
    result.append(elapsed)
    
    return result

#############################################################################################################

threads = 40
nreps=40

n=20
graph = admix_model()
n_samples={"A": n, "B": n, "C":n, "D":n}
seed=123
seq_len=1e5

param_names = ["mutation_rate", "recombination_rate"]

combos = np.array([[1e-8, 1e-8],
                   [1e-9, 1e-8],
                   [1e-8, 1e-9]])

#combos = np.array([[1e-8, 1e-8, 1e-2,0.05,0.05]])

n_combos = len(combos)

n_params = combos.shape[1]

scales = [int(x) for x in [0, 1e3, 2e3, 5e3, 10e3, 20e3, 50e3]]

methods = ["sticcs", "sticcs_dip", "sticcs_tet",
           "tsinfer", "relate", "singer",
           "argweaver", "argweaver_dip",
           "NJ50"]


results = []

for combo_number,combo in enumerate(combos):
    
    print(f"Parameter combination {combo_number + 1} out of {n_combos}...", file=sys.stderr)
    
    mutation_rate, recombination_rate = combo
    
    print(f"Simulating...", file=sys.stderr)
    
    ts_reps = list(msprime.sim_ancestry(samples=n_samples,
                                   demography = msprime.Demography.from_demes(graph),
                                   sequence_length = seq_len,
                                   recombination_rate = recombination_rate, ploidy=1,  #we always simulate haploids. We can make any ploidy in the genotype data later
                                   random_seed=seed, num_replicates=nreps))
    
    # Make a temporary function to do the inferring and getting distances for each
    # This is because pool.map requires a function with only 1 argument
    print(f"Inferring tree sequences and computing distances...", file=sys.stderr)
    
    def _infer_and_get_dist_(i):
        return infer_and_get_dists(ts = ts_reps[i], n = sum(n_samples.values()),
                                  N=100000, rec=recombination_rate, mu=mutation_rate,
                                  n_samples=n_samples, window_sizes=scales,
                                  infer_pol=False, gt_error=False,
                                  missing_rate=0,
                                  max_subtrees_sticcstack=512, seed=seed)
    
    with Pool(threads) as p:
        result_reps = np.array(p.map(_infer_and_get_dist_, range(nreps)))
    
    results.append(result_reps.mean(axis=0))



header = param_names[:]

for method in ["sticcstack", "sticcstack_dip", "sticcstack_tet", "tsinfer", "relate", "singer", "argweaver", "argweaver_dip", "NJ50"]:
    header += [f"{method}_rooted_scale{int(scale)}" for scale in scales]
    header += [f"{method}_unrooted_scale{int(scale)}" for scale in scales]
    header.append(f"{method}_runtime")


output = np.column_stack([combos[:3,], results])

output_name = "/mnt/loki/martin/simon/sticcs_twisst2/sim/admix_n20_weights_accuracy_all.csv"

np.savetxt(output_name, output, fmt='%s', delimiter = ",", header=",".join(header), comments="")


