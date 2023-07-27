from operator import ge
import time
from gtfparse import read_gtf
import numpy as np
import pickle
from maxentpy_updated.maxentpy import maxent_fast
from maxentpy_updated.maxentpy.maxent_fast import load_matrix
import BPP.BP_PPT
#from Bio.Seq import Seq
from typing import Tuple, List, Any, Union
from pathlib import Path
from pyfaidx import Fasta

def pickle_load_helper(filename: str) -> Any:
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    return data
    print("printing sys.path!")

def reverse_complement(s) -> str:
   #old_chars = "ACGT"
    #replace_chars = "TGCA"
    #tab = str.maketrans(old_chars,replace_chars)
    #return s.translate(tab)[::-1]
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(s))

#Motif score lookups
MATRIX5 = load_matrix(5)
MATRIX3 = load_matrix(3)
ESE = pickle_load_helper('scores/ese_dict.py')
ESS = pickle_load_helper('scores/ess_dict.py')

#Distance-based susceptibility scores
NEAR_EXON_DON_GAIN_POS_SCORES = np.load("scores/ne_don_gain.npy")
NEAR_EXON_ACC_GAIN_POS_SCORES = np.load("scores/ne_acc_gain.npy")
NEAR_EXON_GAIN_INTERVAL_SIZE = 500
PSEUDOEXON_LENGTHS = np.load("scores/pseudoexon_lengths_convolved_normalised.npy")
PSEUDOEXON_INTERVAL_SIZE = 200

#ANNOTATIONS
ANNOTATION_DEFAULT_HG38 = "gencode.v38.basic.annotation"
ANNOTATION_DEFAULT_HG19 = "gencode.v38lift37.basic.annotation"
CURRENT_PATH = Path.cwd()
ANNOTATIONS_PATH = Path.cwd() / "annotations"

#Other Constants
FLANK = 50
BPP_SEQ_LENGTH = 50
NEAR_EXON_DON_GAIN_DEFAULT = 14
NEAR_EXON_ACC_GAIN_DEFAULT = 10
FEATURE_SCORES_DEFAULT = np.zeros(23, dtype=float)
FEATURE_SCORES_DEFAULT[4] = NEAR_EXON_DON_GAIN_DEFAULT
FEATURE_SCORES_DEFAULT[5] = NEAR_EXON_ACC_GAIN_DEFAULT

def write_output_header(output_file_path: str, mode: str):

    variant_info_columns = ["chr", "start", "ref", "alt", "gene"]#5
    susceptibility_columns = ["don_region", "acc_region", "bp_region", "esr_region", "ne_don_gain_sus", "ne_acc_gain_sus",
                              "di_don_gain_sus_score", "di_don_gain_dist_score", "di_acc_gain_1_sus_score", "di_acc_gain_1_dist_score",
                              "di_acc_gain_sus_2", "di_esr_gain_sus", "di_bp_gain_sus",
                              "bg_pseudoexon_sample_count", "bg_pseudoexon_coverage_sum",
                              "bg_acc_sai_sample_count", "bg_acc_sai_coverage_sum",
                              "bg_don_sai_sample_count", "bg_don_sai_coverage_sum",
                              "bg_acc_mes_sample_count", "bg_acc_mes_coverage_sum",
                              "bg_don_mes_sample_count", "bg_don_mes_coverage_sum"]#23
    motif_score_columns = ["don_loss_alt_mag", "don_loss_decrease", "acc_loss_alt_mag", "acc_loss_decrease",
                           "don_gain_alt_mag", "don_gain_increase", "acc_gain_alt_mag", "acc_gain_increase",
                           "max_bpp_zscore", "bpzsc_ref_minus_alt", "bpp_tray_ref", "bpp_tray_alt",
                           "enhancer_alt_mag", "enhancer_increase", "silencer_alt_mag", "silencer_increase"]#16
    spliceai_columns = []
    if mode == "spliceai":
        spliceai_columns = ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]#8
    if mode == "spliceai_adapted":
        spliceai_columns = ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL", "DS_AG_REF", "DS_AG_ALT",
                            "DS_AL_REF", "DS_AL_ALT", "DS_DG_REF", "DS_DG_ALT", "DS_DL_REF", "DS_DL_ALT"]#16

    text_file = open(output_file_path, "w")
    #text_file.write("\t".join(variant_info_columns + spliceai_columns + susceptibility_columns + motif_score_columns) + "\n")
    text_file.write("\t".join(variant_info_columns + ["score"]) + "\n")
    text_file.close()

def file_setup_annotation_and_fasta(annotation_arg: str, fasta_name: str) -> Tuple[str]:

    annotation_path = Path("annotations/variant_intersect")
    variant_intersect_filename = annotation_path / ("variant_intersect_" + annotation_arg + ".npy")
    background_acceptor_filename = annotation_path / ("background_acceptor_" + annotation_arg + ".npy")
    background_donor_filename = annotation_path / ("background_donor_" + annotation_arg + ".npy")
    gene_name_to_annotation_indices_filename = annotation_path / ("gene_name_to_annotation_indices_" + annotation_arg + ".pkl")
    annotation_indices_intersect_bed_filename = annotation_path / ("annotation_indices_intersect_" + annotation_arg + ".bed")

    #error handling file existence here with if not filename.exists():

    return (variant_intersect_filename, background_acceptor_filename, background_donor_filename,
            gene_name_to_annotation_indices_filename, annotation_indices_intersect_bed_filename)

def get_esr_scores(seq: str) -> Tuple[float, float]:

    ese = np.zeros(len(seq) - 5)
    ess = np.zeros(len(seq) - 5)
    for i in range(len(seq) -  5):
        ese[i] = ESE.get(seq[i : i+6],0.0)
        ess[i] = ESS.get(seq[i : i+6],0.0)
    return ( np.max(ese), np.min(ess) )

def get_maxent_donor(seq: str) -> Tuple[float, float]:    

    mes = np.full( len(seq)-8, -50.0 )
    for i in range( len(seq)-8 ):
        mes[i] = maxent_fast.score5(seq[i : i+9], matrix=MATRIX5)
        #print(seq[i : i+9], mes[i])
    return np.round(mes, 3)

def get_maxent_acceptor(seq: str) -> Tuple[float, int]:  

    mes = np.full( len(seq) - 22, -50.0 )
    for i in range( len(seq)-22 ):
        mes[i] = maxent_fast.score3(seq[i : i+23], matrix=MATRIX3)
        #print(seq[i : i+23], mes[i])
    return np.round(mes, 3)

def get_bpp_scores(bp_ref_seq: str, bp_alt_seq: str, ag_dinucleotide_index: int) -> Tuple[float, int]:  
    print("ref: ", bp_ref_seq)
    print("alt: ", bp_alt_seq)
    if bp_ref_seq == bp_alt_seq:
        print("equal")

    scores = np.zeros(4)
    if set('ACGT') >= set(bp_ref_seq + bp_alt_seq):
        #bp ref score
        (orinp,zbps,zppt,zsc) = BPP.BP_PPT.bppt_get_BPPTsc(bp_ref_seq ,maxL=-1)
        if bool(len(orinp)):
            scores[0] = zsc[0]

        #check if variant is within TRAY motif for predicted branchpoint 
        bp_location_ref = int(orinp[0].split("\t")[1])
        bp_variant_location = ag_dinucleotide_index - 9
        if (bp_location_ref - 3) < bp_variant_location < (bp_location_ref + 2):
            scores[2] = 1

        #bp alt score
        (orinp,zbps,zppt,zsc) = BPP.BP_PPT.bppt_get_BPPTsc(bp_alt_seq,maxL=-1)
        if bool(len(orinp)):
            scores[0] = max(scores[0], zsc[0])
            scores[1] = scores[0] - zsc[0]

        #check if variant is within TRAY motif for predicted branchpoint 
        bp_location_alt = int(orinp[0].split("\t")[1])
        bp_variant_location = ag_dinucleotide_index - 9
        if (bp_location_alt - 3) < bp_variant_location < (bp_location_alt + 2):
            scores[3] = 1
    return scores

def test_speed(fasta, nloops=1000):
    #better to randomly sample VCF than reuse one variant      
    start_time = time.time()
    #for i in range(nloops):
    #    motif_scores = score_motifs("chr1", 66808, "A", "G", "POS", fasta)
    current_time = time.time()
    elapsed_time = current_time - start_time
    print("score_motif time:", elapsed_time)

    start_time = time.time()
    #for i in range(nloops):
    #    region_scores = score_regions("chr1", "A", "POS", *map(int, annotation_string.split(",")) 
        
    current_time = time.time()
    elapsed_time = current_time - start_time
    print("score_motif time:", elapsed_time)
    #return (result)

def score_variant(chr: str, variant_start: int, ref_allele: str, alt_allele: str, fasta: str, annotation_indices: np.ndarray,
                   variant_intersect_annotation: np.ndarray, background_acceptor: np.ndarray, background_donor: np.ndarray,
                   sai_gain_score_and_pos: List):
    '''
        calculate motif-based scores and susceptibility scores for variant:
    '''

    score_output = np.zeros(16, dtype=float)
    feature_output = FEATURE_SCORES_DEFAULT.copy()
    strand = 1
    strand_flip = 1
    deletion_size = max(0, len(ref_allele) - len(alt_allele))
    insertion_size = max(0, len(alt_allele) - len(ref_allele))
    is_indel = np.sign(deletion_size + insertion_size)

    #retrieve sequence
    if annotation_indices[0] == 1:
        ref_seq = str(fasta[chr][variant_start-81:variant_start+80+abs(len(ref_allele)-len(alt_allele)) ])
        contigs = [ ref_seq[0 : 30 + FLANK], alt_allele, ref_seq[30 + FLANK + len(ref_allele):] ]
        alt_seq = "".join(contigs)
    else:
        strand = 0
        strand_flip = -1
        #ref_allele = str(Seq(ref_allele).reverse_complement())
        #alt_allele = str(Seq(alt_allele).reverse_complement())
        ref_allele = reverse_complement(ref_allele)
        alt_allele = reverse_complement(alt_allele)
        neg_strand_deletion_size = deletion_size
        ref_seq = str(fasta[chr][variant_start - 81 + (strand-1)*deletion_size:variant_start+80].reverse.complement)
        contigs = [ ref_seq[0 : 30 + FLANK -neg_strand_deletion_size], alt_allele, ref_seq[30 + (strand-1)*deletion_size + FLANK + len(ref_allele):] ]
        alt_seq = "".join(contigs)
    
    #raw maxentscan scores
    ref_don_str = ref_seq[FLANK +22 +strand*is_indel +(strand-1)*deletion_size : FLANK +39 +(strand-1)*is_indel +strand*deletion_size]
    alt_don_str = alt_seq[FLANK +22 +strand*is_indel +(strand-1)*deletion_size : FLANK +39 +(strand-1)*is_indel +(strand-1)*deletion_size +insertion_size]
    if set('ACGT') >= set(ref_don_str + alt_allele):
        ref_don_scores = get_maxent_donor(ref_don_str)
        alt_don_scores = get_maxent_donor(alt_don_str)

    ref_acc_str = ref_seq[FLANK +8 +strand*is_indel +(strand-1)*deletion_size : FLANK +53 +(strand-1)*is_indel +strand*deletion_size]
    alt_acc_str = alt_seq[FLANK +8 +strand*is_indel +(strand-1)*deletion_size : FLANK +53 +(strand-1)*is_indel +(strand-1)*deletion_size + insertion_size]
    if set('ACGT') >= set(ref_acc_str + alt_allele):
        ref_acc_scores = get_maxent_acceptor(ref_acc_str)
        alt_acc_scores = get_maxent_acceptor(alt_acc_str)

    #intersect variant and susceptibility features
    feature = np.hstack((variant_intersect_annotation[annotation_indices[1]:annotation_indices[2]+1, :],
                         (variant_intersect_annotation[annotation_indices[1]:annotation_indices[2]+1, 1] - variant_start - is_indel)[:,np.newaxis],
                         (variant_start + deletion_size - variant_intersect_annotation[annotation_indices[1]:annotation_indices[2]+1, 0])[:,np.newaxis]))
    #retain only features where both "variant_end_minus_feature_start" and "feature_end_minus_variant_start" are non-negative
    feature = feature[np.where(np.multiply(np.heaviside(feature[:,-1], 1), np.heaviside(feature[:,-2], 1)))]

    #handle feature intersections by type
    for i in range(feature.shape[0]):
        #within donor
        if feature[i,2] == 0:
            feature_output[0] = 1
            if 'ref_don_scores' in locals():
                alt_offset = np.argmax(alt_don_scores[feature[i,-1-strand]- min(deletion_size, feature[i,-1-strand]):feature[i,-1-strand]+1+insertion_size])
                score_output[0] = alt_don_scores[feature[i,-1-strand] - min(deletion_size,feature[i,-1-strand]) + alt_offset]
                score_output[1] = (ref_don_scores[feature[i,-1-strand]]+100) - (score_output[0]+100)
                alt_don_scores[feature[i,-1-strand] - min(deletion_size,feature[i,-1-strand]) + alt_offset] = -100
                ref_don_scores[feature[i,-1-strand]] = -100
        #within acceptor
        elif feature[i,2] == 1:
            feature_output[1] = 1
            if 'ref_acc_scores' in locals():
                alt_offset = np.argmax(alt_acc_scores[feature[i,-1-strand]- min(deletion_size, feature[i,-1-strand]):feature[i,-1-strand]+1+insertion_size])
                score_output[2] = alt_acc_scores[feature[i,-1-strand] - min(deletion_size,feature[i,-1-strand]) + alt_offset]
                score_output[3] = (ref_acc_scores[feature[i,-1-strand]]+100) - (score_output[2]+100)
                alt_acc_scores[feature[i,-1-strand] - min(deletion_size,feature[i,-1-strand]) + alt_offset] = -100
                ref_acc_scores[feature[i,-1-strand]] = -100
        #bp and esr regions
        elif feature[i,2] in [2,3]:
            feature_output[feature[i,2]] = 1
        #near-exon donor gain
        elif feature[i,2] == 4:
            index = abs(strand * NEAR_EXON_GAIN_INTERVAL_SIZE - int(feature[i,-2]))
            feature_output[4] = max(feature_output[4], np.take(NEAR_EXON_DON_GAIN_POS_SCORES, index, mode='clip') * 100)
        #near-exon acceptor gain
        elif feature[i,2] == 5:
            index = abs(strand * NEAR_EXON_GAIN_INTERVAL_SIZE - int(feature[i, -2]))
            feature_output[5] = max(feature_output[5], np.take(NEAR_EXON_ACC_GAIN_POS_SCORES, index, mode="clip") * 100)
        #donor/acceptor gain troughs
        elif feature[i,2] == 6: feature_output[4] = 1
        elif feature[i,2] == 7: feature_output[5] = 1
        #deep_intronic_donor_gain
        elif feature[i,2] == 8:
            feature_output[6] = feature[i,3]
            index = abs(strand * PSEUDOEXON_INTERVAL_SIZE - int(feature[i,-1]))
            feature_output[7] = max(feature_output[8], np.take(PSEUDOEXON_LENGTHS, index, mode="clip") * 100)
        #deep intronic acceptor gain_1
        elif feature[i,2] == 9:
            feature_output[8] = feature[i,3]
            index = abs(strand * PSEUDOEXON_INTERVAL_SIZE - int(feature[i,-1]))
            feature_output[9] = max(feature_output[9], np.take(PSEUDOEXON_LENGTHS, index, mode="clip") * 100)
        #deep intronic acceptor_gain_2, branchpoint gain, esr_gain
        elif feature[i,2] in [10,11,12]:
            feature_output[feature[i,2]] = feature[i,3]
        #background pseudoexon gain
        elif feature[i,2] == 13:
            feature_output[13] = feature[i,3]
            feature_output[14] = feature[i,4]

    #maxentscan gain and background scores
    background_acceptor_intersect = background_acceptor[annotation_indices[5]:annotation_indices[6]+1]
    background_donor_intersect = background_donor[annotation_indices[3]:annotation_indices[4]+1]

    #if variant_start == 22876154:
    #    print("here")

    if 'ref_don_scores' in locals():
        score_output[4] = np.max(alt_don_scores)
        score_output[5] = (score_output[4]+100) - (np.max(ref_don_scores)+100)
        #print(background_donor_intersect[((background_donor_intersect[:,0]-strand_flip) > variant_start - 25) & ((background_donor_intersect[:,0]-strand_flip) < variant_start + 25), :])
        #test_val = background_donor_intersect[:,0]-strand_flip - variant_start
        if strand == 0:
            background_donor_mes_intersect = background_donor_intersect[(background_donor_intersect[:,0]-strand_flip) == (variant_start +6 - np.argmax(alt_don_scores)), :]
            #if background_donor_mes_intersect.size > 0:
            #    print("d",end='')
        if strand == 1:
            background_donor_mes_intersect = background_donor_intersect[(background_donor_intersect[:,0]-strand_flip) == (variant_start -4 + np.argmax(alt_don_scores)), :]
        if background_donor_mes_intersect.size > 0:
            #    print("d",end='')
            feature_output[21:23] = background_donor_mes_intersect[0,1:]
    if 'ref_acc_scores' in locals():
        score_output[6] = np.max(alt_acc_scores)
        score_output[7] = (score_output[6]+100) - (np.max(ref_acc_scores)+100)
        #print(background_acceptor_intersect[((background_acceptor_intersect[:,0]-strand_flip) > variant_start - 25) & ((background_acceptor_intersect[:,0]-strand_flip) < variant_start + 25), :])
        if strand == 1:
            background_acceptor_mes_intersect = background_acceptor_intersect[(background_acceptor_intersect[:,0]-strand_flip) == (variant_start -4 + np.argmax(alt_acc_scores)), :]
        if strand == 0:
            background_acceptor_mes_intersect = background_acceptor_intersect[(background_acceptor_intersect[:,0]-strand_flip) == (variant_start + 5 - np.argmax(alt_acc_scores)), :]
        if background_acceptor_mes_intersect.size > 0:
            #    print("a",end='')
            feature_output[19:21] = background_acceptor_mes_intersect[0,1:]

    if sai_gain_score_and_pos[0] > 0 :
        background_acceptor_sai_intersect = background_acceptor_intersect[(background_acceptor_intersect[:,0]+strand_flip) == (variant_start + sai_gain_score_and_pos[1]), :]
        if background_acceptor_sai_intersect.size > 0:
            feature_output[15:17] = background_acceptor_sai_intersect[0,1:]
    if sai_gain_score_and_pos[2] > 0:
        background_donor_sai_intersect = background_donor_intersect[(background_donor_intersect[:,0]-strand_flip) == (variant_start + sai_gain_score_and_pos[3]), :]
        if background_donor_sai_intersect.size > 0:
            feature_output[17:19] = background_donor_sai_intersect[0,1:]

    #BPP
    ag_dinucleotide_index = ref_seq[FLANK + (strand-1)*deletion_size + 21 : FLANK + (strand-1)*deletion_size + 71 ].find("AG") + 2
    if ag_dinucleotide_index == -1: ag_dinucleotide_index = 50
    if len(ref_allele) < 10:
    #    print(ref_seq)
    #    print(alt_seq)
    #    if ref_seq == alt_seq:
    #        print("pequals")
        print(ref_seq[ FLANK + (strand-1)*deletion_size - 30 : FLANK + (strand-1)*deletion_size + len(ref_allele) + 21 + ag_dinucleotide_index])
        print(alt_seq[ FLANK + (strand-1)*deletion_size - 30 : FLANK + (strand-1)*deletion_size + len(alt_allele) + 21 + ag_dinucleotide_index])
        score_output[8:12] = get_bpp_scores(ref_seq[ FLANK + (strand-1)*deletion_size - 30 : FLANK + (strand-1)*deletion_size + len(ref_allele) + 21 + ag_dinucleotide_index],
                                            alt_seq[ FLANK + (strand-1)*deletion_size - 30 : FLANK + (strand-1)*deletion_size + len(alt_allele) + 21 + ag_dinucleotide_index],
                                            ag_dinucleotide_index)

    #ESRseq
    esr_ref = get_esr_scores( ref_seq [25 + FLANK +strand*is_indel +(strand-1)*deletion_size: 36 + FLANK + (strand-1)*is_indel +strand*deletion_size] )
    esr_alt = get_esr_scores( alt_seq [25 + FLANK + strand*is_indel +(strand-1)*deletion_size: 36 + FLANK + (strand-1)*is_indel + (strand-1)*deletion_size + insertion_size] )

    score_output[12] = esr_alt[0]
    score_output[13] = esr_alt[0] - esr_ref[0]
    score_output[14] = esr_alt[1]
    score_output[15] = (esr_ref[1]+1) - (esr_alt[1]+1)

    return np.concatenate((feature_output, np.round(score_output, 2)), axis = None) 