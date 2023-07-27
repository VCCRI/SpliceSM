import argparse
import pandas as pd
from gtfparse import read_gtf
import numpy as np
from pyfaidx import Fasta
import utils
import pysam
import logging
from pathlib import Path
from xgboost import XGBClassifier
import sys

np.set_printoptions(suppress=True)

#Other Constants
FEATURE_COUNT_EXCLUDING_SAI = 39
 
def arr_to_str(arr): return '\t'.join([str(num) for num in arr])

def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fasta", required=True,
                    help="path to FASTA")
    ap.add_argument("-a", "--annotation", required=True, 
                     help='specify "hg38" or "hg19", or path to GTF')
    ap.add_argument("-i", "--input", required=True, 
                     help="path to VCF file(s)")
    ap.add_argument("-o", "--output", required=True, 
                     help="output file path and name")
    ap.add_argument("-m", "--mode", required=False, default="spliceai", 
                     help="output file path and name")
    ap.add_argument("-b", "--batch_size", required=False, default="300",
                     help="xgboost batch size")


    args = ap.parse_args()
    if None in [args.fasta, args.input, args.output, args.annotation]:
        logging.error('Usage: run.py [-h] [-i [input]] [-o [output]] [-a [annotation]]')
        exit()
    
    setup_filenames = utils.file_setup_annotation_and_fasta(args.annotation, args.fasta)
    background_acceptor = np.load(setup_filenames[1])
    background_donor = np.load(setup_filenames[2])
    variant_intersect_annotation = np.load(setup_filenames[0])
    gene_name_to_annotation_indices = utils.pickle_load_helper(setup_filenames[3])
    fasta = Fasta(args.fasta, rebuild=False)
    utils.write_output_header(args.output, args.mode)
    model = XGBClassifier()
    batch_size = int(args.batch_size)

    if args.mode == "spliceai_free":
        #handle if spliceai free with "annotation_indices_intersect_bed" 
        #intersect_file_path = Path("temp" + input_file_name + "_intersect.txt")
        print("handle no spliceai")
    else:
        if args.mode == "spliceai_adapted":
            model.load_model("/home/steven/Documents/Work/Sgen/models/model_ssm_fast.txt")
        try:
            vcf = pysam.VariantFile(args.input)
        except:
#probs = model_load_test.predict_proba(X_test_within_donor[X_columns_ssm_fast])
            "VCF file missing"
        else:
            no_score_count = 0
            out_df = pd.DataFrame(index=range(batch_size),columns=["chr", "pos", "ref", "alt", "gene"])
            variant_info_batch = np.zeros([batch_size,5], dtype = "S20")
            batch_counter = 0
            x = np.zeros([batch_size, FEATURE_COUNT_EXCLUDING_SAI + 8], dtype = float)
            for rec in vcf:
                print(rec)
                #check for VCF lines with missing SpliceAI scores
                try:
                    sai_line_info = rec.info["SpliceAI"]
                except:
                    #print("No score: ", rec.chrom, ":", rec.pos, ":", rec.ref, ":", rec.alts)
                    no_score_count+=1
                else:
                    for sai_variant_info_unsplit in sai_line_info:
                        #Skip if:
                        #missing spliceai scores, ie "."
                        #invalid ref/alt alleles or gene name
                        sai_variant_info = sai_variant_info_unsplit.split(sep='|')
                        #only need this line to chop off previous spliceai-adapted fields
                        sai_variant_info = sai_variant_info[:10]
                        annotation_indices = gene_name_to_annotation_indices.get(sai_variant_info[1], "missing")
                        #Complex delins are given lines like... chr2    69152165        .       GCA     AAT     .       .       SpliceAI=AAT|ANTXR1|.|.|.|.|.|.|.|.
                        if sai_variant_info[-1] == "." or set(rec.ref + rec.alts[0]) >= set('ACGT') or annotation_indices == "missing" or "alt" in rec.chrom:
                            no_score_count+=1
                        else:
                            variant_info = [rec.chrom, rec.pos, rec.ref, rec.alts[0]] + sai_variant_info[1:]
                            if len(rec.ref) > 40:
                                
                                no_score_count+=1
                            else:
                                variant_info_batch[batch_counter,:] = np.asarray([rec.chrom, rec.pos, rec.ref, rec.alts[0], sai_variant_info[1]])
                                out_df.iloc[batch_counter, :] = [rec.chrom, rec.pos, rec.ref, rec.alts[0], sai_variant_info[1]]
                                #print("DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL|DS_AG_REF|DS_AG_ALT|DS_AL_REF|DS_AL_ALT|DS_DG_REF|DS_DG_ALT|DS_DL_REF|DS_DL_ALT")
                                #print(variant_info + [gene_name_to_annotation_indices[variant_info[4]][0]])

                                a = np.asarray(sai_variant_info[2:])
                                if variant_info[4] == "AL627309.1":
                                    print("here")
                                b = utils.score_variant(*variant_info[:4], fasta, annotation_indices,
                                                        variant_intersect_annotation, background_acceptor, background_donor,
                                                        [float(sai_variant_info[i]) for i in [2,4,6,8]])




                                x[batch_counter,:] = np.concatenate([np.asarray(sai_variant_info[2:]),
                                                                     utils.score_variant(*variant_info[:4], fasta, gene_name_to_annotation_indices[variant_info[4]],
                                                                                         variant_intersect_annotation, background_acceptor, background_donor,
                                                                                        [float(sai_variant_info[i]) for i in [2,4,6,8]])])

                                #x = x.fillna(0)
                                #x[x == -1] = 0 ###x = x.replace({-1.00: 0.00})
                                
                                batch_counter+=1
                                if batch_counter == batch_size:
                                    batch_counter = 0

                                    print(rec)
                                    with open("output/test_trunc.txt", 'r') as fp:
                                        q = len(fp.readlines())
                                        print('Total lines:', q) # 8

                                    #This grabs the columns needed for SpliceAI free model, in correct order. Later, want columns ordered the same so this is not so convoluted.
                                    #probs = model.predict_proba(x[:, np.r_[39:55, 16:31, 35:39, 31:35, 0:4]])[:,1]
                                    probs = model.predict_proba(x[:, np.r_[31:47, 8:23, 27:31, 23:27, 0:4]])[:,1]
                                    #unknown = x[:, np.r_[31:47, 8:23, 27:31, 23:27, 0:4]]
                                    #print(unknown)
                                    with open(args.output, "a") as fp:
                                        pd.concat([out_df, pd.Series(np.round(probs, 4))], axis=1).to_csv(fp, sep="\t", header=None, index=False)

                                #with open("x_out_test.pkl", "wb") as fp:
                                #    pickle.dump(x[0,:], fp)
                                #print(x[0,:])
                                #include features in output
                                #a = "\t".join([str(i) for i in (variant_info[:5] + x[batch_counter,:].tolist())]) + "\n"
            #final batch (for when n_lines not divisible by batch size)
            if batch_counter > 0:
                probs = model.predict_proba(x[:, np.r_[31:47, 8:23, 27:31, 23:27, 0:4]])[:,1]
                with open(args.output, "a") as fp:
                    pd.concat([out_df, pd.Series(np.round(probs, 4))], axis=1).to_csv(fp, sep="\t", header=None, index=False)

            #print("no_score_count: ", no_score_count)

    # Pysam now working for SpliceAI inputs (where gene is already known, so no need for intersection)
    # Pysam option is much faster. But for code simplicity, it may be best to intersect in all cases
    ############################################################################
    #Bedtools intersect option needed for no SpliceAI.

#    #intersect variants with genes
#    vcf = pybedtools.BedTool(args.input)
#    vcf.intersect(str(setup_filenames[3]), wb=True).saveas(intersect_file)
#        
#    if not os.path.isfile(intersect_file):
#            print("Intersect file missing")
#    try:
#        with open(intersect_file) as fp:
#            for line in fp:
#                #columns d(1-8) are not needed
#                chr, start, d1, ref, alt, d2, d3, d4, d5, d6, d7, annotation_string, gene_id, strand = line.split()
#                start = int(start)
#                gene_name = gene_id_name_dict.get(int(gene_id), "-")
#                sai_info_split = d4.split(",")
##     #retrieve reference allele from fasta to check for mismatch with provided ref allele
##     fasta_ref_allele = s[22 + FLANK, 22 + FLANK + len(ref)].toUpperCase()
##     #check fasta reference matches provided ref allele
##     if not fasta_ref_allele == ref:
    
                #print(chr, str(start), ref, alt)

#                if len(ref) < 45:
#                    output_scores = utils.score_variant(chr, start, ref, alt, strand, fasta, annotation_regions, annotation_transcripts, enst_ids,
#                                                  transcript_boundaries, *map(int, annotation_string.split(",")) )
                        #acc_site_pos = 2
                        #don_site_pos = -100
                        #need to pass relative position of acceptor/donor sites.
                        #complicated when there are different transcript groups, encoded in region_scores info
                        #For simplicity now, just repeat all motif_score calculations fro different transcript groups,
                        #Even though most are redundant.
                        
                        #motif_scores = utils.score_motifs(chr, start, ref, alt, strand, fasta, acc_site_pos, don_site_pos)
    #MOTIF_SCORE_NAMES = ["donor_1", "donor_2", "acceptor_1", "acceptor_2", "donor_pos_1", "donor_pos_2",
    #                     "acceptor_pos_1", "acceptor_pos_2", "enhancer", "silencer", "branchpoint"]
    #INITIAL_SCORES = np.array([[-50.0,-50.0,-50.0,-50.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    #                           [-50.0,-50.0,-50.0,-50.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]])
                    #max_don_score, max_acc_score, max_bpp_zscore
                    #don_ref_minus_alt, acc_ref_minus_alt, ese_ref_minus_alt, ess_ref_minus_alt 
    
                        #motif_score_output = np.round(np.array([max(motif_scores[0,0], motif_scores[1,0]),
                        #                        max(motif_scores[0,2], motif_scores[1,2]),
                        #                        max(motif_scores[0,10], motif_scores[1,10]),
                        #                        motif_scores[0,0]-motif_scores[1,0],
                        #                        motif_scores[0,2]-motif_scores[1,2], 
                        #                        motif_scores[0,8]-motif_scores[1,8], 
                        #                        motif_scores[0,9]-motif_scores[1,9], 
                        #                        motif_scores[0,10]-motif_scores[1,10],
                        #                        motif_scores[0,11], motif_scores[1,11]]) ,2)
                        
                        #For training, use combined feature intersections from any transcript
                        #combined_transcript_scores = region_scores[0][0].flatten()
                        #print(len(region_scores))
                        #print(combined_transcript_scores)
                        #for i in range(1,len(region_scores)):
                        #    combined_transcript_scores = np.maximum(combined_transcript_scores, region_scores[i][0].flatten())
                        #    print(region_scores[i][0].flatten())
                        #    print(combined_transcript_scores)
                        
                        #print(combined_transcript_scores)
                        #out_str = chr + "\t" + str(start) + "\t" + ref + "\t" + alt + "\t" + gene_name + "\t" + strand + "\t" + arr_to_str(combined_transcript_scores) + "\t" + arr_to_str(motif_scores) + "\n"
#                        out_str = chr + "\t" + str(start) + "\t" + ref + "\t" + alt + "\t" + gene_name + "\t" + strand + "\t" + arr_to_str(output_scores) + "\n"
#                        text_file = open(output_file, "a")
#                        n = text_file.write(out_str)
#                        text_file.close()
                    #print(arr_to_str(region_scores))
    
#        except IOError:
#            print("Could not read intersect file")

    
if __name__ == "__main__":
    main()

    #a = pybedtools.BedTool(args.input)
    #a.intersect('annotations/gencode.v36.basic.genes.pos_strand.gtf', wb=True).saveas('pos_intersect.bed')
    #consider using bedtools "cut" here. And saveas(as_text option). And pass to a function that expects an iterator instead of saving to file?
    #also maybe check_genome(**kwargs)?
