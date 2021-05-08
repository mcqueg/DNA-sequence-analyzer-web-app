####################
# import libraries #
####################

import pandas as pd
import streamlit as st
from PIL import Image
from Bio.Seq import Seq


#############
# functions #
#############


def nucleotide_count(seq):
    nucleotide_dict = dict([
            ('A', seq.count('A')),
            ('T', seq.count('T')),
            ('G', seq.count('G')),
            ('C', seq.count('C')),
            ])
    return nucleotide_dict


def get_gc_percentage(seq):
    gc_num = 0
    for i in seq:
        if (i == 'G' or i == 'C'):
            gc_num += 1
        else:
            continue
    gc_percentage = (gc_num/len(seq))*100
    return gc_percentage


def get_reverse_complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_comp = "".join(complement_dict.get(base, base)
                           for base in reversed(seq))
    return reverse_comp


###############
# page header #
###############

image = Image.open("dna-sequence-logo.jpg")

st.image(image, use_column_width=True)

st.write("***")

#######################
# sequence input #
#######################
st.header("Enter DNA Sequence")
seq_example = "ATAGATATTGCGATCCGGATATAGC"
sequence = st.text_area(">DNA Query", seq_example).upper()


#######################################
# Nucleotide count & GC/AT percentage #
#######################################
st.header("Results")

st.subheader("Nucleotide count")
N = nucleotide_count(sequence)
df = pd.DataFrame.from_dict(N, orient='index')
df = df.rename({0: 'count'}, axis='columns')
st.write(df)

st.write("Total base pair count: ", len(sequence))
st.write("% GC: ", get_gc_percentage(sequence))
st.write("% AT: ", (100-get_gc_percentage(sequence)))

#################################################
# reverse complement/ transcription/translation #
#################################################
st.subheader("Reverse complement: ")
reverse_complement = get_reverse_complement(sequence)
st.write(reverse_complement)

st.subheader("Transcription: ")
my_dna = Seq(sequence)
transcribed_strand = my_dna.transcribe()
st.write(transcribed_strand)

st.subheader("Translation: ")
translated_strand = transcribed_strand.translate()
st.write(translated_strand)
