import streamlit as st
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import random
from io import StringIO

# Custom CSS for Enhanced Webpage Styling
st.markdown("""
    <style>
        body {
            background-color: #1e1e2e;
            color: #ffffff;
            font-family: 'Poppins', sans-serif;
        }
        .title {
            text-align: center;
            font-size: 42px;
            font-weight: bold;
            color: #ffcc00;
            margin-bottom: 20px;
            text-shadow: 2px 2px 5px rgba(255, 204, 0, 0.5);
        }
        .header {
            text-align: center;
            font-size: 28px;
            font-weight: bold;
            color: #00e6e6;
            margin-bottom: 15px;
        }
        .stTextArea textarea {
            background-color: #2e2e3e;
            color: #ffffff;
            font-size: 16px;
            border-radius: 10px;
            padding: 15px;
            border: 2px solid #ffcc00;
        }
        .stButton>button {
            background: linear-gradient(135deg, #ff5733, #c70039);
            color: white;
            font-size: 18px;
            font-weight: bold;
            padding: 12px 25px;
            border-radius: 10px;
            transition: 0.3s;
            border: none;
        }
        .stButton>button:hover {
            background: linear-gradient(135deg, #c70039, #ff5733);
            transform: scale(1.1);
            cursor: pointer;
        }
        .footer {
            text-align: center;
            font-size: 14px;
            margin-top: 40px;
            color: #bbbbbb;
        }
    </style>
""", unsafe_allow_html=True)

# Streamlit App Title
st.markdown("<h1 class='title'>  Phylogenetic Tree Generator </h1>", unsafe_allow_html=True)

# User Input Section
st.markdown("<h2 class='header'>Enter Your Sequences 🔬</h2>", unsafe_allow_html=True)
sequences_input = st.text_area("Enter sequences in FASTA format:",
                               """>Seq1\nATCGTACGATCG\n>Seq2\nATGGTACGATCA\n>Seq3\nATCGTACGCTCG\n""")

# Button to Generate Tree
if st.button("Generate Phylogenetic Tree"):
    try:
        # Convert the input into a StringIO object
        fasta_io = StringIO(sequences_input)

        # Read sequences into an alignment object
        alignment = AlignIO.read(fasta_io, "fasta")

        # Compute distance matrix
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)

        # Construct tree using Neighbor-Joining
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(distance_matrix)


        # Generate random colors for branches
        def get_random_color():
            return (random.random(), random.random(), random.random())


        colors = {clade: get_random_color() for clade in tree.find_clades()}

        # Plot the tree
        fig, ax = plt.subplots(figsize=(6, 6))
        Phylo.draw(tree, axes=ax, label_func=lambda x: x.name,
                   branch_labels=lambda x: round(x.branch_length, 2) if x.branch_length else "")

        for clade in tree.find_clades():
            if clade.name:
                x, y = ax.transLimits.transform((clade.branch_length or 0, 0))
                ax.text(x, y, clade.name, fontsize=12, color=colors[clade], verticalalignment='center',
                        fontweight='bold')

        # Display the tree in Streamlit
        st.pyplot(fig)

    except Exception as e:
        st.error(f"Error: {str(e)}")

# Footer
st.markdown("<p class='footer'>Developed for Advanced Molecular Docking & Phylogenetic Analysis 🚀</p>",
            unsafe_allow_html=True)