# Biopython BLAST Automation

This repository serves as an educational resource for learning how to use **Biopython**, a powerful library for bioinformatics tasks in Python. The project demonstrates how to work with sequence data, perform online BLAST (Basic Local Alignment Search Tool) queries, and retrieve results programmatically.

## Repository Contents

- **`ls_orchid.fasta`**: A sample FASTA file containing nucleotide sequences to use as input for BLAST queries.
- **`ls_orchid.gbk`**: A GenBank file corresponding to the sequences in `ls_orchid.fasta` for additional analysis and comparison.
- **`tutoriel.py`**: A Python script that:
  - Loads sequences from `ls_orchid.fasta`.
  - Demonstrates how to perform an online BLAST search with Biopython.
  - Retrieves and formats BLAST results.
- **`test_getfasta.py`**: A test script to verify the handling of FASTA files and sequence extraction functionality.
- **`requirements.txt`**: A list of required Python dependencies, including Biopython, to run the scripts.

## Objectives

- Learn how to use Biopython's **`Bio.Blast`** module to automate BLAST searches.
- Understand how to handle sequence data in both FASTA and GenBank formats.
- Gain experience in parsing and interpreting BLAST results programmatically.

## Prerequisites

- Python 3.8 or later.
- Install dependencies using:
  ```bash
  pip install -r requirements.txt
  ```

## Usage

1. Clone the repository and navigate to its directory.  
2. Install the required dependencies listed in the `requirements.txt` file.  
3. Run the `tutoriel.py` script to perform an online BLAST query using the sequences in the provided FASTA file.  
4. Use the `test_getfasta.py` script to test and verify the handling of sequence data in the FASTA file.  
5. Analyze and interpret the BLAST results in your chosen output format, such as XML or plain text.  

## Learning Outcomes

By exploring this repository, you’ll gain hands-on experience with:  

- Automating bioinformatics workflows using Python.  
- Interfacing with NCBI’s online tools programmatically.  
- Parsing and analyzing bioinformatics data in both FASTA and GenBank formats.  
- Writing and testing bioinformatics scripts effectively.  
