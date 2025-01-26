from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW




def run_blast(sequence, program, database, evalue, id): #TODO
    # Prendre la première séquence pour l'exemple
    sequences = list(SeqIO.parse(sequence, "fasta"))
    sequence_ = sequences[0].seq
    
    # Exécuter BLAST avec les paramètres spécifiés (Pour la première séquence)
    result_handle = NCBIWWW.qblast(program, database, sequence_, expect=evalue)
    
    # Sauvegarder les résultats dans un fichier XML
    with open("blast_results.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    
    print("Les résultats BLAST ont été enregistrés dans 'blast_results.xml'.")
    return result_handle


def parse_blast_results(blast_result_handle):
    # Lire et parser le fichier XML
    blast_records = NCBIXML.parse(blast_result_handle)
    
    # Parcourir les résultats de BLAST
    for blast_record in blast_records:
        print(f"Query: {blast_record.query}")
        
        # Vérifier si des hits ont été trouvés
        if blast_record.alignments:
            for alignment in blast_record.alignments:
                print(f"  Hit: {alignment.hit_id}")
                print(f"  Description: {alignment.hit_def}")
                
                # Vérifier si 'hit_len' est disponible
                if hasattr(alignment, 'hit_len'):
                    print(f"  Length: {alignment.hit_len}")
                else:
                    print("  Length: non disponible")
                
                # Afficher les top alignments
                for hsp in alignment.hsps:
                    print(f"    E-value: {hsp.expect}")
                    print(f"    Score: {hsp.score}")
                    print(f"    Identity: {hsp.identities}/{hsp.align_length} ({(hsp.identities / hsp.align_length) * 100:.2f}%)")
                    print(f"    Query start: {hsp.query_start}, Query end: {hsp.query_end}")
                    print(f"    Hit start: {hsp.sbjct_start}, Hit end: {hsp.sbjct_end}")
                    print("-" * 50)
        else:
            print("  Aucun hit trouvé.")
            print("-" * 50)


# Exemple d'utilisation
if __name__ == "__main__":
    sequences = get_fasta_file()
    if sequences:
        program, database, evalue = get_blast_parameters()
        blast_results = run_blast(sequences, program, database, evalue)
        
        # Lire les résultats BLAST enregistrés pour les analyser
        with open("blast_results.xml", "r") as result_handle:
            parse_blast_results(result_handle)
