from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW


def get_fasta_file():
    # Demander à l'utilisateur de fournir un fichier FASTA
    fasta_file = input("Veuillez entrer le chemin du fichier FASTA : ")
    
    # Vérifier si le fichier existe
    try:
        # Lire le fichier FASTA pour vérifier sa validité
        with open(fasta_file, "r") as file:
            sequences = list(SeqIO.parse(file, "fasta"))
            if not sequences:
                print("Le fichier FASTA ne contient pas de séquences.")
                return None
            print(f"Fichier FASTA chargé avec {len(sequences)} séquences.")
            return sequences
    except FileNotFoundError:
        print("Le fichier n'a pas été trouvé. Veuillez vérifier le chemin.")
        return None

# Exemple d'utilisation
sequences = get_fasta_file()

def run_blast(sequences):
    # Prendre la première séquence pour l'exemple
    sequence = sequences[0].seq
    
    # Exécuter BLAST avec la séquence
    print("Envoi de la séquence à BLAST...")
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    
    # Sauvegarder les résultats dans un fichier XML
    with open("blast_results.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    
    print("Les résultats BLAST ont été enregistrés dans 'blast_results.xml'.")
    return result_handle

# Exemple d'utilisation
if sequences:
    blast_results = run_blast(sequences)

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
with open("blast_results.xml", "r") as result_handle:
    parse_blast_results(result_handle)
