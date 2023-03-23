#VS extraer consenso bam
#Chido5. Cobertura


import pandas as pd
import argparse
import os
import glob
import pysam
from pyfaidx import Fasta

# leer el archivo con las regiones de interés
# la estructrura es el nobre del cromosoma, inicio, final, nombre del gen y dirección
def read_bed_file(file):
    bed_regions = []
    with open(file, "r") as bed:
        for line in bed:
            if line.strip():
                bed_line = line.strip().split("\t")[:5]
                bed_line[1] = int(bed_line[1])
                bed_line[2] = int(bed_line[2])
                bed_regions.append(bed_line)
    return bed_regions

#IUPAC para concenso
def iupac_code(counter):
    iupac_dict = {
        frozenset(['A', 'C']): 'M',
        frozenset(['A', 'G']): 'R',
        frozenset(['A', 'T']): 'W',
        frozenset(['C', 'G']): 'S',
        frozenset(['C', 'T']): 'Y',
        frozenset(['G', 'T']): 'K',
        frozenset(['A', 'C', 'G']): 'V',
        frozenset(['A', 'C', 'T']): 'H',
        frozenset(['A', 'G', 'T']): 'D',
        frozenset(['C', 'G', 'T']): 'B',
        frozenset(['A', 'C', 'G', 'T']): 'N',
    }

    bases = [base for base, count in counter.items() if count == max(counter.values())]
    iupac_base = iupac_dict.get(frozenset(bases), 'N')
    return iupac_base

#ojo: implementar reverso complemento si detecta -


#secuencia consenso base por base en un pileup

def get_consensus_sequence(alignment_file, reference, chr, start, end, min_percentage, min_coverage, min_base_quality=20):
    ref_seq = reference[chr][start:end].seq
    consensus_seq = []
    total_coverage = 0

    for pileupcolumn in alignment_file.pileup(chr, start, end, min_base_quality=min_base_quality):
        if start <= pileupcolumn.pos < end:
            base_counter = {"A": 0, "C": 0, "G": 0, "T": 0}
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    base_counter[base] += 1

            total_bases = sum(base_counter.values())
            total_coverage += total_bases

            if total_bases < min_coverage:
                consensus_seq.append('N')
                continue

            base_frequencies = {base: count / total_bases for base, count in base_counter.items()}
            consensus_base = "".join([base for base, freq in base_frequencies.items() if freq >= min_percentage])

            if not consensus_base:
                consensus_base = iupac_code(base_counter)

            consensus_seq.append(consensus_base)

    avg_coverage = total_coverage / len(ref_seq) if len(ref_seq) != 0 else 0
    # Si la cobertura promedio es menor que la cobertura mínima requerida, devuelve None en lugar de la secuencia consenso
    return (None if avg_coverage < min_coverage else "".join(consensus_seq)), avg_coverage


def save_coverage_table(coverage_data, output_file):
    with open(output_file, "w") as f:
        f.write("Gene\tSample\tCoverage\n")
        for gene, sample_coverage in coverage_data.items():
            for sample, coverage in sample_coverage.items():
                f.write(f"{gene}\t{sample}\t{coverage:.2f}\n")

#ajustar con referencia

def main(args):
    reference = Fasta(args.reference_file)
    bed_regions = read_bed_file(args.bed_file)

    os.makedirs(args.output_folder, exist_ok=True)

    for region in bed_regions:
        chr, start, end, gene, strand = region

        coverage_file = os.path.join(args.output_folder, f"{gene}_coverage.tsv")
        # Agregar la columna 'Retained' a la cabecera del archivo de cobertura
        with open(coverage_file, "w") as f:
            f.write("Gene\tSample\tCoverage\tRetained\n")

        sum_coverage = 0
        sample_count = 0

        for bam_file in glob.glob(os.path.join(args.bam_folder, "*.bam")):
            sample_name = os.path.splitext(os.path.basename(bam_file))[0]
            try:
                alignment_file = pysam.AlignmentFile(bam_file, "rb")
                consensus_seq, avg_coverage = get_consensus_sequence(alignment_file, reference, chr, start, end, args.min_percentage, args.min_coverage, args.min_base_quality)

                # determinar si la secuancia fue retenida o no
                retained = "Yes" if consensus_seq is not None else "No"

                # Si la secuencia consenso no es None,  cobertura promedio es igual o mayor que la cobertura mínima
                if consensus_seq is not None:
                    strand_suffix = "_F" if strand == "+" else "_R"
                    fasta_file = os.path.join(args.output_folder, f"{gene}{strand_suffix}.fasta")

                    with open(fasta_file, "a") as fasta_out:
                        fasta_out.write(f">{sample_name}\n")
                        fasta_out.write(f"{consensus_seq}\n")

                # Escribir la información de cobertura y si la secuancia fue retenida en el archivo de cobertura
                with open(coverage_file, "a") as f:
                    f.write(f"{gene}\t{sample_name}\t{avg_coverage:.2f}\t{retained}\n")

                sum_coverage += avg_coverage
                sample_count += 1

            except ValueError as e:
                print(f"Error procesando {bam_file}: {e}")
                continue

        # Calcular y generar el promedio en el archivo de cobertura para el gen en cada iteración
        average_coverage = sum_coverage / sample_count
        with open(coverage_file, "a") as f:
            f.write("\n")
            f.write(f"{gene}\tPromedio\t{average_coverage:.2f}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extraer secuencias consenso de archivos BAM usando un archivo BED.")
    parser.add_argument("--bed_file", required=True, help="archivo .bed de entrada con regiones a extraer")
    parser.add_argument("--reference_file", required=True, help="Archivo FASTA de genoma referencia")
    parser.add_argument("--bam_folder", required=True, help="Carpeta con contiene archivos BAM")
    parser.add_argument("--output_folder", required=True, help="Carpeta de salida para los archivos FASTA y tablas de cobertura")
    parser.add_argument("--min_percentage", type=float, default=0.7, help="Porcentaje mínimo de lecturas para llamar a una base (dedault: 0.7)")
    # Agregar argumento en español
    parser.add_argument("--min_coverage", type=int, default=10, help="Cobertura mín para considerar una posición en la secuencia consenso (defualt: 10)")
    parser.add_argument("--min_base_quality", type=int, default=20, help="calidad mínima de la base en los reads (default: 20)")

    args = parser.parse_args()
    main(args)

    #pendientes, convertir a F si es -
