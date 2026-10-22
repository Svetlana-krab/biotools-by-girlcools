from biotools.fasta import FastaReader

with FastaReader("data/test.fasta") as reader:
    print(f"Последовательностей: {reader.count()}")
    print(f"Средняя длина: {reader.mean_length():.2f}")
