"""Демонстрация работы FastqReader."""
import sys
from itertools import islice
from biotools.fastq import FastqReader

path = "data/test.fastq"

with FastqReader(path) as reader:
    print("Файл:", path)
    print("\nПервые 2 рида (работает генератор):")
    for read_id, seq, qual in islice(reader, 2):
        print(" ", read_id, seq[:30] + "...", qual[:30] + "...")

    print("\nКоличество ридов:", reader.count())
    print("Средняя длина:", round(reader.mean_length(), 2))

    print("\nСтрою графики...")
    reader.plot_all("plots/")
    print("Готово: смотрите папку plots/")