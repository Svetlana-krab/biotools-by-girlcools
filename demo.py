import os
from biotools.fasta import FastaReader
from biotools.vcf import VcfReader
from biotools.fastq import FastqReader
from itertools import islice
from biotools.sam import SamReader

def main():
    data_dir = "data"

    # Проходимся по всем файлам в папке data
    for filename in os.listdir(data_dir):
        filepath = os.path.join(data_dir, filename)
        
        # Нас интересуют только файлы, а не папки
        if not os.path.isfile(filepath):
            continue

        print(f"\n--- Анализ файла: {filename} ---")
        
        # Смотрим на расширение файла
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            with FastaReader(filepath) as reader:
                print(f"Это FASTA файл.")
                print(f"Количество последовательностей: {reader.count()}")
                print(f"Средняя длина: {reader.mean_length():.2f}")

        elif filename.endswith(".fastq") or filename.endswith(".fq"):
            with FastqReader(filepath) as reader:
                print(f"Это FASTQ файл.")
                print("\nПервые 2 ряда (работает генератор):")
                for read_id, seq, qual in islice(reader, 2):
                    print(" ", read_id, seq[:30] + "...", qual[:30] + "...")
                print(f"\nКоличество ридов: {reader.count()}")
                print(f"Средняя длина: {round(reader.mean_length(), 2)}")
                print("\nСтрою графики...")
                # Проверяем, существует ли папка plots, если нет - создаем
                if not os.path.exists("plots"):
                    os.makedirs("plots")
                reader.plot_all("plots/")
                print("Готово: смотрите папку plots/")

        elif filename.endswith(".vcf"):
            with VcfReader(filepath) as reader:
                print(f"Это VCF файл.")
                print(f"Количество вариантов: {reader.count()}")
                print(f"Хромосомы: {reader.get_chromosomes()}")
                print("Статистика по регионам:")
                print(reader.get_stats_by_region())

        elif filename.endswith(".sam"):
            reader = SamReader(filepath)
            print(f"Это SAM файл.")
            print(f"Количество выравниваний: {reader.count()}")
            print(f"Хромосомы: {reader.get_chromosomes()}")
    
            print("\nСтатистика по хромосомам:")
            print(reader.stats_by_chrom())
    
            # Проверка пересечения (берем первый регион из первой хромосомы для примера)
            chroms = reader.get_chromosomes()
            if chroms:
                first_chrom = chroms[0]
                # Ищем выравнивания с 1 по 1000 позицию
                overlaps = list(reader.intersect(first_chrom, 1, 1000))
                print(f"\nНайдено пересечений в {first_chrom}:1-1000: {len(overlaps)}")
        
        else:
            print("Неизвестный формат файла.")

if __name__ == "__main__":
    main()
