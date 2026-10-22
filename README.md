# BioTools-by-GirlCools
Учебный проект: Python-пакет для чтения и анализа биологических форматов данных (FASTA, FASTQ, SAM, VCF) без использования профильных библиотек (Biopython, pysam и др.).

## Структура проекта
Реализованы базовый каркас, а также парсеры FASTA, FASTQ и VCF.

* `biotools/` — исходный код пакета.
  * `__init__.py` — инициализация пакета.
  * `base.py` — абстрактные классы (`BioFile`, `SequenceFile`, `IntervalFile`, `GenomicDataReader`), задающие интерфейс.
  * `fasta.py` — класс `FastaReader` для потокового чтения FASTA через генераторы, класс `Seq` для анализа последовательностей.
  * `fastq.py` — класс `FastqReader` для потокового чтения FASTQ, подсчета базовой статистики (средняя длина, качество) и построения графиков (аналог FastQC) с помощью `matplotlib` и `pandas`.
  * `vcf.py` — класс `VcfReader` для чтения VCF, парсинга заголовков, подсчета статистики по регионам и фильтрации вариантов по заданному отрезку.
* `data/` — тестовые файлы.
* `demo.py` — скрипт для демонстрации работы ридеров.
* `requirements.txt` — зависимости проекта.
* `.gitignore` — исключения для Git.

## Установка и запуск
1. Клонирование и переход в папку:
   git clone https://github.com/Svetlana-krab/biotools-by-girlcools.git
   cd biotools-by-girlcools
2. Создание виртуального окружения и установка зависимостей:
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
3. Запуск демонстрации
python3 demo.py

Примеры использования
1. FASTA
Чтение последовательностей, подсчет их количества и средней длины.

from biotools.fasta import FastaReader
with FastaReader("data/test.fasta") as reader:
    print(f"Последовательностей: {reader.count()}")
    print(f"Средняя длина: {reader.mean_length():.2f}")
    
2. FASTQ
Чтение ридов и построение трех графиков качества (в стиле FastQC): распределение качества по позициям, содержание нуклеотидов, распределение длин. Поддерживает .gz файлы.

from biotools.fastq import FastqReader
with FastqReader("data/test.fastq") as reader:
    stats = reader.plot_all(output_dir="plots")
    print(f"Всего ридов: {stats.n_reads}")
    
3. VCF
Чтение вариантов, парсинг заголовков, статистика по регионам (через pandas) и поиск в заданном отрезке.

from biotools.vcf import VcfReader
with VcfReader("data/test.vcf") as reader:
    print(f"Количество вариантов: {reader.count()}")
    # Статистика по хромосомам в виде DataFrame
    print(reader.get_stats_by_region())
    # Поиск вариантов в отрезке
    variants = reader.filter_by_region(chrom="chr1", start=100, end=500)
    print(f"Найдено вариантов в отрезке: {len(variants)}")
    
