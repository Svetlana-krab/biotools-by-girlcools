# BioTools-by-GirlCools

Учебный проект: Python-пакет для чтения и анализа биологических форматов данных (FASTA, FASTQ, SAM, VCF) без использования профильных библиотек (Biopython, pysam и др.).

## Структура проекта

На данный момент реализован базовый каркас и парсер FASTA.

* `biotools/` — исходный код пакета.
  * `__init__.py` — инициализация пакета.
  * `base.py` — абстрактные классы (`BioFile`, `SequenceFile`, `IntervalFile`), задающие интерфейс.
  * `fasta.py` — класс `FastaReader` для потокового чтения FASTA через генераторы, класс `Seq` для анализа последовательностей.
* `data/` — тестовые файлы (например, `test.fasta`).
* `demo.py` — скрипт для демонстрации работы `FastaReader`.
* `requirements.txt` — зависимости проекта.
* `.gitignore` — исключения для Git.

## Установка и запуск

1. Клонирование и переход в папку:
   ```bash
   git clone https://github.com/Svetlana-krab/biotools-by-girlcools.git
   cd biotools-by-girlcools
