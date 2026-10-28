# BioTools-by-GirlCools

Учебный проект: Python-пакет для чтения и анализа биологических форматов данных (**FASTA, FASTQ, SAM, VCF**) без использования профильных библиотек (Biopython, pysam и др.).

## Структура проекта

* `biotools/` — исходный код пакета.
  * `__init__.py` — инициализация пакета.
  * `base.py` — абстрактные классы (`BioFile`, `SequenceFile`, `GenomicDataReader`), задающие интерфейс.
  * `fasta.py` — класс `FastaReader` для потокового чтения FASTA и класс `Seq` для анализа последовательностей.
  * `fastq.py` — класс `FastqReader` для потокового чтения FASTQ, подсчета статистики и построения графиков (аналог FastQC).
  * `vcf.py` — класс `VcfReader` для чтения VCF, парсинга заголовков, статистики по регионам и фильтрации.
  * `sam.py` — класс `SamReader` для чтения SAM, анализа выравниваний и поиска пересечений.
* `data/` — тестовые файлы для всех четырех форматов.
* `docs/` — исходники и сгенерированная HTML-документация (Sphinx).
* `demo.py` — универсальный скрипт для демонстрации работы всех ридеров.
* `requirements.txt` — зависимости проекта.
* `.gitignore` — исключения для Git.

## Установка и запуск

1. Клонирование репозитория:

    ```bash
    git clone https://github.com/Svetlana-krab/biotools-by-girlcools.git
    cd biotools-by-girlcools
    ```
2. Создание виртуального окружения и установка зависимостей:
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    ```
3. Запуск демонстрации:
    ```bash
    python demo.py
    ```

## Документация:
**[Открыть документацию на Github Pages](https://svetlana-krab.github.io/biotools-by-girlcools/)**