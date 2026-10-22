"""Модуль для работы с FASTQ-файлами."""
import gzip
import os
import matplotlib.pyplot as plt
import pandas as pd
from .base import SequenceFile

class FastqStats:
    """Накопитель статистики FASTQ (считает на лету, не хранит риды в памяти)."""
    
    def __init__(self):
        self.n_reads = 0
        self.length_counts = {}  # Словарь {длина: количество}
        self._q_sums = []        # Сумма качеств по позициям
        self._q_counts = []      # Сколько ридов дошли до позиции
        self._b_counts = {b: [] for b in "ACGTN"}  # Счетчики нуклеотидов

    def add(self, seq, qual):
        """Добавляет один рид в статистику."""
        self.n_reads += 1
        length = len(seq)
        self.length_counts[length] = self.length_counts.get(length, 0) + 1

        # Если рид длиннее предыдущих, удлиняем наши списки нулями
        if length > len(self._q_sums):
            diff = length - len(self._q_sums)
            self._q_sums.extend([0] * diff)
            self._q_counts.extend([0] * diff)
            for b in self._b_counts:
                self._b_counts[b].extend([0] * diff)

        # Накапливаем данные для каждой буквы в риде
        for i, (s, q) in enumerate(zip(seq.upper(), qual)):
            self._q_sums[i] += ord(q) - 33  # Phred+33
            self._q_counts[i] += 1
            if s in self._b_counts:
                self._b_counts[s][i] += 1

    @property
    def mean_quality(self):
        """Среднее качество для каждой позиции."""
        return [qs / max(qc, 1) for qs, qc in zip(self._q_sums, self._q_counts)]

    @property
    def base_content(self):
        """Процент нуклеотидов на каждой позиции (pandas.DataFrame)."""
        data = {b: [bc / max(qc, 1) * 100 for bc, qc in zip(self._b_counts[b], self._q_counts)] for b in "ACGTN"}
        df = pd.DataFrame(data)
        df.index = range(1, len(df) + 1)
        return df

class FastqReader(SequenceFile):
    """Ридер FASTQ-файлов (генератор записей, подсчет, графики)."""

    def _open(self):
        """Открывает файл (поддерживает .gz)."""
        if str(self.filename).endswith(".gz"):
            return gzip.open(self.filename, "rt")
        return open(self.filename, "r")

    def __enter__(self):
        """Открывает файл (поддерживает .fastq и .fastq.gz)."""
        self._file = self._open()
        return self

    def __iter__(self):
        """Генератор: читает файл по 4 строки и отдает кортеж (id, seq, qual)."""
        if self._file is None or self._file.closed:
            with self._open() as f:
                yield from self._parse(f)
        else:
            self._file.seek(0)
            yield from self._parse(self._file)

    @staticmethod
    def _parse(f):
        """Разбирает открытый файл на записи по 4 строки.

        Пустые строки между записями и в конце файла пропускаются.
        Если запись повреждена, выбрасывается ValueError с номером записи.
        """
        record_no = 0
        while True:
            header = f.readline()
            if not header:                       
                break
            if not header.strip():               
                continue
            record_no += 1
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            if not header.startswith("@"):
                raise ValueError(f"Запись {record_no}: нет '@' в заголовке")
            if not plus.startswith("+"):
                raise ValueError(f"Запись {record_no}: нет '+' в 3-й строке")
            if len(seq) != len(qual):
                raise ValueError(
                    f"Запись {record_no}: длина последовательности ({len(seq)}) "
                    f"не равна длине качества ({len(qual)})"
                )
            yield header[1:].strip(), seq, qual

    def count(self):
        """Считает количество ридов (через генератор)."""
        return sum(1 for _ in self)

    def mean_length(self):
        """Считает среднюю длину ридов (через генератор)."""
        total = n = 0
        for _, seq, _ in self:
            total += len(seq)
            n += 1
        return total / n if n else 0.0

    def compute_stats(self):
        """Собирает статистику за один проход."""
        stats = FastqStats()
        for _, seq, qual in self:
            stats.add(seq, qual)
        return stats

    def plot_all(self, output_dir=None):
        """Строит 3 графика FastQC за один проход и возвращает статистику."""
        stats = self.compute_stats()
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # 1. Среднее качество по позициям
        plt.figure(figsize=(10, 5))
        plt.plot(range(1, len(stats.mean_quality) + 1), stats.mean_quality)
        plt.title("Per base sequence quality")
        plt.xlabel("Позиция")
        plt.ylabel("Качество (Phred)")
        self._finish(output_dir, "per_base_quality.png")

        # 2. Состав нуклеотидов по позициям
        plt.figure(figsize=(10, 5))
        content = stats.base_content
        for b in "ACGTN":
            plt.plot(content.index, content[b], label=b)
        plt.title("Per base sequence content")
        plt.xlabel("Позиция")
        plt.ylabel("Нуклеотид, %")
        plt.legend()
        self._finish(output_dir, "per_base_content.png")

        # 3. Распределение длин ридов
        plt.figure(figsize=(10, 5))
        plt.bar(stats.length_counts.keys(), stats.length_counts.values())
        plt.title("Sequence length distribution")
        plt.xlabel("Длина (bp)")
        plt.ylabel("Количество ридов")
        self._finish(output_dir, "length_distribution.png")

        return stats

    @staticmethod
    def _finish(output_dir, filename):
        """Сохраняет текущий график в папку или показывает его на экране."""
        plt.tight_layout()
        if output_dir:
            plt.savefig(os.path.join(output_dir, filename), dpi=150)
        else:
            plt.show()
        plt.close()