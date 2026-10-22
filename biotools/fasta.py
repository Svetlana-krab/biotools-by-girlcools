from biotools.base import SequenceFile


class Seq:
    """
    Класс для хранения и анализа биологической последовательности.

    Атрибуты:
        NUCLEOTIDES (set): Множество символов нуклеотидов (A, C, G, T, U, N).
        PROTEIN (set): Множество символов 20 стандартных аминокислот.
        _sequence (str): Последовательность в верхнем регистре.
        _header (str): Заголовок FASTA-записи.
    """

    NUCLEOTIDES = set("ACGTUN")
    PROTEIN = set("ACDEFGHIKLMNPQRSTVWY")

    def __init__(self, sequence, header=""):
        """
        Инициализирует объект последовательности.

        Аргументы:
            sequence (str): Последовательность (ДНК, РНК или белок).
            header (str, опционально): Заголовок FASTA-записи. По умолчанию "".

        Исключения:
            ValueError: Если последовательность пустая.
        """
        if not sequence:
            raise ValueError("Последовательность не может быть пустой")
        self._sequence = sequence.upper().strip()
        self._header = header.strip()

    @property
    def sequence(self):
        """
        Возвращает последовательность в верхнем регистре.

        Возвращает:
            str: Сама последовательность.
        """
        return self._sequence

    @property
    def header(self):
        """
        Возвращает заголовок FASTA-записи.

        Возвращает:
            str: Заголовок записи.
        """
        return self._header

    def __len__(self):
        """
        Возвращает длину последовательности.

        Возвращает:
            int: Количество символов.
        """
        return len(self._sequence)

    def get_alphabet(self):
        """
        Определяет тип алфавита последовательности.

        Возвращает:
            str: "nucleotide", если все символы нуклеотиды; 
                 "protein", если все символы аминокислоты;
                 "unknown", если последовательность содержит смешанные символы.
        """
        symbols = set(self._sequence)
        if symbols.issubset(self.NUCLEOTIDES):
            return "nucleotide"
        if symbols.issubset(self.PROTEIN):
            return "protein"
        return "unknown"

    def gc_content(self):
        """
        Вычисляет GC-состав (долю гуанина и цитозина) в процентах.

        Возвращает:
            float: Процентное содержание GC (от 0 до 100).

        Исключения:
            ValueError: Если последовательность не является нуклеотидной.
        """
        if self.get_alphabet() != "nucleotide":
            raise ValueError("Только для нуклеотидов")
        g = self._sequence.count("G")
        c = self._sequence.count("C")
        return (g + c) / len(self) * 100

    def __repr__(self):
        """
        Возвращает строковое представление объекта для отладки.

        Возвращает:
            str: Строка вида Seq(header='...', length=N).
        """
        return f"Seq(header='{self._header}', length={len(self)})"


class FastaReader(SequenceFile):
    """
    Класс для потокового чтения FASTA-файлов с использованием генератора.
    """

    def __iter__(self):
        """
        Генератор: читает файл построчно и возвращает объекты Seq по одному.
        
        Благодаря использованию yield, файл не загружается в память целиком. 
        Метод возвращает указатель файла в начало (seek(0)) перед чтением, 
        чтобы можно было вызывать итератор многократно.

        Yields:
            Seq: Объект последовательности для каждой FASTA-записи.
        """
        self._file.seek(0)
        header = None
        lines = []
        
        for line in self._file:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith(">"):
                if header is not None:
                    seq_str = "".join(lines)
                    if seq_str:
                        yield Seq(seq_str, header)
                header = line[1:].strip()
                lines = []
            else:
                lines.append(line)
                
        if header is not None:
            seq_str = "".join(lines)
            if seq_str:
                yield Seq(seq_str, header)

    def count(self) -> int:
        """
        Вычисляет общее количество последовательностей в файле.

        Возвращает:
            int: Количество записей FASTA.
        """
        return sum(1 for _ in self)

    def mean_length(self) -> float:
        """
        Вычисляет среднюю длину всех последовательностей в файле.

        Возвращает:
            float: Средняя длина последовательностей (0.0, если файл пуст).
        """
        lengths = [len(rec) for rec in self] 
        return sum(lengths) / len(lengths) if lengths else 0.0