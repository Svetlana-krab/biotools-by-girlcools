import pandas
from .base import GenomicDataReader

class SamReader(GenomicDataReader):
    """
    Класс для чтения и анализа файлов формата SAM (Sequence Alignment Map).
    Наследуется от GenomicDataReader.
    """

    def __init__(self, filename):
        """
        Инициализирует ридер SAM.

        Args:
            filename (str): Путь к SAM-файлу.
        """
        super().__init__(filename)
        
    def __iter__(self):
        """
        Возвращает генератор для чтения записей (выравниваний).

        Returns:
            generator: Генератор строк с выравниваниями.
        """
        return self.read_records()
    
    def get_chromosomes(self) -> list:
        """
        Извлекает названия хромосом из заголовков @SQ.

        Returns:
            list: Список названий хромосом.
        """
        chromosomes = []
        for line in self.header("@SQ"):
            fields = line.split("\t")
            for field in fields[1:]:
                if field.startswith("SN:"):
                    chromosomes.append(field.split(":", 1)[1])
        return chromosomes

    def filter_by_region(self, chrom: str, start: int, end: int) -> list:
        """
        Ищет выравнивания в заданном отрезке.

        Args:
            chrom (str): Название хромосомы.
            start (int): Начальная позиция.
            end (int): Конечная позиция.

        Returns:
            list: Список строк с выравниваниями.
        """
        return list(self.intersect(chrom, start, end))

    def read_records(self):
        """
        Генератор для чтения строк с выравниваниями (пропускает заголовки).

        Yields:
            str: Строка с выравниванием.
        """
        with open(self.filename, "r") as file:
            for line in file:
                if line.startswith("@"):
                    continue
                if line.strip() == "":
                    continue
                yield line.rstrip("\n")

    def count(self) -> int:
        """
        Считает общее количество выравниваний в файле.

        Returns:
            int: Количество выравниваний.
        """
        total = 0
        for record in self.read_records():
            total = total + 1
        return total

    def header(self, group=None):
        """
        Читает заголовки файла (строки, начинающиеся с @).

        Args:
            group (str, optional): Конкретная группа заголовков (например, "@SQ"). 
                                   Если не указана, возвращает все группы.

        Returns:
            dict или list: Если group=None, возвращает словарь со всеми группами.
                           Если group указана, возвращает список строк этой группы.
        """
        groups = {}

        with open(self.filename, "r") as file:
            for line in file:
                if not line.startswith("@"):
                    break

                line = line.rstrip("\n")
                fields = line.split("\t")
                group_name = fields[0]

                if group_name not in groups:
                    groups[group_name] = []

                groups[group_name].append(line)

        if group is None:
            return groups

        return groups.get(group, [])

    def stats_by_chrom(self) -> pandas.Series:
        """
        Собирает статистику по количеству выравниваний на каждой хромосоме.
        Пропускает несопоставленные чтения (флаг 4).

        Returns:
            pandas.Series: Серия, где индекс — хромосома, а значение — количество выравниваний.
        """
        counts = {}

        for record in self.read_records():
            fields = record.split("\t")
            flag = int(fields[1])
            chrom = fields[2]

            if flag & 4:
                continue
            if chrom not in counts:
                counts[chrom] = 0
            counts[chrom] = counts[chrom] + 1

        return pandas.Series(counts, dtype="int64", name="alignments")

    def reference_length(self, cigar: str) -> int:
        """
        Вычисляет длину выравнивания по эталону на основе строки CIGAR.

        Args:
            cigar (str): Строка CIGAR (например, "10M2I5M").

        Returns:
            int: Длина участка на эталонном геноме.
        """
        length = 0
        number = ""
        for symbol in cigar:
            if symbol.isdigit():
                number = number + symbol
            else:
                if symbol in "MDN=X":
                    length = length + int(number)
                number = ""
        return length

    def intersect(self, chrom: str, start: int, end: int):
        """
        Генератор, который выдает выравнивания, пересекающиеся с заданным регионом.

        Args:
            chrom (str): Название хромосомы.
            start (int): Начальная позиция региона.
            end (int): Конечная позиция региона.

        Yields:
            str: Строка с выравниванием, которое попало в регион.
        """
        for record in self.read_records():
            fields = record.split("\t")

            flag = int(fields[1])
            record_chrom = fields[2]
            record_start = int(fields[3])
            cigar = fields[5]

            if flag & 4:
                continue
            if record_chrom != chrom or cigar == "*":
                continue

            length = self.reference_length(cigar)
            record_end = record_start + length - 1

            if record_start <= end and record_end >= start:
                yield record

if __name__ == "__main__":
    reader = SamReader("example.sam") 
    print("Количество записей:", reader.count())