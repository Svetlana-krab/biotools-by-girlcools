import pandas
from .base import GenomicDataReader
class SamReader(GenomicDataReader): 
    def __init__(self, filename): 
        super().__init__(filename)
    def __iter__(self):
        return self.read_records()
    
    def get_chromosomes(self) -> list:
        chromosomes = []
        for line in self.header("@SQ"):
            fields = line.split("\t")
            for field in fields[1:]:
                if field.startswith("SN:"):
                    chromosomes.append(field.split(":", 1)[1])
        return chromosomes

    def filter_by_region(self, chrom: str, start: int, end: int) -> list:
        return list(self.intersect(chrom, start, end))

    def read_records(self):
        with open(self.filename, "r") as file:
            for line in file:
                if line.startswith("@"):
                    continue
                if line.strip() == "":
                    continue
                yield line.rstrip("\n")

    def count(self):
        total = 0
        for record in self.read_records():
            total = total + 1
        return total
    def header(self, group=None):
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

    def stats_by_chrom(self):
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

    def reference_length(self, cigar):
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

    def intersect(self, chrom, start, end):
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

    print("Количество групп:", len(reader.header()))
    print("Строк @SQ:", len(reader.header("@SQ")))
    print("Строк @PG:", len(reader.header("@PG")))

    print("Статистика по хромосомам:")
    print(reader.stats_by_chrom())

    print("Длина по эталону:", reader.reference_length("5M2I3M2D4M"))

    overlap_count = 0

    for record in reader.intersect("1", 100, 30427671):
        overlap_count = overlap_count + 1

    print("Записей с пересечением:", overlap_count)