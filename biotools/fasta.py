from biotools.base import SequenceFile

class Seq:
    NUCLEOTIDES = set("ACGTUN")
    PROTEIN = set("ACDEFGHIKLMNPQRSTVWY")

    def __init__(self, sequence, header=""):
        if not sequence:
            raise ValueError("Последовательность не может быть пустой")
        self._sequence = sequence.upper().strip()
        self._header = header.strip()

    @property
    def sequence(self):
        return self._sequence

    @property
    def header(self):
        return self._header

    def __len__(self):
        return len(self._sequence)

    def get_alphabet(self):
        symbols = set(self._sequence)
        if symbols.issubset(self.NUCLEOTIDES):
            return "nucleotide"
        if symbols.issubset(self.PROTEIN):
            return "protein"
        return "unknown"

    def gc_content(self):
        if self.get_alphabet() != "nucleotide":
            raise ValueError("Только для нуклеотидов")
        g = self._sequence.count("G")
        c = self._sequence.count("C")
        return (g + c) / len(self) * 100

    def __repr__(self):
        return f"Seq(header='{self._header}', length={len(self)})"


class FastaReader(SequenceFile):
    def __iter__(self):
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
        return sum(1 for _ in self)

    def mean_length(self) -> float:
        lengths = [len(rec) for rec in self] 
        return sum(lengths) / len(lengths) if lengths else 0.0