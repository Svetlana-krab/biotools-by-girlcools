from abc import ABC, abstractmethod

class BioFile(ABC):
    def __init__(self, filename: str):
        self.filename = filename
        self._file = None
        
    def __enter__(self):
        # ВОТ ЗДЕСЬ открывается файл для FastaReader
        self._file = open(self.filename, "r", encoding="utf-8")
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._file:
            self._file.close()
            
    @abstractmethod
    def __iter__(self):
        pass
        
    @abstractmethod
    def count(self) -> int:
        pass

class SequenceFile(BioFile):
    @abstractmethod
    def mean_length(self) -> float:
        pass
