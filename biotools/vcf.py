import pandas as pd
from .base import GenomicDataReader

class VcfReader(GenomicDataReader):
    """
    Класс для чтения и анализа файлов формата VCF (Variant Call Format).
    Наследуется от GenomicDataReader.
    """

    def __init__(self, filename: str):
        """
        Инициализирует ридер VCF.

        Args:
            filename (str): Путь к VCF-файлу.
        """
        super().__init__(filename)
        self._header_lines = []
        self._header_dict = {}
        self._header_columns = []
        self._is_parsed_header = False

    def __enter__(self):
        """
        Открывает файл и парсит заголовок при входе в контекстный менеджер (with).

        Returns:
            VcfReader: Текущий объект ридера.
        """
        super().__enter__()
        self._parse_header()
        return self

    def _parse_header(self):
        """
        Читает шапку файла (строки, начинающиеся с #).
        Сохраняет информацию из ##INFO, ##FORMAT и названия колонок из #CHROM.
        """
        self._header_lines = []
        self._header_dict = {}
        self._header_columns = []
        
        self._file.seek(0)
        
        for line in self._file:
            if line.startswith('##'):
                self._header_lines.append(line.strip())
                if '=' in line:
                    key = line.split('=')[0].replace('##', '')
                    self._header_dict.setdefault(key, []).append(line.strip())
            elif line.startswith('#'):
                self._header_columns = line.strip().replace('#', '').split('\t')
                break
                
        self._is_parsed_header = True

    def __iter__(self):
        """
        Генератор для чтения вариантов.
        Пропускает строки заголовков и пустые строки.

        Yields:
            dict: Словарь с данными одного варианта.
        """
        if not self._is_parsed_header:
            self._parse_header()
            
        self._file.seek(0)
        
        for line in self._file:
            if line.startswith('#') or not line.strip():
                continue
            
            variant = self._parse_variant_line(line)
            if variant:
                yield variant

    def _parse_variant_line(self, line: str) -> dict:
        """
        Превращает строку с данными в словарь.

        Args:
            line (str): Строка из VCF-файла.

        Returns:
            dict: Словарь с ключами CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO.
                  Возвращает None, если строка некорректна (меньше 8 колонок).
        """
        parts = line.strip().split()
        
        if len(parts) < 8:
            return None
            
        return {
            'CHROM': parts[0],
            'POS': int(parts[1]) if parts[1].isdigit() else 0,
            'ID': parts[2],
            'REF': parts[3],
            'ALT': parts[4],
            'QUAL': float(parts[5]) if parts[5] != '.' else 0.0,
            'FILTER': parts[6],
            'INFO': parts[7],
            'FORMAT': parts[8] if len(parts) > 8 else None,
            'SAMPLES': parts[9:] if len(parts) > 9 else []
        }

    def get_header_info(self) -> dict:
        """
        Возвращает информацию из заголовков файла.

        Returns:
            dict: Словарь, где ключи — это группы заголовков (INFO, FORMAT и т.д.).
        """
        return self._header_dict

    def get_header_columns(self) -> list:
        """
        Возвращает названия колонок из строки #CHROM.

        Returns:
            list: Список названий колонок.
        """
        return self._header_columns

    def count(self) -> int:
        """
        Считает общее количество вариантов в файле.

        Returns:
            int: Количество вариантов.
        """
        return sum(1 for _ in self)

    def get_chromosomes(self) -> list:
        """
        Возвращает список уникальных хромосом, встречающихся в файле.

        Returns:
            list: Список названий хромосом.
        """
        chroms = set()
        for variant in self:
            chroms.add(variant['CHROM'])
        return list(chroms)

    def get_stats_by_region(self) -> pd.DataFrame:
        """
        Строит статистику 'количество вариантов - хромосома' с помощью Pandas.

        Returns:
            pd.DataFrame: Датафрейм с колонками 'chrom' (хромосома) и 'count' (количество).
        """
        data = []
        for variant in self:
            data.append({'chrom': variant['CHROM'], 'pos': variant['POS']})
        
        df = pd.DataFrame(data)
        if df.empty:
            return df
        
        stats = df.groupby('chrom').agg(count=('pos', 'count')).reset_index()
        return stats

    def filter_by_region(self, chrom: str, start: int, end: int) -> list:
        """
        Ищет варианты в заданном геномном отрезке (аналог bedtools intersect).

        Args:
            chrom (str): Название хромосомы.
            start (int): Начальная позиция.
            end (int): Конечная позиция.

        Returns:
            list: Список словарей с вариантами, попадающими в отрезок.
        """
        results = []
        for variant in self:
            if variant['CHROM'] == chrom and start <= variant['POS'] <= end:
                results.append(variant)
        return results