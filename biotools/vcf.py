import pandas as pd
from .base import GenomicDataReader

class VcfReader(GenomicDataReader):
    def __init__(self, filename: str):
        super().__init__(filename)
        self._header_lines = []
        self._header_dict = {} # Для хранения инфо из ##INFO, ##FORMAT и т.д.
        self._header_columns = []
        self._is_parsed_header = False

    def __enter__(self):
        super().__enter__()
        self._parse_header() # Парсим заголовок при открытии
        return self

    def _parse_header(self):
        """Читает шапку файла и запоминает, где начинаются данные"""
        self._header_lines = []
        self._header_dict = {}
        
        # Сохраняем текущую позицию, чтобы потом вернуться (если нужно)
        # Но так как мы используем генераторы, лучше просто запомнить строки
        while True:
            pos = self._file.tell()
            line = self._file.readline()
            if not line:
                break
            if line.startswith('##'):
                self._header_lines.append(line.strip())
                # Пример парсинга: ##INFO=<ID=DP,Number=1...>
                if '=' in line:
                    key = line.split('=')[0].replace('##', '')
                    self._header_dict.setdefault(key, []).append(line.strip())
            elif line.startswith('#'):
                self._header_columns = line.strip().replace('#', '').split('\t')
                # Файл начался, возвращаем указатель назад, чтобы итератор сработал
                self._file.seek(pos)
                break
        self._is_parsed_header = True

    def __iter__(self):
        """Генератор для чтения данных (вариантов)"""
        if not self._is_parsed_header:
            self._parse_header()
        
        # Пропускаем заголовок (он уже прочитан в _parse_header, 
        # но т.к. мы сделали seek, итератор начнет с заголовка. Пропускаем его)
        for line in self._file:
            if line.startswith('#'):
                continue
            yield self._parse_variant_line(line)

    def _parse_variant_line(self, line: str) -> dict:
        """Превращает строку в словарь для удобства"""
        parts = line.strip().split('\t')
        if len(parts) < 8:
            return None
        
        return {
            'CHROM': parts[0],
            'POS': int(parts[1]),
            'ID': parts[2],
            'REF': parts[3],
            'ALT': parts[4],
            'QUAL': float(parts[5]) if parts[5] != '.' else None,
            'FILTER': parts[6],
            'INFO': parts[7],
            'FORMAT': parts[8] if len(parts) > 8 else None,
            'SAMPLES': parts[9:] if len(parts) > 9 else []
        }

    # --- Методы, требуемые заданием (скриншот 3) ---

    def get_header_info(self):
        """Получение заголовка и информации по группам заголовков"""
        return self._header_dict

    def get_header_columns(self):
        return self._header_columns

    def count(self) -> int:
        """Получение количества вариантов"""
        return sum(1 for _ in self)

    def get_chromosomes(self) -> list:
        """Возвращает уникальные хромосомы"""
        chroms = set()
        for variant in self:
            chroms.add(variant['CHROM'])
        return list(chroms)

    def get_stats_by_region(self) -> pd.DataFrame:
        """Статистика 'количество вариантов - регион' (через pandas)"""
        data = []
        for variant in self:
            data.append({'chrom': variant['CHROM'], 'pos': variant['POS']})
        
        df = pd.DataFrame(data)
        if df.empty:
            return df
        
        # Группируем по хромосоме
        stats = df.groupby('chrom').agg(count=('pos', 'count')).reset_index()
        return stats

    def filter_by_region(self, chrom: str, start: int, end: int) -> list:
        """Аналог bedtools intersect: поиск вариантов в отрезке"""
        results = []
        for variant in self:
            if variant['CHROM'] == chrom and start <= variant['POS'] <= end:
                results.append(variant)
        return results