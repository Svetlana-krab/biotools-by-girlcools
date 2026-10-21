import pandas as pd
from .base import GenomicDataReader

class VcfReader(GenomicDataReader):
    def __init__(self, filename: str):
        super().__init__(filename)
        self._header_lines = []
        self._header_dict = {}
        self._header_columns = []
        self._is_parsed_header = False

    def __enter__(self):
        super().__enter__()
        self._parse_header()  # Парсим заголовок при открытии
        return self

    def _parse_header(self):
        """Читает шапку файла (строки, начинающиеся с #)"""
        self._header_lines = []
        self._header_dict = {}
        self._header_columns = []
        
        self._file.seek(0)  # Всегда начинаем с начала файла
        
        for line in self._file:
            if line.startswith('##'):
                self._header_lines.append(line.strip())
                if '=' in line:
                    key = line.split('=')[0].replace('##', '')
                    self._header_dict.setdefault(key, []).append(line.strip())
            elif line.startswith('#'):
                # Это строка #CHROM, запоминаем названия колонок
                self._header_columns = line.strip().replace('#', '').split('\t')
                break  # Заголовок закончился, дальше пойдут данные
                
        self._is_parsed_header = True

    def __iter__(self):
        """Генератор для чтения вариантов"""
        if not self._is_parsed_header:
            self._parse_header()
            
        self._file.seek(0)  # ВАЖНО: сбрасываем указатель в начало перед каждым чтением
        
        for line in self._file:
            if line.startswith('#') or not line.strip():
                continue  # Пропускаем заголовки и пустые строки
            
            variant = self._parse_variant_line(line)
            if variant:  # Если строка успешно распарсилась
                yield variant

    def _parse_variant_line(self, line: str) -> dict:
        """Превращает строку с данными в словарь"""
        # split() без аргументов разбивает строку по любым пробелам (табам или пробелам)
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

    # --- Методы, требуемые заданием ---

    def get_header_info(self):
        """Получение заголовка и информации по группам заголовков"""
        return self._header_dict

    def get_header_columns(self):
        """Получение названий колонок"""
        return self._header_columns

    def count(self) -> int:
        """Получение количества вариантов"""
        return sum(1 for _ in self)

    def get_chromosomes(self) -> list:
        """Возвращает список уникальных хромосом"""
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
        
        # Группируем по хромосоме и считаем количество
        stats = df.groupby('chrom').agg(count=('pos', 'count')).reset_index()
        return stats

    def filter_by_region(self, chrom: str, start: int, end: int) -> list:
        """Поиск вариантов в заданном отрезке (аналог bedtools intersect)"""
        results = []
        for variant in self:
            if variant['CHROM'] == chrom and start <= variant['POS'] <= end:
                results.append(variant)
        return results