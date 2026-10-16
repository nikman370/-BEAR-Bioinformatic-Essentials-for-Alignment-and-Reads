import os

from helpers.bio_files import read_multiline_fasta
from helpers.fastq import check_gc, check_length, check_quality


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Фильтрует FASTQ-файл по GC-составу, длине ридов и среднему качеству Phred33.

    Прошедшие все фильтры последовательности сохраняются в отдельный файл
    внутри создаваемой папки 'filtered/'.

    Аргументы:
        input_fastq (str):
            Путь к исходному входному FASTQ-файлу.
        output_fastq (str):
            Имя выходного файла, в который запишутся отфильтрованные риды.
        gc_bounds (tuple[float, float] | float, optional):
            Интервал GC-состава в процентах. Если передан кортеж (min, max),
            включаются обе границы. Если передано одно число, оно считается
            верхней границей (от 0 до числа). По умолчанию (0, 100).
        length_bounds (tuple[int, int] | int, optional):
            Интервал длины рида в нуклеотидах. Принцип аналогичен gc_bounds.
            По умолчанию (0, 2**32).
        quality_threshold (float, optional):
            Минимальный порог среднего качества рида по шкале Phred33.
            Риды со средним значением ниже порога отбрасываются. По умолчанию 0.

    Возвращает:
        None
    """

    # Создаем папку "filtered", если её еще нет
    os.makedirs("filtered", exist_ok=True)

    # путь к выходному файлу внутри папки filtered
    output_path = os.path.join("filtered", output_fastq)

    # Открываем входной файл для чтения, а выходной — для записи
    with (
        open(input_fastq, "r", encoding="utf-8") as infile,
        open(output_path, "w", encoding="utf-8") as outfile,
    ):
        while True:
            # Пытаемся прочитать первую строку (ID рида)
            line1 = infile.readline()

            # Если строка пустая, значит мы дошли до конца файла — выходим из цикла
            if not line1:
                break

            # Дочитываем оставшиеся 3 строки для текущего рида
            line2 = infile.readline()  # Последовательность
            line3 = infile.readline()  # Плюс (+)
            line4 = infile.readline()  # Строка качества

            # Метод .strip() убирает невидимый символ переноса строки '\n' в конце, чтобы он не мешал нам считать длину и GC-состав
            seq = line2.strip()
            qual = line4.strip()

            # Передаем данные в наши функции-проверяторы
            is_gc_ok = check_gc(seq, gc_bounds)
            is_len_ok = check_length(seq, length_bounds)
            is_qual_ok = check_quality(qual, quality_threshold)

            # Если рид прошел абсолютно ВСЕ проверки, записываем его в выходной файл
            if is_gc_ok and is_len_ok and is_qual_ok:
                outfile.write(line1)
                outfile.write(line2)
                outfile.write(line3)
                outfile.write(line4)


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str | None = None
) -> None:
    """
    Конвертирует FASTA-файл с многострочными последовательностями в однострочный.

    Arguments:
        input_fasta (str):
            Путь к входному FASTA-файлу.
        output_fasta (str | None, optional):
            Путь к выходному файлу. Если не указан, файл сохраняется
            с суффиксом '_oneline.fasta' в той же директории. По умолчанию None.

    """
    if output_fasta is None:
        base, ext = os.path.splitext(input_fasta)
        output_fasta = base + "_oneline" + ext

    fasta_data = read_multiline_fasta(input_fasta)

    with open(output_fasta, "w", encoding="utf-8") as f:
        for header, sequence in fasta_data:
            f.write(header + "\n")
            f.write(sequence + "\n")


"""
# Для тестов и ручного запуска
if __name__ == "__main__":
    print("Запуск фильтрации FASTQ...")

    # Вызываем функцию с вашим тестовым файлом
    filter_fastq(
        input_fastq="example_fastq.fastq",   # тут вместо "example_fastq.fastq" - напишите имя СВОЕГО файла, который должен лежать в папке вместе с bear.py
        output_fastq="output.fastq",
        gc_bounds=(50, 70),  # Оставим риды с GC от 30% до 70% - тут вы можете поменять 
        length_bounds=(50, 150),  # Длина от 50 до 150 нуклеотидов - тут вы можете поменять 
        quality_threshold=20,  # Среднее качество не ниже Q20 - тут вы можете поменять: Q10 (Качество 10) — Вероятность ошибки 10% ; Q20 (Качество 20) — Вероятность ошибки 1% ; Q30 (Качество 30) — Вероятность ошибки 0.1% ; Q40 (Качество 40) — Вероятность ошибки 0.01%.
    )

    print("Фильтрация успешно завершена! Проверьте папку 'filtered'.")
"""
